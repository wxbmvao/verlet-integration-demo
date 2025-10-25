const canvas = document.getElementById('canvas');
const ctx = canvas.getContext('2d');

// scale, 1 pixel = 0.01 meters (1cm per pixel)
const PIXELS_PER_METER = 100;
const METER_TO_PIXEL = PIXELS_PER_METER;
const PIXEL_TO_METER = 1 / PIXELS_PER_METER;

// material properties database (real values)
const MATERIALS = {
  steel: {
    density: 7850, // kg/m^3
    restitution: 0.4, // coefficient of restitution
    staticFriction: 0.74,
    dynamicFriction: 0.57,
    dragCoef: 0.47, // for circle
  },
  wood: {
    density: 600,
    restitution: 0.55,
    staticFriction: 0.45,
    dynamicFriction: 0.35,
    dragCoef: 1.05,
  },
  rubber: {
    density: 1100,
    restitution: 0.85,
    staticFriction: 1.15,
    dynamicFriction: 1.05,
    dragCoef: 0.47,
  },
  ice: {
    density: 917,
    restitution: 0.1,
    staticFriction: 0.05,
    dynamicFriction: 0.03,
    dragCoef: 0.47,
  },
  concrete: {
    density: 2400,
    restitution: 0.3,
    staticFriction: 0.62,
    dynamicFriction: 0.52,
    dragCoef: 1.15,
  },
};

// 2d vector math
class Vec2 {
  constructor(x = 0, y = 0) {
    this.x = x;
    this.y = y;
  }

  add(v) {
    return new Vec2(this.x + v.x, this.y + v.y);
  }
  sub(v) {
    return new Vec2(this.x - v.x, this.y - v.y);
  }
  mul(s) {
    return new Vec2(this.x * s, this.y * s);
  }
  dot(v) {
    return this.x * v.x + this.y * v.y;
  }
  cross(v) {
    return this.x * v.y - this.y * v.x;
  }
  length() {
    return Math.sqrt(this.x * this.x + this.y * this.y);
  }
  lengthSq() {
    return this.x * this.x + this.y * this.y;
  }

  normalize() {
    const len = this.length();
    return len > 0 ? this.mul(1 / len) : new Vec2(0, 0);
  }

  rotate(angle) {
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    return new Vec2(this.x * c - this.y * s, this.x * s + this.y * c);
  }

  perp() {
    return new Vec2(-this.y, this.x);
  }
}

// rigid body
class RigidBody {
  constructor(
    type,
    positionPixels,
    sizeMeters,
    density,
    restitution,
    staticFriction,
    dynamicFriction
  ) {
    this.type = type;

    // position in meters for physics calculations
    this.position = new Vec2(
      positionPixels.x * PIXEL_TO_METER,
      positionPixels.y * PIXEL_TO_METER
    );

    this.velocity = new Vec2(0, 0); // m/s
    this.force = new Vec2(0, 0); // newtons
    this.angle = 0; // radians
    this.angularVelocity = 0; // rad/s
    this.torque = 0; // N x m

    this.sizeMeters = sizeMeters;
    this.density = density; // kg/m^3
    this.restitution = restitution;
    this.staticFriction = staticFriction;
    this.dynamicFriction = dynamicFriction;

    // calculate mass and moment of inertia using real formulas
    if (type === 'circle') {
      // circle, m = ρ x (4/3)πr^3 (but we're doing 2d so πr^2)
      const area = Math.PI * sizeMeters * sizeMeters;
      this.mass = area * density;
      // moment of inertia for solid disk: I = (1/2)mr^2
      this.inertia = 0.5 * this.mass * sizeMeters * sizeMeters;
      this.radius = sizeMeters;
      this.dragCoef = 0.47; // sphere drag coefficient
    } else if (type === 'box') {
      // box, m = ρ x (2r)^2
      const side = sizeMeters * 2;
      this.mass = side * side * density;
      // moment of inertia for square: I = (1/6)m(a^2)
      this.inertia = (this.mass * side * side) / 6;
      this.dragCoef = 1.05; // box drag coefficient
    } else if (type === 'polygon') {
      // approximate as circle for mass
      const area = 2.5 * sizeMeters * sizeMeters;
      this.mass = area * density;
      this.inertia = 0.4 * this.mass * sizeMeters * sizeMeters;
      this.dragCoef = 1.0;
    }

    this.invMass = this.mass > 0 ? 1 / this.mass : 0;
    this.invInertia = this.inertia > 0 ? 1 / this.inertia : 0;

    this.color = `hsl(${Math.random() * 360}, 45%, 50%)`;

    // set up vertices for polygons (in meters)
    if (type === 'box') {
      const s = sizeMeters;
      this.localVertices = [
        new Vec2(-s, -s),
        new Vec2(s, -s),
        new Vec2(s, s),
        new Vec2(-s, s),
      ];
    } else if (type === 'polygon') {
      const sides = 5;
      this.localVertices = [];
      for (let i = 0; i < sides; i++) {
        const ang = (Math.PI * 2 * i) / sides - Math.PI / 2;
        this.localVertices.push(
          new Vec2(Math.cos(ang) * sizeMeters, Math.sin(ang) * sizeMeters)
        );
      }
    }
  }

  getVertices() {
    if (!this.localVertices) return [];
    return this.localVertices.map((v) =>
      v.rotate(this.angle).add(this.position)
    );
  }

  getAABB() {
    if (this.type === 'circle') {
      return {
        minX: this.position.x - this.radius,
        maxX: this.position.x + this.radius,
        minY: this.position.y - this.radius,
        maxY: this.position.y + this.radius,
      };
    } else {
      const verts = this.getVertices();
      return {
        minX: Math.min(...verts.map((v) => v.x)),
        maxX: Math.max(...verts.map((v) => v.x)),
        minY: Math.min(...verts.map((v) => v.y)),
        maxY: Math.max(...verts.map((v) => v.y)),
      };
    }
  }

  applyForce(force, point) {
    this.force = this.force.add(force);
    if (point) {
      const r = point.sub(this.position);
      this.torque += r.cross(force);
    }
  }

  integrate(dt) {
    if (this.invMass === 0) return;

    // f = ma, so a = f/m
    const accel = this.force.mul(this.invMass);

    // realistic air drag, F_d = 0.5 x ρ x v^2 x C_d x A
    const airDensity = parseFloat(document.getElementById('airDensity').value);
    const speed = this.velocity.length();

    if (speed > 0.01) {
      // cross sectional area
      let area;
      if (this.type === 'circle') {
        area = Math.PI * this.radius * this.radius;
      } else {
        area = 4 * this.sizeMeters * this.sizeMeters;
      }

      // drag force magnitude
      const dragMag = 0.5 * airDensity * speed * speed * this.dragCoef * area;
      const dragForce = this.velocity.normalize().mul(-dragMag);

      // apply drag acceleration
      const dragAccel = dragForce.mul(this.invMass);
      this.velocity = this.velocity.add(dragAccel.mul(dt));
    }

    // semi implicit euler integration (more stable than explicit euler)
    this.velocity = this.velocity.add(accel.mul(dt));
    this.position = this.position.add(this.velocity.mul(dt));

    // rotational dynamics
    if (document.getElementById('rotationEnabled').checked) {
      const angularAccel = this.torque * this.invInertia;
      this.angularVelocity += angularAccel * dt;

      // realistic angular damping from air resistance
      const angularDrag =
        0.02 *
        this.angularVelocity *
        this.angularVelocity *
        Math.sign(this.angularVelocity);
      this.angularVelocity -= angularDrag * dt;

      this.angle += this.angularVelocity * dt;
    }

    this.force = new Vec2(0, 0);
    this.torque = 0;
  }

  draw() {
    // convert meters to pixels for rendering
    const px = this.position.x * METER_TO_PIXEL;
    const py = this.position.y * METER_TO_PIXEL;

    ctx.save();
    ctx.translate(px, py);
    ctx.rotate(this.angle);

    ctx.fillStyle = this.color;
    ctx.strokeStyle = '#3d4d49';
    ctx.lineWidth = 2;

    if (this.type === 'circle') {
      const radiusPx = this.radius * METER_TO_PIXEL;
      ctx.beginPath();
      ctx.arc(0, 0, radiusPx, 0, Math.PI * 2);
      ctx.fill();
      ctx.stroke();

      // small dot to show rotation
      ctx.fillStyle = '#2d3e3a';
      ctx.beginPath();
      ctx.arc(radiusPx * 0.5, 0, 3, 0, Math.PI * 2);
      ctx.fill();
    } else {
      const verts = this.localVertices.map((v) => v.mul(METER_TO_PIXEL));
      ctx.beginPath();
      ctx.moveTo(verts[0].x, verts[0].y);
      for (let i = 1; i < verts.length; i++) {
        ctx.lineTo(verts[i].x, verts[i].y);
      }
      ctx.closePath();
      ctx.fill();
      ctx.stroke();
    }

    ctx.restore();

    // velocity vector
    if (
      document.getElementById('showVectors').checked &&
      this.velocity.length() > 0.1
    ) {
      const vel = this.velocity.mul(METER_TO_PIXEL * 0.1);
      ctx.beginPath();
      ctx.moveTo(px, py);
      ctx.lineTo(px + vel.x, py + vel.y);
      ctx.strokeStyle = '#5a7a6f';
      ctx.lineWidth = 2;
      ctx.stroke();

      const angle = Math.atan2(vel.y, vel.x);
      const headlen = 8;
      ctx.beginPath();
      ctx.moveTo(px + vel.x, py + vel.y);
      ctx.lineTo(
        px + vel.x - headlen * Math.cos(angle - Math.PI / 6),
        py + vel.y - headlen * Math.sin(angle - Math.PI / 6)
      );
      ctx.moveTo(px + vel.x, py + vel.y);
      ctx.lineTo(
        px + vel.x - headlen * Math.cos(angle + Math.PI / 6),
        py + vel.y - headlen * Math.sin(angle + Math.PI / 6)
      );
      ctx.stroke();
    }
  }
}

// collision detection functions
function circleVsCircle(a, b) {
  const diff = b.position.sub(a.position);
  const dist = diff.length();
  const radiusSum = a.radius + b.radius;

  if (dist < radiusSum && dist > 0.0001) {
    const normal = diff.normalize();
    const penetration = radiusSum - dist;
    const contactPoint = a.position.add(normal.mul(a.radius));
    return { normal, penetration, point: contactPoint };
  }
  return null;
}

function polygonVsPolygon(a, b) {
  const vertsA = a.getVertices();
  const vertsB = b.getVertices();

  let minOverlap = Infinity;
  let smallestAxis = null;

  function projectPolygon(vertices, axis) {
    let min = Infinity,
      max = -Infinity;
    for (const v of vertices) {
      const proj = v.dot(axis);
      min = Math.min(min, proj);
      max = Math.max(max, proj);
    }
    return { min, max };
  }

  function checkAxes(verts1, verts2) {
    for (let i = 0; i < verts1.length; i++) {
      const v1 = verts1[i];
      const v2 = verts1[(i + 1) % verts1.length];
      const edge = v2.sub(v1);
      const axis = edge.perp().normalize();

      const projA = projectPolygon(verts1, axis);
      const projB = projectPolygon(verts2, axis);

      if (projA.max < projB.min || projB.max < projA.min) {
        return false;
      }

      const overlap = Math.min(projA.max - projB.min, projB.max - projA.min);
      if (overlap < minOverlap) {
        minOverlap = overlap;
        smallestAxis = axis;
      }
    }
    return true;
  }

  if (!checkAxes(vertsA, vertsB) || !checkAxes(vertsB, vertsA)) {
    return null;
  }

  const diff = b.position.sub(a.position);
  if (diff.dot(smallestAxis) < 0) {
    smallestAxis = smallestAxis.mul(-1);
  }

  const contactPoint = a.position.add(smallestAxis.mul(minOverlap / 2));
  return { normal: smallestAxis, penetration: minOverlap, point: contactPoint };
}

function circleVsPolygon(circle, polygon) {
  const vertices = polygon.getVertices();
  let minDist = Infinity;
  let closestPoint = null;

  for (let i = 0; i < vertices.length; i++) {
    const v1 = vertices[i];
    const v2 = vertices[(i + 1) % vertices.length];
    const edge = v2.sub(v1);
    const toCircle = circle.position.sub(v1);

    const edgeLenSq = edge.lengthSq();
    const t =
      edgeLenSq > 0
        ? Math.max(0, Math.min(1, toCircle.dot(edge) / edgeLenSq))
        : 0;
    const closest = v1.add(edge.mul(t));
    const dist = circle.position.sub(closest).length();

    if (dist < minDist) {
      minDist = dist;
      closestPoint = closest;
    }
  }

  if (minDist < circle.radius) {
    const normal = circle.position.sub(closestPoint).normalize();
    const penetration = circle.radius - minDist;
    return { normal, penetration, point: closestPoint };
  }

  return null;
}

function detectCollision(a, b) {
  if (a.type === 'circle' && b.type === 'circle') {
    return circleVsCircle(a, b);
  } else if (a.type === 'circle') {
    return circleVsPolygon(a, b);
  } else if (b.type === 'circle') {
    const contact = circleVsPolygon(b, a);
    if (contact) contact.normal = contact.normal.mul(-1);
    return contact;
  } else {
    return polygonVsPolygon(a, b);
  }
}

// realistic impulse based collision resolution
function resolveCollision(a, b, contact) {
  if (!contact) return;

  const { normal, penetration, point } = contact;

  // position correction using baumgarte stabilization
  const slop = 0.0001; // very small to prevent jitter
  const percent = 0.8; // correction percentage
  const correction =
    (Math.max(penetration - slop, 0) / (a.invMass + b.invMass)) * percent;
  const correctionVec = normal.mul(correction);

  a.position = a.position.sub(correctionVec.mul(a.invMass));
  b.position = b.position.add(correctionVec.mul(b.invMass));

  // calculate relative velocity at contact point
  const ra = point.sub(a.position);
  const rb = point.sub(b.position);

  const va = a.velocity.add(
    new Vec2(-ra.y * a.angularVelocity, ra.x * a.angularVelocity)
  );
  const vb = b.velocity.add(
    new Vec2(-rb.y * b.angularVelocity, rb.x * b.angularVelocity)
  );
  const relVel = vb.sub(va);

  const velAlongNormal = relVel.dot(normal);
  if (velAlongNormal > 0) return;

  // coefficient of restitution (combined using minimum)
  const e = Math.min(a.restitution, b.restitution);

  // calculate impulse scalar with rotation
  const raCrossN = ra.cross(normal);
  const rbCrossN = rb.cross(normal);

  const invMassSum =
    a.invMass +
    b.invMass +
    raCrossN * raCrossN * a.invInertia +
    rbCrossN * rbCrossN * b.invInertia;

  let j = (-(1 + e) * velAlongNormal) / invMassSum;
  const impulse = normal.mul(j);

  // apply normal impulse
  a.velocity = a.velocity.sub(impulse.mul(a.invMass));
  b.velocity = b.velocity.add(impulse.mul(b.invMass));
  a.angularVelocity -= raCrossN * j * a.invInertia;
  b.angularVelocity += rbCrossN * j * b.invInertia;

  // coulomb friction model
  const tangent = relVel.sub(normal.mul(velAlongNormal));
  if (tangent.lengthSq() > 0.0001) {
    const t = tangent.normalize();
    const velAlongTangent = relVel.dot(t);

    const raCrossT = ra.cross(t);
    const rbCrossT = rb.cross(t);
    const invMassSumT =
      a.invMass +
      b.invMass +
      raCrossT * raCrossT * a.invInertia +
      rbCrossT * rbCrossT * b.invInertia;

    let jt = -velAlongTangent / invMassSumT;

    // coulomb friction, use static or dynamic friction
    const muS = Math.sqrt(a.staticFriction * b.staticFriction);
    const muK = Math.sqrt(a.dynamicFriction * b.dynamicFriction);

    let frictionImpulse;
    if (Math.abs(jt) < j * muS) {
      // static friction
      frictionImpulse = t.mul(jt);
    } else {
      // kinetic friction
      frictionImpulse = t.mul(-j * muK);
    }

    a.velocity = a.velocity.sub(frictionImpulse.mul(a.invMass));
    b.velocity = b.velocity.add(frictionImpulse.mul(b.invMass));
    a.angularVelocity -= ra.cross(frictionImpulse) * a.invInertia;
    b.angularVelocity += rb.cross(frictionImpulse) * b.invInertia;
  }

  return point;
}

// main app
const app = {
  bodies: [],
  paused: false,
  lastTime: performance.now(),
  frameCount: 0,
  fps: 60,
  collisionCount: 0,
  contactPoints: [],
  accumulator: 0,
  fixedDt: 1 / 120, // 120hz physics for accuracy
  substeps: 10,

  updateMaterial() {
    const material = document.getElementById('material').value;
    if (material === 'custom') return;

    const props = MATERIALS[material];
    document.getElementById('density').value = props.density;
    document.getElementById('restitution').value = props.restitution;
    document.getElementById('staticFriction').value = props.staticFriction;
    document.getElementById('dynamicFriction').value = props.dynamicFriction;
  },

  addObject(x, y, vx = 0, vy = 0) {
    const type = document.getElementById('shapeType').value;
    const sizeMeters = parseFloat(document.getElementById('size').value);
    const density = parseFloat(document.getElementById('density').value);
    const restitution = parseFloat(
      document.getElementById('restitution').value
    );
    const staticFriction = parseFloat(
      document.getElementById('staticFriction').value
    );
    const dynamicFriction = parseFloat(
      document.getElementById('dynamicFriction').value
    );

    if (x === undefined) {
      x = Math.random() * (canvas.width - 200) + 100;
      y = 100;
    }

    const body = new RigidBody(
      type,
      new Vec2(x, y),
      sizeMeters,
      density,
      restitution,
      staticFriction,
      dynamicFriction
    );
    body.velocity = new Vec2(vx * PIXEL_TO_METER, vy * PIXEL_TO_METER);
    body.angularVelocity = (Math.random() - 0.5) * 2;
    this.bodies.push(body);
  },

  addMultiple() {
    for (let i = 0; i < 5; i++) {
      setTimeout(() => this.addObject(), i * 100);
    }
  },

  createStack() {
    const sizeMeters = parseFloat(document.getElementById('size').value);
    const type = document.getElementById('shapeType').value;
    const density = parseFloat(document.getElementById('density').value);
    const restitution = parseFloat(
      document.getElementById('restitution').value
    );
    const staticFriction = parseFloat(
      document.getElementById('staticFriction').value
    );
    const dynamicFriction = parseFloat(
      document.getElementById('dynamicFriction').value
    );

    const x = canvas.width / 2;
    const layers = 6;
    const spacing = sizeMeters * 2.1 * METER_TO_PIXEL;

    for (let i = 0; i < layers; i++) {
      const y = canvas.height - 50 - i * spacing;
      const body = new RigidBody(
        type,
        new Vec2(x, y),
        sizeMeters,
        density,
        restitution,
        staticFriction,
        dynamicFriction
      );
      body.angle = 0;
      body.angularVelocity = 0;
      this.bodies.push(body);
    }
  },

  createPyramid() {
    const sizeMeters = parseFloat(document.getElementById('size').value);
    const type = 'box';
    const density = parseFloat(document.getElementById('density').value);
    const restitution = parseFloat(
      document.getElementById('restitution').value
    );
    const staticFriction = parseFloat(
      document.getElementById('staticFriction').value
    );
    const dynamicFriction = parseFloat(
      document.getElementById('dynamicFriction').value
    );

    const startX = canvas.width / 2;
    const layers = 5;
    const spacing = sizeMeters * 2.1 * METER_TO_PIXEL;

    for (let layer = 0; layer < layers; layer++) {
      const y = canvas.height - 50 - layer * spacing;
      const boxesInLayer = layers - layer;
      const startXForLayer = startX - ((boxesInLayer - 1) * spacing) / 2;

      for (let i = 0; i < boxesInLayer; i++) {
        const x = startXForLayer + i * spacing;
        const body = new RigidBody(
          type,
          new Vec2(x, y),
          sizeMeters,
          density,
          restitution,
          staticFriction,
          dynamicFriction
        );
        body.angle = 0;
        body.angularVelocity = 0;
        this.bodies.push(body);
      }
    }
  },

  shootObject() {
    const sizeMeters = parseFloat(document.getElementById('size').value);
    const type = document.getElementById('shapeType').value;
    const density = parseFloat(document.getElementById('density').value);
    const restitution = parseFloat(
      document.getElementById('restitution').value
    );
    const staticFriction = parseFloat(
      document.getElementById('staticFriction').value
    );
    const dynamicFriction = parseFloat(
      document.getElementById('dynamicFriction').value
    );

    const body = new RigidBody(
      type,
      new Vec2(50, canvas.height / 2),
      sizeMeters,
      density,
      restitution,
      staticFriction,
      dynamicFriction
    );
    body.velocity = new Vec2(15, (Math.random() - 0.5) * 3);
    body.angularVelocity = (Math.random() - 0.5) * 10;
    this.bodies.push(body);
  },

  togglePause() {
    this.paused = !this.paused;
  },

  clearAll() {
    this.bodies = [];
  },

  physicsStep(dt) {
    // gravity in m/s^2
    const gravityAccel = parseFloat(document.getElementById('gravity').value);
    const gravity = new Vec2(0, gravityAccel);

    // apply gravity force
    this.bodies.forEach((body) => {
      const gravityForce = gravity.mul(body.mass);
      body.applyForce(gravityForce);
    });

    // integrate motion
    this.bodies.forEach((body) => body.integrate(dt));

    // boundary collisions (canvas edges)
    const canvasWidthMeters = canvas.width * PIXEL_TO_METER;
    const canvasHeightMeters = canvas.height * PIXEL_TO_METER;

    this.bodies.forEach((body) => {
      const aabb = body.getAABB();

      // ground
      if (aabb.maxY > canvasHeightMeters) {
        body.position.y -= aabb.maxY - canvasHeightMeters;
        if (body.velocity.y > 0) {
          body.velocity.y *= -body.restitution;
          body.velocity.x *= 1 - body.dynamicFriction * 0.5;
          body.angularVelocity *= 1 - body.dynamicFriction * 0.5;
        }
      }

      // ceiling
      if (aabb.minY < 0) {
        body.position.y -= aabb.minY;
        if (body.velocity.y < 0) {
          body.velocity.y *= -body.restitution;
        }
      }

      // right wall
      if (aabb.maxX > canvasWidthMeters) {
        body.position.x -= aabb.maxX - canvasWidthMeters;
        if (body.velocity.x > 0) {
          body.velocity.x *= -body.restitution;
          body.angularVelocity *= 1 - body.dynamicFriction * 0.5;
        }
      }

      // left wall
      if (aabb.minX < 0) {
        body.position.x -= aabb.minX;
        if (body.velocity.x < 0) {
          body.velocity.x *= -body.restitution;
          body.angularVelocity *= 1 - body.dynamicFriction * 0.5;
        }
      }
    });

    // collision detection and resolution
    this.collisionCount = 0;
    this.contactPoints = [];

    if (document.getElementById('collisionsEnabled').checked) {
      // iterations for constraint solving
      for (let iter = 0; iter < this.substeps; iter++) {
        for (let i = 0; i < this.bodies.length; i++) {
          for (let j = i + 1; j < this.bodies.length; j++) {
            const contact = detectCollision(this.bodies[i], this.bodies[j]);
            if (contact) {
              const point = resolveCollision(
                this.bodies[i],
                this.bodies[j],
                contact
              );
              if (point && iter === 0) {
                this.contactPoints.push(point);
                this.collisionCount++;
              }
            }
          }
        }
      }
    }
  },

  update(currentTime) {
    const deltaTime = (currentTime - this.lastTime) / 1000;
    this.lastTime = currentTime;

    // frames
    this.frameCount++;
    if (this.frameCount % 30 === 0) {
      this.fps = Math.round(1 / deltaTime);
    }

    // fixed timestep physics
    if (!this.paused) {
      this.accumulator += Math.min(deltaTime, 0.1);

      while (this.accumulator >= this.fixedDt) {
        this.physicsStep(this.fixedDt);
        this.accumulator -= this.fixedDt;
      }
    }

    // render
    ctx.fillStyle = '#f5f9f7';
    ctx.fillRect(0, 0, canvas.width, canvas.height);

    this.bodies.forEach((body) => body.draw());

    // show contact points
    if (document.getElementById('showContacts').checked) {
      this.contactPoints.forEach((point) => {
        const px = point.x * METER_TO_PIXEL;
        const py = point.y * METER_TO_PIXEL;
        ctx.fillStyle = '#ff6b6b';
        ctx.beginPath();
        ctx.arc(px, py, 4, 0, Math.PI * 2);
        ctx.fill();
      });
    }

    // update stats
    document.getElementById('objectCount').textContent = this.bodies.length;
    document.getElementById('fps').textContent = this.fps;
    document.getElementById('substeps').textContent = this.substeps;
    document.getElementById('collisions').textContent = this.collisionCount;

    // calculate total energy
    const gravityAccel = parseFloat(document.getElementById('gravity').value);
    const canvasHeightMeters = canvas.height * PIXEL_TO_METER;

    let totalEnergy = 0;
    this.bodies.forEach((body) => {
      // kinetic energy: KE = 0.5 * m * v²
      const ke = 0.5 * body.mass * body.velocity.lengthSq();
      // rotational energy: RE = 0.5 * I * ω²
      const re =
        0.5 * body.inertia * body.angularVelocity * body.angularVelocity;
      // potential energy: PE = m * g * h
      const pe =
        body.mass * gravityAccel * (canvasHeightMeters - body.position.y);
      totalEnergy += ke + re + pe;
    });
    document.getElementById('energy').textContent = totalEnergy.toFixed(2);

    requestAnimationFrame((time) => this.update(time));
  },
};

// mouse interaction
let mouseDown = false;
let mousePos = new Vec2(0, 0);

canvas.addEventListener('mousedown', (e) => {
  mouseDown = true;
  const rect = canvas.getBoundingClientRect();
  mousePos = new Vec2(e.clientX - rect.left, e.clientY - rect.top);
});

canvas.addEventListener('mousemove', (e) => {
  const rect = canvas.getBoundingClientRect();
  mousePos = new Vec2(e.clientX - rect.left, e.clientY - rect.top);
});

canvas.addEventListener('mouseup', (e) => {
  if (mouseDown) {
    const rect = canvas.getBoundingClientRect();
    const endPos = new Vec2(e.clientX - rect.left, e.clientY - rect.top);
    const velocity = endPos.sub(mousePos);

    // only add if there was actual dragging
    const dragDistance = velocity.length();
    if (dragDistance > 5) {
      app.addObject(mousePos.x, mousePos.y, velocity.x * 0.5, velocity.y * 0.5);
    } else {
      // just a click, not a drag
      app.addObject(mousePos.x, mousePos.y);
    }
  }
  mouseDown = false;
});

// init
app.update(performance.now());
