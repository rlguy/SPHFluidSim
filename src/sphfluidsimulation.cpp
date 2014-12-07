#include "sphfluidsimulation.h"

SPHFluidSimulation::SPHFluidSimulation()
{
    h = 1.0;
}

SPHFluidSimulation::SPHFluidSimulation(double smoothingRadius)
{
    h = smoothingRadius;
    grid = SpatialGrid(h);
    gravityForce = glm::vec3(0.0, -gravityMagnitude, 0.0);
}

void SPHFluidSimulation::setBounds(double _xmin, double _xmax,
                                   double _ymin, double _ymax,
                                   double _zmin, double _zmax) {
    xmin = _xmin; xmax = _xmax;
    ymin = _ymin; ymax = _ymax;
    zmin = _zmin; zmax = _zmax;
}

void SPHFluidSimulation::setDampingConstant(double c) {
    motionDampingCoefficient = c;
}

void SPHFluidSimulation::addFluidParticles(std::vector<glm::vec3> points) {
    for (uint i=0; i<points.size(); i++) {
        addFluidParticle(points[i], glm::vec3(0.0, 0.0, 0.0));
    }
}

int SPHFluidSimulation::addObstacleParticles(std::vector<glm::vec3> points) {
    SPHObstacle *obs = new SPHObstacle();
    obs->id = getUniqueObstacleID();

    SPHParticle *p;
    for (uint i=0; i<points.size(); i++) {
        p = addObstacleParticle(points[i]);
        obs->particles.push_back(p);
    }
    obstacles.push_back(obs);

    std::pair<int,SPHObstacle*> pair(obs->id, obs);
    obstaclesByID.insert(pair);

    return obs->id;
}

int SPHFluidSimulation::getUniqueObstacleID() {
    int id = currentObstacleID;
    currentObstacleID++;
    return id;
}

SPHParticle* SPHFluidSimulation::createSPHParticle(glm::vec3 pos, glm::vec3 velocity) {
    SPHParticle *s = new SPHParticle();

    s->position = pos;
    s->velocity = velocity;
    s->velocityAtHalfTimeStep = glm::vec3(0.0, 0.0, 0.0);
    s->isHalfTimeStepVelocityInitialized = false;
    s->XSPHVelocity = glm::vec3(0.0, 0.0, 0.0);
    s->soundSpeed = 0;
    s->acceleration = glm::vec3(0.0, 0.0, 0.0);

    // Create pressure offset from initial pressure
    // uniform densities will cause uniform pressure of 0, meaning no acceleration
    // of system
    s->density = initialDensity;
    s->densityVelocity = 0.0;

    // mass of sphere
    s->mass = particleMass;

    // initial pressure will be calculated once all particles are in place
    s->pressure = 0.0;

    return s;
}

SPHParticle* SPHFluidSimulation::createSPHObstacleParticle(glm::vec3 pos) {
    SPHParticle *p = createSPHParticle(pos, glm::vec3(0.0, 0.0, 0.0));
    return p;
}

void SPHFluidSimulation::addFluidParticle(glm::vec3 pos, glm::vec3 velocity) {
    SPHParticle *sp = createSPHParticle(pos, velocity);
    sp->isObstacle = false;

    sp->gridID = grid.insertPoint(sp->position);
    std::pair<int,SPHParticle*> pair(sp->gridID, sp);
    particlesByGridID.insert(pair);

    fluidParticles.push_back(sp);
    allParticles.push_back(sp);
}

SPHParticle* SPHFluidSimulation::addObstacleParticle(glm::vec3 pos) {
    SPHParticle *sp = createSPHObstacleParticle(pos);
    sp->isObstacle = true;

    sp->gridID = grid.insertPoint(sp->position);
    std::pair<int,SPHParticle*> pair(sp->gridID, sp);
    particlesByGridID.insert(pair);

    obstacleParticles.push_back(sp);
    allParticles.push_back(sp);

    return sp;
}

inline double SPHFluidSimulation::evaluateSpeedOfSound(SPHParticle *sp) {
    double sqr = ratioOfSpecificHeats*(sp->pressure)/sp->density;
    if (sqr < 0) {
        sqr = -sqr;
    }
    return sqrt(sqr);
}

inline double SPHFluidSimulation::evaluateSpeedOfSoundSquared(SPHParticle *sp) {
    if (sp->density < 0.00001) {
        return 0.0;
    }
    return ratioOfSpecificHeats*(sp->pressure)/sp->density;
}


void SPHFluidSimulation::updateGrid() {
    SPHParticle *sp;
    for (uint i=0; i<allParticles.size(); i++) {
        sp = allParticles[i];
        grid.movePoint(sp->gridID, sp->position);
    }

    grid.update();
}

double SPHFluidSimulation::calculateTimeStep() {
    double maxvsq = 0.0;         // max velocity squared
    double maxcsq = 0.0;         // max speed of sound squared
    double maxasq = 0.0;         // max accelleration squared
    SPHParticle *sp;
    for (uint i=0; i<fluidParticles.size(); i++) {
        sp = fluidParticles[i];
        double vsq = glm::dot(sp->velocity, sp->velocity);
        double asq = glm::dot(sp->acceleration, sp->acceleration);
        double csq = evaluateSpeedOfSoundSquared(sp);
        if (vsq > maxvsq) { maxvsq = vsq; }
        if (csq > maxcsq) { maxcsq = csq; }
        if (asq > maxasq) { maxasq = asq; }
    }

    double maxv = sqrt(maxvsq);
    double maxc = sqrt(maxcsq);
    double maxa = sqrt(maxasq);

    double vStep = courantSafetyFactor*h / fmax(1.0, maxv);
    double cStep = courantSafetyFactor*h / maxc;
    double aStep = sqrt(h/maxa);
    double tempMin = fmin(vStep, cStep);

    //qDebug() << maxv << maxa << maxc;

    return fmax(minTimeStep, fmin(tempMin, aStep));
}

void SPHFluidSimulation::updateNearestNeighbours() {
    SPHParticle *sp;
    for (uint i=0; i<fluidParticles.size(); i++) {
        sp = fluidParticles[i];
        sp->neighbours.clear();
        std::vector<int> refs = grid.getIDsInRadiusOfPoint(sp->gridID, h);
        for (uint j=0; j<refs.size(); j++) {
            sp->neighbours.push_back(particlesByGridID[refs[j]]);
        }
    }
}


void SPHFluidSimulation::enforceFluidParticlePositionBounds(SPHParticle *p) {
    if (p->position.x < xmin) {
        p->position = glm::vec3(xmin, p->position.y, p->position.z);
        p->velocity = glm::vec3(-p->velocity.x, p->velocity.y, p->velocity.z);
    } else if (p->position.x > xmax) {
        p->position = glm::vec3(xmax, p->position.y, p->position.z);
        p->velocity = glm::vec3(-p->velocity.x, p->velocity.y, p->velocity.z);
    }

    if (p->position.y < ymin) {
        p->position = glm::vec3(p->position.x, ymin, p->position.z);
        p->velocity = glm::vec3(p->velocity.x, -p->velocity.y, p->velocity.z);
    } else if (p->position.y > ymax) {
        p->position = glm::vec3(p->position.x, ymax, p->position.z);
        p->velocity = glm::vec3(p->velocity.x, -p->velocity.y, p->velocity.z);
    }

    if (p->position.z < zmin) {
        p->position = glm::vec3(p->position.x, p->position.y, zmin);
        p->velocity = glm::vec3(p->velocity.x, p->velocity.y, -p->velocity.z);
    } else if (p->position.z > zmax) {
        p->position = glm::vec3(p->position.x, p->position.y, zmax);
        p->velocity = glm::vec3(p->velocity.x, p->velocity.y, -p->velocity.z);
    }
}

void SPHFluidSimulation::updateFluidPositionAndDensity(double dt) {
    SPHParticle *p;
    for (uint i=0; i<fluidParticles.size(); i++) {
        p = fluidParticles[i];

        // calculate velocity at half timestep interval for leapfrog integration
        if (p->isHalfTimeStepVelocityInitialized) {
            p->velocityAtHalfTimeStep += (float)dt * p->acceleration;
        } else {
            p->velocityAtHalfTimeStep = p->velocity + (float)(0.5*dt)*p->acceleration;
            p->isHalfTimeStepVelocityInitialized = true;
        }

        // new position calculated with half time step velocity plus xsph variant
        glm::vec3 velocity;
        p->position += (float)dt * velocity;

        // update sph velocity by advancing half time step velocty by 1/2 interval
        p->velocity = p->velocityAtHalfTimeStep + (float)(0.5*dt) * p->acceleration;
        if (glm::length(p->velocity) > maximumVelocity) {
            glm::vec3 unit = glm::normalize(p->velocity);
            p->velocity = (float)maximumVelocity*unit;
        }

        enforceFluidParticlePositionBounds(p);

    }

}

void SPHFluidSimulation::updateBoundaryForces() {
    SPHParticle *sp;
    double r = boundaryForceRadius;
    double minf = minBoundaryForce;
    double maxf = maxBoundaryForce;

    for (uint i=0; i<fluidParticles.size(); i++) {
        sp = fluidParticles[i];
        glm::vec3 p = sp->position;

        if (p.x < xmin + r) {
            double dist = fmax(0.0, p.x - xmin);
            double force = utils::lerp(maxf, minf, dist/r);
            sp->acceleration += glm::vec3(force/sp->mass, 0.0, 0.0);
        }

        if (p.x > xmax - r) {
            double dist = fmax(0.0, xmax - p.x);
            double force = utils::lerp(maxf, minf, dist/r);
            sp->acceleration += glm::vec3(-force/sp->mass, 0.0, 0.0);
        }

        if (p.y < ymin + r) {
            double dist = fmax(0.0, p.y - ymin);
            double force = utils::lerp(maxf, minf, dist/r);
            sp->acceleration += glm::vec3(0.0, force/sp->mass, 0.0);
        }

        if (p.y > ymax - r) {
            double dist = fmax(0.0, ymax - p.y);
            double force = utils::lerp(maxf, minf, dist/r);
            sp->acceleration += glm::vec3(0.0, -force/sp->mass, 0.0);
        }

        if (p.z < zmin + r) {
            double dist = fmax(0.0, p.z - zmin);
            double force = utils::lerp(maxf, minf, dist/r);
            sp->acceleration += glm::vec3(0.0, 0.0, force/sp->mass);
        }

        if (p.z > zmax - r) {
            double dist = fmax(0.0, zmax - p.z);
            double force = utils::lerp(maxf, minf, dist/r);
            sp->acceleration += glm::vec3(0.0, 0.0, -force/sp->mass);
        }
    }
}

void SPHFluidSimulation::update(float dt) {
    int numSteps = 0;

    StopWatch t1 = StopWatch();
    StopWatch t2 = StopWatch();

    t1.start();
    double timeLeft = dt;
    while (timeLeft > 0.0) {

        t2.start();
        updateGrid();
        updateNearestNeighbours();
        t2.stop();

        updateBoundaryForces();



        // calculate next time step
        double timeStep = calculateTimeStep();
        timeLeft -= timeStep;
        if (timeLeft < 0.0) {
            timeStep = timeStep + timeLeft;
            timeLeft = 0.0;
        }

        numSteps += 1;

        updateFluidPositionAndDensity((double)timeStep);
    }
    //qDebug() << numSteps;
    t1.stop();

    //qDebug() << "update:" << t1.getTime() << "neighbours:" << t2.getTime() <<
    //            "pct:" << (t2.getTime()/t1.getTime())*100.0;

}

void SPHFluidSimulation::draw() {
    //grid.draw();

    glColor3f(0.0, 0.3, 1.0);
    glPointSize(5.0);
    glBegin(GL_POINTS);
    for (uint i=0; i<fluidParticles.size(); i++) {
        glm::vec3 p = fluidParticles[i]->position;
        glVertex3f(p.x, p.y, p.z);
    }

    //glColor4f(0.0, 0.0, 0.0, 0.1);
    for (uint i=0; i<obstacleParticles.size(); i++) {
        //glm::vec3 p = obstacleParticles[i]->position;
        //glVertex3f(p.x, p.y, p.z);
    }

    glEnd();

    double w = xmax - xmin;
    double h = ymax - ymin;
    double d = zmax - zmin;
    glColor3f(0.0, 0.0, 0.0);
    utils::drawWireframeCube(glm::vec3(w/2, h/2, d/2), w, h, d);
}










