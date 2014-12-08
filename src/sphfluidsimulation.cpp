#include "sphfluidsimulation.h"


SPHFluidSimulation::SPHFluidSimulation()
{
    h = 1.0;
}

SPHFluidSimulation::SPHFluidSimulation(double smoothingRadius)
{
    h = smoothingRadius;
    grid = SpatialGrid(h);

    initSimulationConstants();
    initKernelConstants();
}

void SPHFluidSimulation::initSimulationConstants() {
    using namespace luabridge;

    lua_State* L = luaL_newstate();
    luaL_openlibs(L);

    if (luaL_dofile(L, "scripts/fluid_config.lua") != 0) {
        qDebug() << "Error loading script";
        exit(1);
    }
    LuaRef t = getGlobal(L, "settings");

    hsq = h*h;
    ratioOfSpecificHeats        = t["ratioOfSpecificHeats"].cast<double>();
    pressureCoefficient         = t["pressureCoefficient"].cast<double>();
    initialDensity              = t["initialDensity"].cast<double>();
    viscosityCoefficient        = t["viscosityCoefficient"].cast<double>();
    particleMass                = t["particleMass"].cast<double>();
    maximumVelocity             = t["maximumVelocity"].cast<double>();
    maximumAcceleration         = t["maximumAcceleration"].cast<double>();
    motionDampingCoefficient    = t["motionDampingCoefficient"].cast<double>();
    boundaryDampingCoefficient  = t["boundaryDampingCoefficient"].cast<double>();
    gravityMagnitude            = t["gravityMagnitude"].cast<double>();
    isMotionDampingEnabled      = t["isMotionDampingEnabled"].cast<bool>();

    gravityForce = glm::vec3(0.0, -gravityMagnitude, 0.0);
}

void SPHFluidSimulation::initKernelConstants() {
    double pi = 3.1415926535897;

    poly6Coefficient = 315.0/(64.0*pi*powf(h, 9.0));
    spikeyGradCoefficient = -45.0/(pi*powf(h, 6.0));
    viscosityLaplacianCoefficient = 45.0/(pi*powf(h, 6.0f));
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

std::vector<SPHParticle*> SPHFluidSimulation::getFluidParticles() {
    return fluidParticles;
}

float SPHFluidSimulation::getParticleSize() {
    return h;
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

void SPHFluidSimulation::updateFluidDensityAndPressure() {
    // once we find a particle's density, we can find it's pressure
    SPHParticle *pi, *pj;
    glm::vec3 r;
    for (uint i=0; i<fluidParticles.size(); i++) {
        pi = fluidParticles[i];
        double density = 0.0;

        for (uint j=0; j<pi->neighbours.size(); j++) {
            pj = pi->neighbours[j];
            r = pi->position - pj->position;
            double distsq = glm::dot(r, r);
            double diff = hsq - distsq;
            density += pj->mass*poly6Coefficient*diff*diff*diff;
        }

        pi->density = fmax(density, initialDensity);  // less than initial density
                                                      // produces negative pressures
        pi->pressure = pressureCoefficient*(pi->density - initialDensity);


    }
}

void SPHFluidSimulation::updateFluidAcceleration() {
    // we know particle density and pressure, so acceleration can be found
    SPHParticle *pi, *pj;
    glm::vec3 acc;
    glm::vec3 r;
    glm::vec3 vdiff;
    for (uint i=0; i<fluidParticles.size(); i++) {
        pi = fluidParticles[i];
        acc = glm::vec3(0.0, 0.0, 0.0);

        for (uint j=0; j<pi->neighbours.size(); j++) {
            pj = pi->neighbours[j];
            r = pi->position - pj->position;
            double dist = glm::length(r);

            if (dist == 0.0) { continue; }
            float inv = 1/dist;
            r = inv*r;

            // acceleration due to pressure
            float diff = h - dist;
            float spikey = spikeyGradCoefficient*diff*diff;
            float massRatio = pj->mass/pi->mass;
            float pterm = (pi->pressure + pj->pressure) / (2*pi->density*pj->density);
            acc -= (float)(massRatio*pterm*spikey)*r;

            // acceleration due to viscosity
            float lap = viscosityLaplacianCoefficient*diff;
            vdiff = pj->velocity - pi->velocity;
            acc += (float)(viscosityCoefficient*massRatio*(1/pj->density)*lap)*vdiff;
        }

        // acceleration due to gravity
        acc += gravityForce;

        // acceleration due to simulation bounds
        acc += calculateBoundaryAcceleration(pi);

        // motion damping;
        double mag = glm::length(acc);
        if (isMotionDampingEnabled) {
            glm::vec3 damp = pi->velocity * (float) motionDampingCoefficient;
            if (glm::length(damp) > mag) {
                acc = glm::vec3(0.0, 0.0, 0.0);
            } else {
                acc -= damp;
            }
        }


        if (mag > maximumAcceleration) {
            acc = (acc / (float)mag) * (float)maximumAcceleration;
        }

        pi->acceleration = acc;

    }
}

void SPHFluidSimulation::enforceFluidParticlePositionBounds(SPHParticle *p) {

    float d = boundaryDampingCoefficient;
    if (p->position.x < xmin) {
        p->position = glm::vec3(xmin, p->position.y, p->position.z);
        p->velocity = glm::vec3(-d*p->velocity.x, p->velocity.y, p->velocity.z);
    } else if (p->position.x > xmax) {
        p->position = glm::vec3(xmax, p->position.y, p->position.z);
        p->velocity = glm::vec3(-d*p->velocity.x, p->velocity.y, p->velocity.z);
    }

    if (p->position.y < ymin) {
        p->position = glm::vec3(p->position.x, ymin, p->position.z);
        p->velocity = glm::vec3(p->velocity.x, -d*p->velocity.y, p->velocity.z);
    } else if (p->position.y > ymax) {
        p->position = glm::vec3(p->position.x, ymax, p->position.z);
        p->velocity = glm::vec3(p->velocity.x, -d*p->velocity.y, p->velocity.z);
    }

    if (p->position.z < zmin) {
        p->position = glm::vec3(p->position.x, p->position.y, zmin);
        p->velocity = glm::vec3(p->velocity.x, p->velocity.y, -d*p->velocity.z);
    } else if (p->position.z > zmax) {
        p->position = glm::vec3(p->position.x, p->position.y, zmax);
        p->velocity = glm::vec3(p->velocity.x, p->velocity.y, -d*p->velocity.z);
    }
}

void SPHFluidSimulation::updateFluidPosition(double dt) {
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

        // new position calculated with half time step for leap frog integration
        p->position += (float)dt * p->velocityAtHalfTimeStep;

        // update sph velocity by advancing half time step velocty by 1/2 interval
        p->velocity = p->velocityAtHalfTimeStep + (float)(0.5*dt) * p->acceleration;
        if (glm::length(p->velocity) > maximumVelocity) {
            glm::vec3 unit = glm::normalize(p->velocity);
            p->velocity = (float)maximumVelocity*unit;
        }

        enforceFluidParticlePositionBounds(p);
    }

}

glm::vec3 SPHFluidSimulation::calculateBoundaryAcceleration(SPHParticle *sp) {
    double r = boundaryForceRadius;
    double minf = minBoundaryForce;
    double maxf = maxBoundaryForce;

    glm::vec3 p = sp->position;
    glm::vec3 acceleration = glm::vec3(0.0, 0.0, 0.0);

    if (p.x < xmin + r) {
        double dist = fmax(0.0, p.x - xmin);
        double force = utils::lerp(maxf, minf, dist/r);
        acceleration += glm::vec3(force/sp->mass, 0.0, 0.0);
    } else if (p.x > xmax - r) {
        double dist = fmax(0.0, xmax - p.x);
        double force = utils::lerp(maxf, minf, dist/r);
        acceleration += glm::vec3(-force/sp->mass, 0.0, 0.0);
    }

    if (p.y < ymin + r) {
        double dist = fmax(0.0, p.y - ymin);
        double force = utils::lerp(maxf, minf, dist/r);
        acceleration += glm::vec3(0.0, force/sp->mass, 0.0);
    } else if (p.y > ymax - r) {
        double dist = fmax(0.0, ymax - p.y);
        double force = utils::lerp(maxf, minf, dist/r);
        acceleration += glm::vec3(0.0, -force/sp->mass, 0.0);
    }

    if (p.z < zmin + r) {
        double dist = fmax(0.0, p.z - zmin);
        double force = utils::lerp(maxf, minf, dist/r);
        acceleration += glm::vec3(0.0, 0.0, force/sp->mass);
    } else if (p.z > zmax - r) {
        double dist = fmax(0.0, zmax - p.z);
        double force = utils::lerp(maxf, minf, dist/r);
        acceleration += glm::vec3(0.0, 0.0, -force/sp->mass);
    }

    return acceleration;
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

        updateFluidDensityAndPressure();
        updateFluidAcceleration();

        // calculate next time step
        double timeStep = calculateTimeStep();
        timeLeft -= timeStep;
        if (timeLeft < 0.0) {
            timeStep = timeStep + timeLeft;
            timeLeft = 0.0;
        }
        numSteps += 1;

        updateFluidPosition((double)timeStep);
    }
    //qDebug() << numSteps;
    t1.stop();

    //qDebug() << "update:" << t1.getTime() << "neighbours:" << t2.getTime() <<
    //           "pct:" << (t2.getTime()/t1.getTime())*100.0;

}

void SPHFluidSimulation::draw() {
    //grid.draw();

    /*
    glColor3f(0.0, 0.3, 1.0);
    glPointSize(5.0);
    glBegin(GL_POINTS);
    for (uint i=0; i<fluidParticles.size(); i++) {
        glm::vec3 p = fluidParticles[i]->position;
        glVertex3f(p.x, p.y, p.z);
    }
    */

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










