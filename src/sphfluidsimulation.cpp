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
    initializeBoundaryParticles();

    cameraPosition = glm::vec3(0.0, 0.0, 0.0);
    fluidGradient = Gradients::getSkyblueGradient();
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
    ratioOfSpecificHeats             = t["ratioOfSpecificHeats"].cast<double>();
    pressureCoefficient              = t["pressureCoefficient"].cast<double>();
    initialDensity                   = t["initialDensity"].cast<double>();
    viscosityCoefficient             = t["viscosityCoefficient"].cast<double>();
    particleMass                     = t["particleMass"].cast<double>();
    maximumVelocity                  = t["maximumVelocity"].cast<double>();
    maximumAcceleration              = t["maximumAcceleration"].cast<double>();
    motionDampingCoefficient         = t["motionDampingCoefficient"].cast<double>();
    boundaryDampingCoefficient       = t["boundaryDampingCoefficient"].cast<double>();
    gravityMagnitude                 = t["gravityMagnitude"].cast<double>();
    isMotionDampingEnabled           = t["isMotionDampingEnabled"].cast<bool>();
    isBoundaryParticlesEnabled       = t["isBoundaryParticlesEnabled"].cast<bool>();
    isHiddenBoundaryParticlesEnabled = t["isHiddenBoundaryParticlesEnabled"].cast<bool>();
    displayConsoleOutput             = t["displaySimulationConsoleOutput"].cast<bool>();

    // graphics
    minColorDensity                = t["minColorDensity"].cast<double>();
    maxColorDensity                = t["maxColorDensity"].cast<double>();
    maxColorVelocity               = t["maxColorVelocity"].cast<double>();
    maxColorAcceleration           = t["maxColorAcceleration"].cast<double>();
    colorArrivalRadius             = t["colorArrivalRadius"].cast<double>();
    stuckToBoundaryRadius          = t["stuckToBoundaryRadius"].cast<double>();
    stuckToBoundaryAlphaVelocity   = t["stuckToBoundaryAlphaVelocity"].cast<double>();

    gravityForce = glm::vec3(0.0, -gravityMagnitude, 0.0);
}

void SPHFluidSimulation::updateFluidConstants() {
    initSimulationConstants();
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
    initializeBoundaryParticles();

    isEnforcingFluidParticlePositionBoundsThisTimeStep = false;
}

void SPHFluidSimulation::setDampingConstant(double c) {
    motionDampingCoefficient = c;
}

void SPHFluidSimulation::setTexture(GLuint *tex) {
    texture = tex;
    isTextureInitialized = true;
}

void SPHFluidSimulation::setCamera(camera3d *cam) {
    camera = cam;
    isCameraInitialized = true;
}

QString SPHFluidSimulation::getTimingData() {
    double total = neighbourSearchTime + simulationTime +
                   graphicsUpdateTime + graphicsDrawTime;

    QString str = QString::number(total) + " " +
                  QString::number(neighbourSearchTime) + " " +
                  QString::number(simulationTime) + " " +
                  QString::number(graphicsUpdateTime) + " " +
                  QString::number(graphicsDrawTime);

    return str;
}

std::vector<SPHParticle*> SPHFluidSimulation::getFluidParticles() {
    return fluidParticles;
}

std::vector<SPHParticle*> SPHFluidSimulation::getObstacleParticles() {
    return obstacleParticles;
}

std::vector<SPHParticle*> SPHFluidSimulation::getAllParticles() {
    return allParticles;
}

float SPHFluidSimulation::getParticleSize() {
    return h;
}

float SPHFluidSimulation::getInitialDensity() {
    return initialDensity;
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

void SPHFluidSimulation::removeObstacle(int id) {
    if (obstaclesByID.find(id) == obstaclesByID.end()) {
        return;
    }
    isSPHParticleRemoved = true;

    SPHObstacle *o = obstaclesByID[id];
    obstaclesByID.erase(id);

    SPHParticle *p;
    for (uint i=0; i<o->particles.size(); i++) {
        p = o->particles[i];
        grid.removePoint(p->gridID);
        particlesByGridID.erase(p->gridID);
        p->isMarkedForRemoval = true;
    }

    for (uint i=0; i<obstacles.size(); i++) {
        SPHObstacle *op = obstacles[i];
        if (op->id == o->id) {
            obstacles.erase(obstacles.begin() + i);
            break;
        }
    }

    delete o;
}

void SPHFluidSimulation::setObstaclePosition(int id, glm::vec3 pos) {
    if (obstaclesByID.find(id) == obstaclesByID.end()) {
        return;
    }

    SPHObstacle *o = obstaclesByID[id];
    translateObstacle(id, pos - o->position);
}

void SPHFluidSimulation::translateObstacle(int id, glm::vec3 trans) {
    if (obstaclesByID.find(id) == obstaclesByID.end()) {
        return;
    }
    SPHObstacle *o = obstaclesByID[id];

    SPHParticle *sp;
    for (uint i=0; i<o->particles.size(); i++) {
        sp = o->particles[i];
        sp->position += trans;
    }

    o->position += trans;
}

void SPHFluidSimulation::rotateObstacle(int id, Quaternion q) {
    if (obstaclesByID.find(id) == obstaclesByID.end()) {
        return;
    }

    SPHObstacle *o = obstaclesByID[id];
    glm::mat4 rot = q.getRotationMatrix();

    SPHParticle *sp;
    glm::vec3 op = o->position;
    glm::vec4 p;
    for (uint i=0; i<o->particles.size(); i++) {
        sp = o->particles[i];
        p = glm::vec4(sp->position.x - op.x,
                      sp->position.y - op.y,
                      sp->position.z - op.z, 1.0);
        p = rot*p;

        sp->position = glm::vec3(p.x + op.x, p.y + op.y, p.z + op.z);
    }

}

int SPHFluidSimulation::getUniqueObstacleID() {
    int id = currentObstacleID;
    currentObstacleID++;
    return id;
}

void SPHFluidSimulation::initializeBoundaryParticles() {
    if (!isBoundaryParticlesEnabled) {
        return;
    }

    if (isBoundaryObstacleInitialized) {
        removeObstacle(boundaryObstacleID);
    }

    std::vector<glm::vec3> obsPoints;
    glm::vec3 xdir = glm::vec3(1.0, 0.0, 0.0);
    glm::vec3 ydir = glm::vec3(0.0, 1.0, 0.0);
    glm::vec3 zdir = glm::vec3(0.0, 0.0, 1.0);
    float xwidth = xmax - xmin;
    float ywidth = ymax - ymin;
    float zwidth = zmax - zmin;
    float pad = h;
    int layers = 1;
    bool isStaggered = true;

    // x-z (ymin) plane
    glm::vec3 o = glm::vec3(xmin + 0.5*xwidth, 0.0, zmin + 0.5*zwidth);
    std::vector<glm::vec3> points = utils::createPointPanel(xwidth, zwidth, pad,
                                                            layers, xdir, zdir, isStaggered);
    points = utils::translatePoints(points, o);
    obsPoints = utils::mergePoints(obsPoints, points);

    // x-z (y-max) plane
    o = glm::vec3(xmin + 0.5*xwidth, ymin + ywidth, zmin + 0.5*zwidth);
    points = utils::createPointPanel(xwidth, zwidth, pad, layers, xdir, zdir, isStaggered);
    points = utils::translatePoints(points, o);
    obsPoints = utils::mergePoints(obsPoints, points);

    // x-y (z-min) plane
    o = glm::vec3(xmin + 0.5*xwidth, ymin + 0.5*xwidth, 0.0);
    points = utils::createPointPanel(xwidth, ywidth, pad, layers, xdir, ydir, isStaggered);
    points = utils::translatePoints(points, o);
    obsPoints = utils::mergePoints(obsPoints, points);

    // x-y (z-max) plane
    o = glm::vec3(xmin + 0.5*xwidth, ymin + 0.5*xwidth, zmin + zwidth);
    points = utils::createPointPanel(xwidth, ywidth, pad, layers, xdir, ydir, isStaggered);
    points = utils::translatePoints(points, o);
    obsPoints = utils::mergePoints(obsPoints, points);

    // y-z (x-min) plane
    o = glm::vec3(0.0 ,ymin + 0.5*xwidth, zmin + 0.5*zwidth);
    points = utils::createPointPanel(ywidth, zwidth, pad, layers, ydir, zdir, isStaggered);
    points = utils::translatePoints(points, o);
    obsPoints = utils::mergePoints(obsPoints, points);

    // y-z (x-max) plane
    o = glm::vec3(xmin + xwidth, ymin + 0.5*xwidth, zmin + 0.5*zwidth);
    points = utils::createPointPanel(ywidth, zwidth, pad, layers, ydir, zdir, isStaggered);
    points = utils::translatePoints(points, o);
    obsPoints = utils::mergePoints(obsPoints, points);

    boundaryObstacleID = addObstacleParticles(obsPoints);

    SPHObstacle *obs = obstaclesByID[boundaryObstacleID];
    obs->isVisible = false;

    for (uint i=0; i<obs->particles.size(); i++) {
        obs->particles[i]->isVisible = false;
    }

    isBoundaryObstacleInitialized = true;
}

SPHParticle* SPHFluidSimulation::createSPHParticle(glm::vec3 pos, glm::vec3 velocity) {
    SPHParticle *s = new SPHParticle();

    s->position = pos;
    s->prevPosition = pos;
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

    // graphics
    s->color = glm::vec3(1.0, 1.0, 1.0);
    s->colorDensity = s->density;

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

void SPHFluidSimulation::removeSPHParticlesMarkedForRemoval() {
    if (allParticles.size() == 0 || !isSPHParticleRemoved) {
        return;
    }

    SPHParticle *p;
    for (int i=(int)fluidParticles.size() - 1; i>=0; i--) {
        if (fluidParticles[i]->isMarkedForRemoval) {
            p = fluidParticles[i];
            fluidParticles.erase(fluidParticles.begin() + i);
        }
    }

    for (int i=(int)obstacleParticles.size() - 1; i>=0; i--) {
        if (obstacleParticles[i]->isMarkedForRemoval) {
            p = obstacleParticles[i];
            obstacleParticles.erase(obstacleParticles.begin() + i);
        }
    }

    for (int i=(int)allParticles.size() - 1; i>=0; i--) {
        if (allParticles[i]->isMarkedForRemoval) {
            p = allParticles[i];
            allParticles.erase(allParticles.begin() + i);
            delete p;
        }
    }

    isSPHParticleRemoved = false;
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
    for (uint i=0; i<allParticles.size(); i++) {
        sp = allParticles[i];
        sp->neighbours.clear();
        std::vector<int> refs = grid.getIDsInRadiusOfPoint(sp->gridID, h);
        for (uint j=0; j<refs.size(); j++) {
            sp->neighbours.push_back(particlesByGridID[refs[j]]);
        }
    }
}

void SPHFluidSimulation::updateObstacleVelocity(double dt) {
    SPHObstacle *obs;
    SPHParticle *sp;
    glm::vec3 trans;
    for (uint i=0; i<obstacles.size(); i++) {
        obs = obstacles[i];

        for (uint j=0; j<obs->particles.size(); j++) {
            sp = obs->particles[j];

            if (sp->position == sp->prevPosition) {
                sp->velocity = glm::vec3(0.0, 0.0, 0.0);
            }

            trans = sp->position - sp->prevPosition;
            double dist = glm::length(trans);
            double eps = 0.00000001;
            if (dist > eps) {
                float speed = fmin((dist/dt), maximumVelocity);
                sp->velocity = (trans / (float)dist) * speed;
            } else {
                sp->velocity = glm::vec3(0.0, 0.0, 0.0);
            }

            sp->prevPosition = sp->position;
        }
    }
}

void SPHFluidSimulation::updateFluidDensityAndPressure() {
    // once we find a particle's density, we can find it's pressure
    SPHParticle *pi, *pj;
    glm::vec3 r;
    for (uint i=0; i<allParticles.size(); i++) {
        pi = allParticles[i];
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
            if (!pj->isObstacle) {
                float lap = viscosityLaplacianCoefficient*diff;
                vdiff = pj->velocity - pi->velocity;
                acc += (float)(viscosityCoefficient*massRatio*(1/pj->density)*lap)*vdiff;
            }
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
    if (!isEnforcingFluidParticlePositionBoundsThisTimeStep) {
        isEnforcingFluidParticlePositionBoundsThisTimeStep = true;
        return;
    }

    double eps = 0.001;
    float d = boundaryDampingCoefficient;
    if (p->position.x < xmin) {
        p->position = glm::vec3(xmin + eps, p->position.y, p->position.z);
        p->velocity = glm::vec3(-d*p->velocity.x, p->velocity.y, p->velocity.z);
    } else if (p->position.x > xmax) {
        p->position = glm::vec3(xmax - eps, p->position.y, p->position.z);
        p->velocity = glm::vec3(-d*p->velocity.x, p->velocity.y, p->velocity.z);
    }

    if (p->position.y < ymin) {
        p->position = glm::vec3(p->position.x, ymin + eps, p->position.z);
        p->velocity = glm::vec3(p->velocity.x, -d*p->velocity.y, p->velocity.z);
    } else if (p->position.y > ymax) {
        p->position = glm::vec3(p->position.x, ymax - eps, p->position.z);
        p->velocity = glm::vec3(p->velocity.x, -d*p->velocity.y, p->velocity.z);
    }

    if (p->position.z < zmin) {
        p->position = glm::vec3(p->position.x, p->position.y, zmin + eps);
        p->velocity = glm::vec3(p->velocity.x, p->velocity.y, -d*p->velocity.z);
    } else if (p->position.z > zmax) {
        p->position = glm::vec3(p->position.x, p->position.y, zmax - eps);
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

void SPHFluidSimulation::updateZSortingDistance() {
    if (!isCameraInitialized || !isTextureInitialized) {
        return;
    }

    glm::vec3 r;
    glm::vec3 cpos = camera->position;
    for (uint i=0; i<allParticles.size(); i++) {
        r = cpos - allParticles[i]->position;
        allParticles[i]->zdistance = glm::dot(r, r);
    }
}

bool SPHFluidSimulation::isFluidParticleStuckToBoundary(SPHParticle *sp) {
    double r = stuckToBoundaryRadius;
    bool isStuck = false;
    glm::vec3 p = sp->position;

    if (p.x < xmin + r || p.x > xmax - r ||
        p.y < ymin + r || p.y > ymax - r ||
        p.z < zmin + r || p.z > zmax - r) {
        isStuck = true;
    }

    return isStuck;
}

void SPHFluidSimulation::updateFluidParticleColorDensity(double dt, SPHParticle *sp) {
    // update color density
    // find color density target
    double targetDensity = sp->density;
    double currentDensity = sp->colorDensity;
    double desired = targetDensity - currentDensity;
    if (fabs(desired) < colorArrivalRadius) {
        double r = fabs(desired) / colorArrivalRadius;
        double newMag = utils::lerp(0, maxColorVelocity, r);
        if (desired > 0) {
            desired = newMag;
        } else {
            desired = -newMag;
        }
    } else {
        if (desired > 0) {
            desired = maxColorVelocity;
        } else {
            desired = -maxColorVelocity;
        }
    }

    // color density acceleration
    double acc = desired - sp->colorVelocity;
    if (fabs(acc) > maxColorAcceleration) {
        if (acc > 0) {
            acc = maxColorAcceleration;
        } else {
            acc = -maxColorAcceleration;
        }
    }

    // update color density velocity and position from acceleration
    sp->colorVelocity += acc*dt;
    sp->colorDensity += sp->colorVelocity*dt;
    sp->colorDensity = fmax(0.0, sp->colorDensity);
}

glm::vec3 SPHFluidSimulation::calculateFluidParticleColor(SPHParticle *sp) {
    int idxOffset = 50;
    int maxIdx = fluidGradient.size() - 1;
    int minIdx = fmin(0 + idxOffset, maxIdx);
    double minD = minColorDensity;
    double maxD = maxColorDensity;
    double invDiff = 1/(maxD - minD);
    glm::vec3 targetColor;
    std::array<double, 3> cmin;
    std::array<double, 3> cmax;

    // update color based upon color density
    double d = sp->colorDensity;
    double r = 1 - (d - minD) * invDiff;
    r = fmax(0.0, r);
    r = fmin(1.0, r);

    // find target color
    double lerpVal = utils::lerp(minIdx, maxIdx, r);
    int minCIdx = floor(lerpVal);
    int maxCIdx = ceil(lerpVal);
    if (minCIdx == maxCIdx) {
        cmin = fluidGradient[minCIdx];
        targetColor = glm::vec3(cmin[0], cmin[1], cmin[2]);
    } else {
        cmin = fluidGradient[minCIdx];
        cmax = fluidGradient[maxCIdx];
        double f = lerpVal - minCIdx;

        targetColor = glm::vec3(utils::lerp(cmin[0], cmax[0], f),
                                utils::lerp(cmin[1], cmax[1], f),
                                utils::lerp(cmin[2], cmax[2], f));
    }

    return targetColor;
}

void SPHFluidSimulation::updateFluidParticleAlpha(double dt, SPHParticle *sp) {
    sp->isStuckInBoundary = isFluidParticleStuckToBoundary(sp);

    if (sp->isStuckInBoundary && sp->boundaryAlphaValue > 0.0) {
        sp->boundaryAlphaValue -= stuckToBoundaryAlphaVelocity*dt;
        sp->boundaryAlphaValue = fmax(0.0, sp->boundaryAlphaValue);
        sp->alpha = utils::smoothstep(sp->boundaryAlphaValue);
    } else if (!sp->isStuckInBoundary && sp->boundaryAlphaValue < 1.0) {
        sp->boundaryAlphaValue += stuckToBoundaryAlphaVelocity*dt;
        sp->boundaryAlphaValue = fmin(1.0, sp->boundaryAlphaValue);
        sp->alpha = utils::smoothstep(sp->boundaryAlphaValue);
    }
}

void SPHFluidSimulation::updateFluidColor(double dt) {
    if (dt == 0) { return; }

    SPHParticle *sp;

    // find target color index in gradient. Lerp between color indices.
    // Move color to target based upon smoothed value of particle density.
    // Particle density is smoothed by keeping track of acceleration and velocity
    // of smoothed value changes
    for (uint i=0; i<fluidParticles.size(); i++) {
        sp = fluidParticles[i];

        updateFluidParticleColorDensity(dt, sp);
        updateFluidParticleAlpha(dt, sp);
        sp->color = calculateFluidParticleColor(sp);

        // update fluid alpha.

    }
}

bool compareByZDistance(const SPHParticle *p1, const SPHParticle *p2) {
    return p1->zdistance > p2->zdistance;
}

void SPHFluidSimulation::updateGraphics(double dt) {
    updateZSortingDistance();
    updateFluidColor(dt);

    std::sort(allParticles.begin(), allParticles.end(), compareByZDistance);
}

void SPHFluidSimulation::update(float dt) {
    updateFluidConstants();
    removeSPHParticlesMarkedForRemoval();

    StopWatch neighbourTimer = StopWatch();
    StopWatch simulationTimer = StopWatch();
    StopWatch graphicsTimer = StopWatch();

    int numSteps = 0;
    double timeLeft = dt;
    while (timeLeft > 0.0) {

        neighbourTimer.start();
        updateGrid();
        updateNearestNeighbours();
        neighbourTimer.stop();

        simulationTimer.start();
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
        updateObstacleVelocity((double)timeStep);
        simulationTimer.stop();
    }
    graphicsTimer.start();
    updateGraphics(dt);
    graphicsTimer.stop();

    neighbourSearchTime = neighbourTimer.getTime();
    simulationTime = simulationTimer.getTime();
    graphicsUpdateTime = graphicsTimer.getTime();

    double total = neighbourSearchTime + simulationTime;
    if (displayConsoleOutput) {
        qDebug() << "total: " << QString::number(total)
                    + " pct neighbour: " +
                    QString::number((100*(neighbourSearchTime/total))) +
                    " #particles: " + QString::number(allParticles.size());
    }

}

void SPHFluidSimulation::drawBounds() {
    double w = xmax - xmin;
    double h = ymax - ymin;
    double d = zmax - zmin;
    glColor3f(0.0, 0.0, 0.0);
    utils::drawWireframeCube(glm::vec3(xmin + w/2, ymin + h/2, zmin + d/2),
                             w, h, d);
}

void SPHFluidSimulation::draw() {
    StopWatch drawTimer = StopWatch();
    drawTimer.start();

    //grid.draw();

    SPHParticle *sp;
    float size = 0.5*h;
    for (uint i=0; i<allParticles.size(); i++) {
        sp = allParticles[i];
        glm::vec3 p = sp->position;

        if (sp->isObstacle) {
            glColor3f(0.5, 0.5, 0.5);
        } else {
            if (isHiddenBoundaryParticlesEnabled) {
                glColor4f(sp->color.x, sp->color.y, sp->color.z, sp->alpha);
            } else {
                glColor3f(sp->color.x, sp->color.y, sp->color.z);
            }
        }

        if (sp->isVisible) {
            utils::drawBillboard(camera, texture, p, size);
        }
    }

    drawBounds();

    drawTimer.stop();
    graphicsDrawTime = drawTimer.getTime();

}










