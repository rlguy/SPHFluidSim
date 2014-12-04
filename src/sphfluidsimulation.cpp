#include "sphfluidsimulation.h"

SPHFluidSimulation::SPHFluidSimulation()
{
    h = 1.0;
}

SPHFluidSimulation::SPHFluidSimulation(double smoothingRadius)
{
    initSmoothingRadius(smoothingRadius);
    grid = SpatialGrid(2*h);
    gravityForce = glm::vec3(0.0, -9.8, 0.0);

    double B = 200*initialDensity*glm::length(gravityForce)*maxDepth / ratioOfSpecificHeats;
    pressureStateCoefficient = B;
}

void SPHFluidSimulation::initSmoothingRadius(double hh) {
    h = hh;
    hsq = h*h;
    double hpow6 = h*h*h*h*h*h;
    double hpow9 = h*h*h*h*h*h*h*h*h;
    double pi = 3.1415926535897;

    poly6Coefficient = 315.0/(64.0*pi*hpow9);
    spikeyGradCoefficient = -45.0/(pi*hpow6);
}

void SPHFluidSimulation::addFluidParticles(std::vector<glm::vec3> points) {
    for (uint i=0; i<points.size(); i++) {
        addFluidParticle(points[i], glm::vec3(0.0, 0.0, 0.0));
    }
}

SPHParticle* SPHFluidSimulation::createSPHParticle(glm::vec3 pos, glm::vec3 velocity) {
    SPHParticle *s = new SPHParticle();

    s->position = pos;
    s->velocity = velocity;
    s->acceleration = glm::vec3(0.0, 0.0, 0.0);

    // Create pressure offset from initial pressure
    // uniform densities will cause uniform pressure of 0, meaning no acceleration
    // of system
    double pressureOffset = 2*((double)rand() / (double)RAND_MAX) - 1;
    s->density = initialDensity + pressureOffset;

    // mass of sphere
    double radius = h * physicalRadiusFactor;
    s->mass = (4.0/3.0)*3.141592653*radius*radius*radius*s->density;

    // initial pressure will be calculated once all particles are in place
    s->pressure = 0.0;

    s->gridID = grid.insertPoint(s->position);

    return s;
}

void SPHFluidSimulation::addFluidParticle(glm::vec3 pos, glm::vec3 velocity) {
    SPHParticle *sp = createSPHParticle(pos, velocity);
    fluidParticles.push_back(sp);

    std::pair<int,SPHParticle*> pair(sp->gridID, sp);
    particlesByGridID.insert(pair);
}


inline double SPHFluidSimulation::evaluatePoly6Kernel(double rsq) {
    double diff = (hsq - rsq);
    return poly6Coefficient*diff*diff*diff;
}

inline glm::vec3 SPHFluidSimulation::evaluateSpikeyGradKernel(
                                 double r, glm::vec3 pi, glm::vec3 pj) {
    return (float)(spikeyGradCoefficient*(h-r))*(pi-pj);
}

inline double SPHFluidSimulation::evaluatePressureState(SPHParticle *sp) {
    return pressureStateCoefficient*(
                pow(sp->density/initialDensity, ratioOfSpecificHeats)-1.0);
}

inline double SPHFluidSimulation::evaluateSpeedOfSoundSquared(SPHParticle *sp) {
    return ratioOfSpecificHeats*sp->pressure/sp->density;
}

void SPHFluidSimulation::updatePressure() {
    SPHParticle *sp;
    for (uint i=0; i<fluidParticles.size(); i++) {
        sp = fluidParticles[i];
        sp->pressure = evaluatePressureState(sp);
    }
}

void SPHFluidSimulation::updateGrid() {
    SPHParticle *sp;
    for (uint i=0; i<fluidParticles.size(); i++) {
        sp = fluidParticles[i];
        grid.movePoint(sp->gridID, sp->position);
    }
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

    return fmin(tempMin, aStep);
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

void SPHFluidSimulation::update(float dt) {
    int numSteps = 0;

    double timeLeft = dt;
    while (timeLeft > 0.0) {
        updatePressure();
        updateGrid();
        updateNearestNeighbours();

        // calculate next time step
        double timeStep = calculateTimeStep();
        timeLeft -= timeStep;
        if (timeLeft < 0.0) {
            timeStep = timeStep + timeLeft;
            timeLeft = 0.0;
        }
        numSteps += 1;
    }

}

void SPHFluidSimulation::draw() {
    grid.draw();
}










