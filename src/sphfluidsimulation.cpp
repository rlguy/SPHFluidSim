#include "sphfluidsimulation.h"

SPHFluidSimulation::SPHFluidSimulation()
{
    h = 1.0;
}

SPHFluidSimulation::SPHFluidSimulation(double smoothingRadius)
{
    initSmoothingRadius(smoothingRadius);
    grid = SpatialGrid(2*h);
    gravityForce = glm::vec3(0.0, -3.0, 0.0);

    double B = 200*initialDensity*glm::length(gravityForce)*maxDepth/ratioOfSpecificHeats;
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
    s->velocityAtHalfTimeStep = glm::vec3(0.0, 0.0, 0.0);
    s->isHalfTimeStepVelocityInitialized = false;
    s->XSPHVelocity = glm::vec3(0.0, 0.0, 0.0);
    s->soundSpeed = 0;
    s->acceleration = glm::vec3(0.0, 0.0, 0.0);

    // Create pressure offset from initial pressure
    // uniform densities will cause uniform pressure of 0, meaning no acceleration
    // of system
    double pressureOffset = 2*((double)rand() / (double)RAND_MAX) - 1;
    s->density = initialDensity + pressureOffset;
    s->densityVelocity = 0.0;

    // mass of sphere
    double radius = h * physicalRadiusFactor;
    s->mass = (4.0/3.0)*3.141592653*radius*radius*radius*s->density;

    // initial pressure will be calculated once all particles are in place
    s->pressure = 0.0;

    return s;
}

void SPHFluidSimulation::addFluidParticle(glm::vec3 pos, glm::vec3 velocity) {
    SPHParticle *sp = createSPHParticle(pos, velocity);

    sp->gridID = grid.insertPoint(sp->position);
    std::pair<int,SPHParticle*> pair(sp->gridID, sp);
    particlesByGridID.insert(pair);

    fluidParticles.push_back(sp);
}


inline double SPHFluidSimulation::evaluatePoly6Kernel(double rsq) {
    double diff = (hsq - rsq);
    return poly6Coefficient*diff*diff*diff;
}

inline glm::vec3 SPHFluidSimulation::evaluateSpikeyGradKernel(
                                 double r, glm::vec3 pi, glm::vec3 pj) {
    return (float)(spikeyGradCoefficient*(h-r)*(h-r))*(pi-pj);
}

inline glm::vec3 SPHFluidSimulation::evaluateGradKernel(
                                 double r, glm::vec3 pi, glm::vec3 pj) {
    if (r <= h) {
        return (float)(3*r*(3*r-4*h)/(4*3.141592653*h*h*h*h*h*h))*(pi-pj);
    } else if (r <= 2*h) {
        return (float)(-3*(r-2*h)*(r-2*h)/(4*3.141592653*h*h*h*h*h*h))*(pi-pj);
    }

    return glm::vec3(0.0, 0.0, 0.0);
}

inline double SPHFluidSimulation::evaluateKernel(double r) {
    double q = r/h;
    if (r <= h) {
        return (1/(3.141592653*h*h*h))*(1 - (3/2)*q*q + (3/4)*q*q*q);
    } else if (r <= 2*h) {
        return (1/(3.141592653*h*h*h))*(0.25*(2-q)*(2-q)*(2-q));
    }

    return 0.0;
}

inline double SPHFluidSimulation::evaluatePressureState(SPHParticle *sp) {
    return pressureStateCoefficient*(
                pow(sp->density/initialDensity, ratioOfSpecificHeats)-1.0);
}

inline double SPHFluidSimulation::evaluateSpeedOfSound(SPHParticle *sp) {
    double sqr = ratioOfSpecificHeats*(sp->pressure)/sp->density;
    if (sqr < 0) {
        sqr = -sqr;
    }
    return sqrt(sqr);
}

double SPHFluidSimulation::calculateViscosityTerm(SPHParticle *pi, SPHParticle *pj) {
    glm::vec3 pdiff = pi->position - pj->position;
    glm::vec3 vdiff = pi->velocity - pj->velocity;
    double dot = glm::dot(pdiff, vdiff);

    if (dot > 0.0) {
        return 0.0;
    }

    double u = h*dot / (glm::dot(pdiff, pdiff) + 0.01*h*h);
    double meanc = 0.5 * (pi->soundSpeed + pj->soundSpeed);
    double meanp = 0.5 * (pi->density + pj->density);

    return -(viscosityAlpha*meanc*u + viscosityBeta*u*u) /  meanp;
}

void SPHFluidSimulation::updatePressure() {
    SPHParticle *sp;
    for (uint i=0; i<fluidParticles.size(); i++) {
        sp = fluidParticles[i];
        sp->pressure = evaluatePressureState(sp);
    }
}

void SPHFluidSimulation::updateSpeedOfSound() {
    SPHParticle *sp;
    for (uint i=0; i<fluidParticles.size(); i++) {
        sp = fluidParticles[i];
        sp->soundSpeed = evaluateSpeedOfSound(sp);
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
        double csq = sp->soundSpeed * sp->soundSpeed;

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

    return fmax(minTimeStep, fmin(tempMin, aStep));
}

void SPHFluidSimulation::updateNearestNeighbours() {
    SPHParticle *sp;
    for (uint i=0; i<fluidParticles.size(); i++) {
        sp = fluidParticles[i];
        sp->neighbours.clear();
        std::vector<int> refs = grid.getIDsInRadiusOfPoint(sp->gridID, 2*h);
        for (uint j=0; j<refs.size(); j++) {
            sp->neighbours.push_back(particlesByGridID[refs[j]]);
        }
    }
}

void SPHFluidSimulation::updateAccelerationAndDensityRateOfChange() {
    SPHParticle *pi;
    SPHParticle *pj;
    glm::vec3 acc;
    double dp;
    glm::vec3 w;
    for (uint i=0; i<fluidParticles.size(); i++) {
        pi = fluidParticles[i];
        double t1 = pi->pressure/(pi->density*pi->density);
        acc = glm::vec3(0.0, 0.0, 0.0);
        dp = 0.0;

        for (uint j=0; j<pi->neighbours.size(); j++) {
            pj = pi->neighbours[j];

            // kernel is the same for acceleration and density rate of change
            double radius = glm::length(pi->position - pj->position);
            w = evaluateGradKernel(radius, pi->position, pj->position);

            // sum acceleration
            double coef = pj->mass*(t1 + pj->pressure/(pj->density*pj->density));
            if (isViscosityEnabled) {
                coef += calculateViscosityTerm(pi, pj);
            }
            acc += (float)coef*w;

            // sum density rate of change
            dp += glm::dot((float)pj->mass*(pi->velocity - pj->velocity), w);
        }
        pi->acceleration = -acc + gravityForce;
        //pi->acceleration = -acc;
        pi->densityVelocity = dp;
    }
}

void SPHFluidSimulation::updateXSPHVelocity() {
    SPHParticle *pi;
    SPHParticle *pj;
    glm::vec3 vsum = glm::vec3(0.0, 0.0, 0.0);
    glm::vec3 pdiff;
    for (uint i=0; i<fluidParticles.size(); i++) {
        pi = fluidParticles[i];

        for (uint j=0; j<pi->neighbours.size(); j++) {
            pj = pi->neighbours[j];

            pdiff = pi->position - pj->position;
            double w = evaluateKernel(glm::length(pdiff));
            double pavg = 0.5 * (pi->density + pj->density);

            vsum += (float)(w * pj->mass / pavg)*(pi->velocity - pj->velocity);
        }

        pi->XSPHVelocity = (float)XSPHCoefficient*vsum;
    }
}

void SPHFluidSimulation::updatePositionAndDensity(double dt) {
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
        if (isXSPHEnabled) {
            velocity = p->velocityAtHalfTimeStep + p->XSPHVelocity;
        } else {
            velocity = p->velocityAtHalfTimeStep;
        }
        p->position += (float)dt * velocity;

        // update sph velocity by advancing half time step velocty by 1/2 interval
        p->velocity = p->velocityAtHalfTimeStep + (float)(0.5*dt) * p->acceleration;

        // update density
        p->density += (float)dt * p->densityVelocity;
        if (p->density < 0.0) {
            p->density = 0.0;
        }


        // debug
        if (p->position.y < 0.1) {
            p->position = glm::vec3(p->position.x, 0.1, p->position.z);
            p->velocity = glm::vec3(p->velocity.x, -p->velocity.y, p->velocity.z);
        }


        /*
        p->velocity += p->acceleration*(float)dt;
        p->position += p->velocity * (float)dt;
        p->density += p->densityVelocity*(float)dt;
        */


        //qDebug() <<

    }
}

void SPHFluidSimulation::update(float dt) {
    int numSteps = 0;

    double timeLeft = dt;
    while (timeLeft > 0.0) {
        updatePressure();
        updateSpeedOfSound();
        updateGrid();
        updateNearestNeighbours();

        updateAccelerationAndDensityRateOfChange();
        updateXSPHVelocity();

        // calculate next time step
        double timeStep = calculateTimeStep();
        timeLeft -= timeStep;
        if (timeLeft < 0.0) {
            timeStep = timeStep + timeLeft;
            timeLeft = 0.0;
        }

        numSteps += 1;

        updatePositionAndDensity((double)timeStep);
    }
    qDebug() << numSteps;

}

void SPHFluidSimulation::draw() {
    grid.draw();

    glColor3f(0.0, 0.3, 1.0);
    glPointSize(3.0);
    glBegin(GL_POINTS);
    for (uint i=0; i<fluidParticles.size(); i++) {
        glm::vec3 p = fluidParticles[i]->position;
        glVertex3f(p.x, p.y, p.z);
    }
    glEnd();
}










