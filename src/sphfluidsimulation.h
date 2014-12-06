#ifndef SPHFLUIDSIMULATION_H
#define SPHFLUIDSIMULATION_H

#include <QDebug>
#include <vector>
#include <time.h>
#include <cmath>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include "glm/glm.hpp"
#include "spatialgrid.h"
#include "sphparticle.h"
#include "stopwatch.h"

class SPHFluidSimulation
{
public:
    SPHFluidSimulation();
    SPHFluidSimulation(double smoothingRadius);
    void addFluidParticles(std::vector<glm::vec3> points);
    void addFluidParticle(glm::vec3 pos, glm::vec3 velocity);
    void update(float dt);
    void draw();

private:
    inline double evaluatePoly6Kernel(double rSquared);
    inline glm::vec3 evaluateSpikeyGradKernel(double r, glm::vec3 pi, glm::vec3 pj);
    inline glm::vec3 evaluateGradKernel(double r, glm::vec3 pi, glm::vec3 pj);
    inline double evaluateKernel(double r);
    inline double evaluatePressureState(SPHParticle *sp);
    inline double evaluateSpeedOfSound(SPHParticle *sp);
    inline double evaluateSpeedOfSoundSquared(SPHParticle *sp);
    double calculateViscosityTerm(SPHParticle *pi, SPHParticle *pj);
    void initSmoothingRadius(double h);
    SPHParticle* createSPHParticle(glm::vec3 pos, glm::vec3 velocity);

    void updatePressure();
    void updateSpeedOfSound();
    void updateGrid();
    void updateNearestNeighbours();
    void updateAccelerationAndDensityRateOfChange();
    void updateXSPHVelocity();
    void updatePositionAndDensity(double dt);
    double calculateTimeStep();

    double h;                              // smoothing radius
    double hsq;                            // radius squared
    double physicalRadiusFactor = 0.4;
    double initialDensity = 1000.0;        // kg/m^3
    double ratioOfSpecificHeats = 1.0;
    double maxDepth = 4.0;
    glm::vec3 gravityForce;
    double poly6Coefficient;
    double spikeyGradCoefficient;
    double pressureStateCoefficient;       // kg/m^2
    double XSPHCoefficient = 0.25;
    double courantSafetyFactor = 1.0;
    double minTimeStep = 1.0/240.0;
    bool isXSPHEnabled = false;
    bool isViscosityEnabled = false;
    double viscosityAlpha = 1;
    double viscosityBeta = 2;

    double maximumAcceleration = 250.0;
    double maximumVelocity = 100.0;

    SpatialGrid grid;
    std::vector<SPHParticle*> fluidParticles;
    std::unordered_map<int,SPHParticle*> particlesByGridID;

};

#endif // SPHFLUIDSIMULATION_H
