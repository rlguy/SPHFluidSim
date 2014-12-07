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
#include "sphobstacle.h"
#include "stopwatch.h"

class SPHFluidSimulation
{
public:
    SPHFluidSimulation();
    SPHFluidSimulation(double smoothingRadius);
    void update(float dt);
    void draw();

    void addFluidParticles(std::vector<glm::vec3> points);
    void addFluidParticle(glm::vec3 pos, glm::vec3 velocity);
    int addObstacleParticles(std::vector<glm::vec3> points);

    void setBounds(double xmin, double xmax,
                   double ymin, double ymax,
                   double zmin, double zmax);
    void setDampingConstant(double c);

private:
    inline double evaluatePoly6Kernel(double rSquared);
    inline glm::vec3 evaluateSpikeyGradKernel(double r, glm::vec3 pi, glm::vec3 pj);
    inline glm::vec3 evaluateGradKernel(double r, glm::vec3 pi, glm::vec3 pj);
    inline double evaluateKernel(double r);
    inline double evaluatePressureState(SPHParticle *sp);
    inline double evaluateObstaclePressureState(SPHParticle *sp);
    inline double evaluateSpeedOfSound(SPHParticle *sp);
    inline double evaluateSpeedOfSoundSquared(SPHParticle *sp);
    double calculateViscosityTerm(SPHParticle *pi, SPHParticle *pj);
    void initSmoothingRadius(double h);
    SPHParticle* createSPHParticle(glm::vec3 pos, glm::vec3 velocity);
    SPHParticle* createSPHObstacleParticle(glm::vec3 pos);

    SPHParticle* addObstacleParticle(glm::vec3 pos);
    int getUniqueObstacleID();
    int currentObstacleID = 0;

    void updatePressure();
    void updateSpeedOfSound();
    void updateGrid();
    void updateNearestNeighbours();
    void updateFluidAccelerationAndDensityRateOfChange();
    void updateBoundaryForces();
    void updateObstacleDensityRateOfChange();
    void updateFluidDensity();
    void updateXSPHVelocity();
    void updateFluidPositionAndDensity(double dt);
    void enforceFluidParticlePositionBounds(SPHParticle *p);
    void updateObstacleVelocityAndDensity(double dt);
    double calculateTimeStep();

    double h;                              // smoothing radius
    double hsq;                            // radius squared
    double poly6Coefficient;
    double spikeyGradCoefficient;
    double particleMass = 0.3;
    double initialDensity = 1000.0;        // kg/m^3
    double ratioOfSpecificHeats = 1.0;
    glm::vec3 gravityForce;
    double gravityMagnitude = 2.0;
    double pressureStateCoefficient = 100000;       // kg/m^2

    double courantSafetyFactor = 1.0;
    double minTimeStep = 1.0/240.0;

    bool isDensityRateOfChangeEnabled = true;
    bool isDensityInitializationStaggered = true;

    bool isXSPHEnabled = false;
    double XSPHCoefficient = 0.25;

    bool isViscosityEnabled = false;
    double viscosityAlpha = 1;
    double viscosityBeta = 2;

    bool isMotionDampingEnabled = true;
    double motionDampingCoefficient = 8.0;

    double maximumAcceleration = 250.0;
    double maximumVelocity = 100.0;
    double maximumDensityVelocity = 100.0;
    double minimumDensity = 0.1;

    // boundary
    double boundaryForceRadius = 0.1;
    double minBoundaryForce = 0.0;
    double maxBoundaryForce = 0.8;
    double xmin = 0.0;
    double xmax = 1.0;
    double ymin = 0.0;
    double ymax = 1.0;
    double zmin = 0.0;
    double zmax = 1.0;

    SpatialGrid grid;
    std::vector<SPHParticle*> fluidParticles;
    std::vector<SPHParticle*> obstacleParticles;
    std::vector<SPHParticle*> allParticles;
    std::vector<SPHObstacle*> obstacles;
    std::unordered_map<int,SPHParticle*> particlesByGridID;
    std::unordered_map<int,SPHObstacle*> obstaclesByID;

};

#endif // SPHFLUIDSIMULATION_H
