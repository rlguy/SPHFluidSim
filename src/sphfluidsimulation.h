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

    inline double evaluateSpeedOfSound(SPHParticle *sp);
    inline double evaluateSpeedOfSoundSquared(SPHParticle *sp);
    SPHParticle* createSPHParticle(glm::vec3 pos, glm::vec3 velocity);
    SPHParticle* createSPHObstacleParticle(glm::vec3 pos);

    SPHParticle* addObstacleParticle(glm::vec3 pos);
    int getUniqueObstacleID();
    int currentObstacleID = 0;

    void updateGrid();
    void updateNearestNeighbours();
    void updateBoundaryForces();
    void updateFluidPositionAndDensity(double dt);
    void enforceFluidParticlePositionBounds(SPHParticle *p);
    double calculateTimeStep();

    double h;                              // smoothing radius
    double hsq;                            // radius squared
    glm::vec3 gravityForce;
    double gravityMagnitude = 2.0;
    double initialDensity = 1000.0;
    double particleMass = 1.0;
    double motionDampingCoefficient = 1.0;
    double ratioOfSpecificHeats = 1.0;
    double courantSafetyFactor = 1.0;
    double minTimeStep = 1.0/240.0;

    // constraints
    double maximumVelocity = 100.0;

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
