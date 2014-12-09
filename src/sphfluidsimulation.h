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
#include "quaternion.h"

extern "C" {
# include "lua/lua.h"
# include "lua/lauxlib.h"
# include "lua/lualib.h"
}
#include "LuaBridge/LuaBridge.h"

class SPHFluidSimulation
{
public:
    SPHFluidSimulation();
    SPHFluidSimulation(double smoothingRadius);
    void update(float dt);
    void draw();

    // user functions
    void addFluidParticles(std::vector<glm::vec3> points);
    void addFluidParticle(glm::vec3 pos, glm::vec3 velocity);
    int addObstacleParticles(std::vector<glm::vec3> points);
    void removeObstacle(int id);

    void setBounds(double xmin, double xmax,
                   double ymin, double ymax,
                   double zmin, double zmax);
    void setDampingConstant(double c);
    std::vector<SPHParticle*> getFluidParticles();
    std::vector<SPHParticle*> getObstacleParticles();
    float getParticleSize();
    float getInitialDensity();

    void setObstaclePosition(int id, glm::vec3 pos);
    void translateObstacle(int id, glm::vec3 trans);
    void rotateObstacle(int id, Quaternion q);

private:

    // init
    void initSimulationConstants();
    void initKernelConstants();

    SPHParticle* createSPHParticle(glm::vec3 pos, glm::vec3 velocity);
    SPHParticle* createSPHObstacleParticle(glm::vec3 pos);
    SPHParticle* addObstacleParticle(glm::vec3 pos);
    int getUniqueObstacleID();
    int currentObstacleID = 0;

    // simulation
    inline double evaluateSpeedOfSound(SPHParticle *sp);
    inline double evaluateSpeedOfSoundSquared(SPHParticle *sp);
    void removeSPHParticlesMarkedForRemoval();
    void updateFluidConstants();
    void updateObstacleVelocity(double dt);
    void updateGrid();
    void updateNearestNeighbours();
    void updateFluidDensityAndPressure();
    glm::vec3 calculateBoundaryAcceleration(SPHParticle *sp);
    void updateFluidAcceleration();
    void updateFluidPosition(double dt);
    void enforceFluidParticlePositionBounds(SPHParticle *p);
    bool isEnforcingFluidParticlePositionBoundsThisTimeStep = false;
    double calculateTimeStep();

    // simulation constants
    double h;                              // smoothing radius
    double hsq;                            // radius squared
    glm::vec3 gravityForce;
    double courantSafetyFactor = 1.0;
    double minTimeStep = 1.0/240.0;
    double gravityMagnitude;
    double initialDensity;
    double pressureCoefficient;
    double particleMass;
    bool isMotionDampingEnabled = false;
    double motionDampingCoefficient;
    double boundaryDampingCoefficient;
    double ratioOfSpecificHeats;
    double viscosityCoefficient;
    double maximumVelocity;
    double maximumAcceleration;
    bool displaySimulationConsoleOutput;

    // kernel constants
    double poly6Coefficient;
    double spikeyGradCoefficient;
    double viscosityLaplacianCoefficient;

    // boundary constraints
    double boundaryForceRadius = 0.1;
    double minBoundaryForce = 0.0;
    double maxBoundaryForce = 0.0;
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
    bool isSPHParticleRemoved = false;
};

#endif // SPHFLUIDSIMULATION_H
