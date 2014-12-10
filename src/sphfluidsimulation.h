#ifndef SPHFLUIDSIMULATION_H
#define SPHFLUIDSIMULATION_H

#include <QDebug>
#include <QString>
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
#include "gradients.h"
#include "camera3d.h"

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
    void setTexture(GLuint *tex);
    void setCamera(camera3d *cam);
    QString getTimingData();
    std::vector<SPHParticle*> getFluidParticles();
    std::vector<SPHParticle*> getObstacleParticles();
    std::vector<SPHParticle*> getAllParticles();
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
    void initializeBoundaryParticles();
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

    // graphics
    void updateGraphics(double dt);
    void updateZSortingDistance();
    void updateFluidColor(double dt);
    void updateFluidParticleColorDensity(double dt, SPHParticle *sp);
    void updateFluidParticleAlpha(double dt, SPHParticle *sp);
    glm::vec3 calculateFluidParticleColor(SPHParticle *sp);
    bool isFluidParticleStuckToBoundary(SPHParticle *sp);
    std::vector<std::array<double, 3>> fluidGradient;
    double maxColorVelocity = 1.0;
    double maxColorAcceleration = 1.0;
    double minColorDensity = 0.0;
    double maxColorDensity = 100.0;
    double colorArrivalRadius = 0.5;
    double stuckToBoundaryRadius = 0.01;
    double stuckToBoundaryAlphaVelocity = 1.0;
    bool isTextureInitialized = false;
    GLuint *texture;
    bool isCameraInitialized = false;
    camera3d *camera;

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
    int boundaryObstacleID;
    bool isBoundaryObstacleInitialized = false;
    bool isHiddenBoundaryParticlesEnabled = true;
    bool isBoundaryParticlesEnabled = false;

    SpatialGrid grid;
    std::vector<SPHParticle*> fluidParticles;
    std::vector<SPHParticle*> obstacleParticles;
    std::vector<SPHParticle*> allParticles;
    std::vector<SPHObstacle*> obstacles;
    std::unordered_map<int,SPHParticle*> particlesByGridID;
    std::unordered_map<int,SPHObstacle*> obstaclesByID;
    bool isSPHParticleRemoved = false;
    glm::vec3 cameraPosition;

    // timing metrics
    double neighbourSearchTime = 0.0;
    double simulationTime = 0.0;
    double graphicsUpdateTime = 0.0;
    double graphicsDrawTime = 0.0;

};

#endif // SPHFLUIDSIMULATION_H
