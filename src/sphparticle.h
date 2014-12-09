#ifndef SPHPARTICLE_H
#define SPHPARTICLE_H

#endif // SPHPARTICLE_H

#include <vector>
#include "glm/glm.hpp"


struct SPHParticle {
    glm::vec3 position;
    glm::vec3 prevPosition;   // only used for obstacles
    glm::vec3 velocity;
    glm::vec3 velocityAtHalfTimeStep;
    glm::vec3 XSPHVelocity;
    glm::vec3 acceleration;
    double soundSpeed;
    double mass;
    double density;
    double densityVelocity;
    double pressure;
    std::vector<SPHParticle*> neighbours;
    int gridID;  // used for spatial grid lookup
    bool isHalfTimeStepVelocityInitialized;
    bool isObstacle;
};


