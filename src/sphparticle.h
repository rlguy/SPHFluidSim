#ifndef SPHPARTICLE_H
#define SPHPARTICLE_H

#endif // SPHPARTICLE_H

#include <vector>
#include "glm/glm.hpp"


struct SPHParticle {
    glm::vec3 position;
    glm::vec3 velocity;
    glm::vec3 acceleration;
    double mass;
    double density;
    double pressure;
    std::vector<SPHParticle*> neighbours;
    int gridID;  // used for spatial grid lookup
};


