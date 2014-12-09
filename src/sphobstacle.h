#ifndef SPHOBSTACLE_H
#define SPHOBSTACLE_H

#endif // SPHOBSTACLE_H

#include <vector>
#include "glm/glm.hpp"

struct SPHObstacle {
    glm::vec3 position;
    std::vector<SPHParticle*> particles;
    bool isVisible = true;
    int id;
};
