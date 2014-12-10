#ifndef UTILS_H
#define UTILS_H

#include <GL/glu.h>
#include <glm/glm.hpp>
#include <vector>
#include <cmath>
#include "camera3d.h"
#include "quaternion.h"

namespace utils {
    void drawBillboard(camera3d *camera, GLuint *tex, glm::vec3 p, float width);
    void drawGrid();
    void drawWireframeCube(glm::vec3 pos, float size);
    void drawWireframeCube(glm::vec3 pos, float width, float height, float depth);
    void drawCircle(glm::vec3 pos, float r, glm::vec3 axis);
    void drawHalfCircle(glm::vec3 pos, float r,
                        glm::vec3 axis, glm::vec3 plane);
    void getHalfCirclePoints(glm::vec3 pos, float r,
                             glm::vec3 axis, glm::vec3 plane,
                             std::vector<glm::vec3> *points);

    float pointToLineDistance(glm::vec3 p, glm::vec3 o, glm::vec3 dir);
    glm::vec3 closestPointOnLineFromPoint(glm::vec3 p, glm::vec3 o, glm::vec3 dir);
    float lerp(float x1, float x2, float t);
    float smoothstep(float t);
    float easeInOut(float t, float t1, float t2);
    glm::vec3 linePlaneIntersection(glm::vec3 lOrigin, glm::vec3 lDir,
                                    glm::vec3 pOrigin, glm::vec3 pNormal);

    std::vector<glm::vec3> createPointPanel(float width, float height,
                                            float spacing, int numLayers,
                                            glm::vec3 w, glm::vec3 h, bool isStaggered);
    std::vector<glm::vec3> translatePoints(std::vector<glm::vec3> points, glm::vec3 trans);
    std::vector<glm::vec3> rotatePoints(std::vector<glm::vec3> points, Quaternion q);
    std::vector<glm::vec3> mergePoints(std::vector<glm::vec3> points1,
                                       std::vector<glm::vec3> points2);
    void drawPoints(std::vector<glm::vec3> points);

}

#endif // UTILS_H
