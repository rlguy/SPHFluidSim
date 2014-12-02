/* CSC 473 Summer 2014
 * Assignment #1
 *
 * Ryan Guy
 * V00484803
 *
 */

#ifndef CAMERA3D_H
#define CAMERA3D_H
#define GLM_FORCE_RADIANS
#include "glm/glm.hpp"
#include "glm/gtx/rotate_vector.hpp"
#include <GL/glu.h>

class camera3d
{
public:
    camera3d();
    camera3d(glm::vec3 pos, glm::vec3 dir,
             float scrWidth, float scrHeight,
             float fovy, float nearPlane, float farPlane);

    glm::mat4 getViewMatrix();
    glm::mat4 getProjectionMatrix();
    glm::vec3 worldToScreenCoordinates(glm::vec3 p);
    glm::vec3 getPosition();
    glm::vec3 castRayFromScreen(double mx, double my);

    glm::vec3 position;
    glm::vec3 direction;
    glm::vec3 up;          // assume <0,1,0> on initialization
    float screenWidth;
    float screenHeight;

    float getFieldOfView();
    float getAspectRatio();
    void set();
    void unset();
    void resize(float width, float height);
    void moveForward(float val);    // move along direction vector
    void moveBackward(float val);
    void moveRight(float val);      // move side to side
    void moveLeft(float val);
    void moveUp(float val);         // move along up vector
    void moveDown(float val);
    void rotateRight(float rad);    // rotaties about up vector
    void rotateLeft(float rad);
    void rotateUp(float rad);       // rotate about right/left vector
    void rotateDown(float rad);
    void rollRight(float rad);      // roll about direction vector
    void rollLeft(float rad);
    void initializeOrientation();
    void setRotation(glm::mat4 rotMatrix);
    void setPosition(glm::vec3 pos);

private:
    float fov;
    float nearDist;
    float farDist;
};

#endif // CAMERA3D_H
