/* CSC 473 Summer 2014
 * Assignment #1
 *
 * Ryan Guy
 * V00484803
 *
 */

#include "camera3d.h"
#include <QDebug>

camera3d::camera3d() {
}

// Camera is initialized with a position and direction
camera3d::camera3d(glm::vec3 pos, glm::vec3 dir,
                   float scrWidth, float scrHeight,
                   float fovy, float nearPlane, float farPlane)
{
    position = pos;
    direction = dir;
    screenWidth = scrWidth; screenHeight = scrHeight;
    fov = fovy;
    nearDist = nearPlane;
    farDist = farPlane;

    // find up direction
    glm::vec3 right = glm::cross(direction, glm::vec3(0.0, 1.0, 0.0));
    up = glm::normalize(glm::cross(right, direction));
}

// returns fov in radians
// fov is in y direction
float camera3d::getFieldOfView() {
    return fov * (3.141592653/180.0);
}

float camera3d::getAspectRatio() {
    return screenWidth/screenHeight;
}

// set up the projection matrix
void camera3d::set() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glViewport(0, 0, screenWidth, screenHeight);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(fov, screenWidth/screenHeight, nearDist, farDist);

    float px = position.x, py = position.y, pz = position.z;
    gluLookAt(px, py, pz,
              px + direction.x, py + direction.y, pz + direction.z,
              up.x, up.y, up.z);
}

// restore matrix to identity
void camera3d::unset() {
    glPopMatrix();
    glLoadIdentity();
}

void camera3d::resize(float width, float height) {
    screenWidth = width;
    screenHeight = height;
}

glm::vec3 camera3d::getPosition() {
    return position;
}

// convert screen mouse coordinates to a direction in world coordinates
glm::vec3 camera3d::castRayFromScreen(double mx, double my) {
    // Does not work correctly when camera is rolled
    /*
    // normalize mouse coordinates to [-1,1] (clip space)
    float x = (2.0*mx / screenWidth) - 1.0;
    float y = 1.0f - (2.0*my) / screenHeight;
    float z = -1.0f;
    glm::vec4 rayClip = glm::vec4(x, y, z, 1.0);

    // to eye coordinates
    glm::vec4 rayEye = rayClip * glm::inverse(getProjectionMatrix());
    rayEye = glm::vec4(rayEye.x, rayEye.y, 1.0, 0.0);

    // to world coordinates
    glm::vec4 rayWorld = rayEye * glm::inverse(getViewMatrix());
    glm::vec3 rayDir = glm::normalize(glm::vec3(rayWorld));
    */

    // Following code section is from:
    // http://www.bfilipek.com/2012/06/select-mouse-opengl.html
    double matModelView[16], matProjection[16];
    int viewport[4];

    // get matrix and viewport:
    glGetDoublev(GL_MODELVIEW_MATRIX, matModelView);
    glGetDoublev(GL_PROJECTION_MATRIX, matProjection);
    glGetIntegerv(GL_VIEWPORT, viewport);

    mx = mx;
    my = viewport[3] - my;

    double x1, y1, z1 = 0;
    double x2, y2, z2 = 0;

    // get point on the 'near' plane (third param is set to 0.0)
    gluUnProject(mx, my, 0.0, matModelView, matProjection,
             viewport, &x1, &y1, &z1);

    // get point on the 'far' plane (third param is set to 1.0)
    gluUnProject(mx, my, 1.0, matModelView, matProjection,
             viewport, &x2, &y2, &z2);

    glm::vec3 start = glm::vec3(x1, y1, z1);
    glm::vec3 end = glm::vec3(x2, y2, z2);
    glm::vec3 dir = glm::normalize(end - start);

    return dir;
}

// projection matrix in row order
glm::mat4 camera3d::getProjectionMatrix() {
    float aspect = screenWidth / screenHeight;
    float n = nearDist;
    float f = farDist;
    float fv = 2 * fov * (3.141592653/180.0);
    return glm::transpose(glm::mat4(
                            fv/aspect, 0 ,0 ,0,
                            0, fv, 0, 0,
                            0, 0, (f+n)/(n-f), (2*f*n)/(n-f),
                            0, 0, -1, 0));
}

// view matrix in row order
glm::mat4 camera3d::getViewMatrix() {
    glm::vec3 right = glm::cross(direction, glm::vec3(0.0, 1.0, 0.0));
    return glm::transpose(glm::mat4(
                              right.x, up.x, direction.x, position.x,
                              right.y, up.y, direction.y, position.y,
                              right.z, up.z, direction.z, position.z,
                              0.0, 0.0, 0.0, 1.0
                         ));
}

// convert world position p to screen coordinates
glm::vec3 camera3d::worldToScreenCoordinates(glm::vec3 p) {
    GLdouble scrX, scrY, scrZ;
    int viewport[4];
    double mvmatrix[16],projmatrix[16];

    // identity
    for (int i=0 ;i<16; i++) {
        mvmatrix[i] = 0.0;
    }
    mvmatrix[0] = 1.0;
    mvmatrix[5] = 1.0;
    mvmatrix[10] = 1.0;
    mvmatrix[15] = 1.0;

    glGetDoublev(GL_PROJECTION_MATRIX, projmatrix);
    glGetIntegerv(GL_VIEWPORT,viewport);

    gluProject((GLdouble)p.x,(GLdouble)p.y,(GLdouble)p.z,
               mvmatrix, projmatrix, viewport,
               &scrX,&scrY,&scrZ);

    return glm::vec3((float)scrX, (float)viewport[3] - (float)scrY - 1.0, 0.0);
}

void camera3d::moveForward(float val) {
    position = position + val * direction;
}

void camera3d::moveBackward(float val) {
    position = position - val * direction;
}

void camera3d::moveRight(float val) {
    glm::vec3 right = glm::cross(direction, up);
    position = position + val * right;
}

void camera3d::moveLeft(float val) {
    glm::vec3 left = -glm::cross(direction, up);
    position = position + val * left;
}

void camera3d::moveUp(float val) {
    position = position + val * up;
}

void camera3d::moveDown(float val) {
    position = position - val * up;
}

void camera3d::rotateRight(float rad) {
    direction = glm::rotate(direction, rad, -up);
}

void camera3d::rotateLeft(float rad) {
    direction = glm::rotate(direction, -rad, -up);
}

void camera3d::rotateUp(float rad) {
    glm::vec3 right = glm::cross(direction, up);
    up = glm::rotate(up, rad, right);
    direction = glm::rotate(direction, rad, right);
}

void camera3d::rotateDown(float rad) {
    glm::vec3 right = glm::cross(direction, up);
    up = glm::rotate(up, -rad, right);
    direction = glm::rotate(direction, -rad, right);
}

void camera3d::rollRight(float rad) {
    up = glm::rotate(up, rad, -direction);
}

void camera3d::rollLeft(float rad) {
    up = glm::rotate(up, -rad, -direction);
}

void camera3d::initializeOrientation() {
    position = glm::vec3(0.0, 0.0, 0.0);
    direction = glm::vec3(1.0, 0.0, 0.0);
    up = glm::vec3(0.0, 1.0, 0.0);
}

void camera3d::setRotation(glm::mat4 rotMatrix) {
    initializeOrientation();
    glm::vec4 dir4 = rotMatrix * glm::vec4(direction.x,
                                           direction.y,
                                           direction.z, 1.0);
    glm::vec4 up4 = rotMatrix * glm::vec4(up.x, up.y, up.z, 1.0);

    direction = glm::vec3(dir4.x, dir4.y, dir4.z);
    up = glm::vec3(up4.x, up4.y, up4.z);
}

void camera3d::setPosition(glm::vec3 pos) {
    position = pos;
}







