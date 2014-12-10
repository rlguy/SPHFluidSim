/* CSC 473 Summer 2014
 * Assignment #1
 *
 * Ryan Guy
 * V00484803
 *
 */

/*
    CSC 486 Assignment 1
    Ryan Guy
    V00484803
*/

/****************************************************************************
**
** Copyright (C) 2013 Digia Plc and/or its subsidiary(-ies).
** Contact: http://www.qt-project.org/legal
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** You may use this file under the terms of the BSD license as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of Digia Plc and its Subsidiary(-ies) nor the names
**     of its contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
** $QT_END_LICENSE$
**
****************************************************************************/

#include <QtWidgets>
#include <QtOpenGL>
#include <QCursor>
#include <GL/glu.h>
#include <algorithm>
#include "glwidget.h"
#include "glm/glm.hpp"
#include "quaternion.h"
#include "gradients.h"

#ifndef GL_MULTISAMPLE
#define GL_MULTISAMPLE  0x809D
#endif

//! [0]
GLWidget::GLWidget(QWidget *parent)
    : QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
{

    screenWidth = 1280;
    screenHeight = 720;

    // update timer
    float updatesPerSecond = 60;
    updateTimer = new QTimer(this);
    connect(updateTimer, SIGNAL(timeout()), this, SLOT(updateSimulation()));
    updateTimer->start(1000.0/updatesPerSecond);

    deltaTimer = new QTime();
    deltaTimer->start();

    // Initialize camera
    glm::vec3 pos = glm::vec3(14.0, 15.0, 33.0);
    glm::vec3 dir = glm::normalize(glm::vec3(-pos.x, -pos.y, -pos.z));
    camera = camera3d(pos, dir, screenWidth, screenHeight,
                      60.0, 0.5, 100.0);

    // for updating camera movement
    isMovingForward = false;
    isMovingBackward = false;
    isMovingRight = false;
    isMovingLeft = false;
    isMovingUp = false;
    isMovingDown = false;
    isRotatingRight = false;
    isRotatingLeft = false;
    isRotatingUp = false;
    isRotatingDown = false;
    isRollingRight = false;
    isRollingLeft = false;

    // simulation system
    minDeltaTimeModifier = 0.125;
    maxDeltaTimeModifier = 1.0;
    deltaTimeModifier = maxDeltaTimeModifier;
    runningTime = 0.0;
    currentFrame = 0;

    initializeSimulation();

}

GLWidget::~GLWidget()
{
}

QSize GLWidget::minimumSizeHint() const
{
    return QSize(50, 50);
}

QSize GLWidget::sizeHint() const
{
    return QSize(screenWidth, screenHeight);
}

void GLWidget::initializeGL()
{

    QPixmap piximg = QPixmap("images/blankball.png");
    if (piximg.isNull()) {
        qDebug() << "pixmap is null";
        exit(1);
    }
    texture[0] = bindTexture(piximg, GL_TEXTURE_2D);

    static const GLfloat lightPos[4] = { 20.0f, 20.0f, 20.0f, 1.0f };
    glLightfv(GL_LIGHT0, GL_POSITION, lightPos);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_NORMALIZE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

}


// move and rotate camera
void GLWidget::updateCameraMovement(float dt) {
    float camSpeed = 10.0;
    float camRotSpeed = 2.0;

    if (isMovingForward) { camera.moveForward(camSpeed * dt); }
    if (isMovingBackward) { camera.moveBackward(camSpeed * dt); }
    if (isMovingRight) { camera.moveRight(camSpeed * dt); }
    if (isMovingLeft) { camera.moveLeft(camSpeed * dt); }
    if (isMovingUp) { camera.moveUp(camSpeed * dt); }
    if (isMovingDown) { camera.moveDown(camSpeed * dt); }
    if (isRotatingRight) { camera.rotateRight(camRotSpeed * dt); }
    if (isRotatingLeft) { camera.rotateLeft(camRotSpeed * dt); }
    if (isRotatingUp) { camera.rotateUp(camRotSpeed * dt); }
    if (isRotatingDown) { camera.rotateDown(camRotSpeed * dt); }
    if (isRollingRight) { camera.rollRight(camRotSpeed * dt); }
    if (isRollingLeft) { camera.rollLeft(camRotSpeed * dt); }
}

void GLWidget::initializeSimulation() {
    using namespace luabridge;

    lua_State* L = luaL_newstate();
    luaL_openlibs(L);

    if (luaL_dofile(L, "scripts/fluid_config.lua") != 0) {
        qDebug() << "Error loading script";
        exit(1);
    }
    LuaRef t = getGlobal(L, "settings");

    simulationFPS = t["fps"].cast<double>();
    isSimulationPaused = t["isSimulationPaused"].cast<bool>();


    double radius = t["smoothingRadius"].cast<double>();
    fluidSim = SPHFluidSimulation(radius);

    LuaRef bounds = t["initialBounds"];
    double minx = bounds["minx"].cast<double>();
    double maxx = bounds["maxx"].cast<double>();
    double miny = bounds["miny"].cast<double>();
    double maxy = bounds["maxy"].cast<double>();
    double minz = bounds["minz"].cast<double>();
    double maxz = bounds["maxz"].cast<double>();
    double rx = 1.0;
    double ry = 0.8;
    double rz = 1.0;

    fluidSim.setBounds(minx, maxx, miny, maxy, minz, maxz);
    std::vector<glm::vec3> points;
    double n = t["numParticles"].cast<int>();
    for (int i=0; i<n; i++) {
        float x = minx + rx*((float)rand()/RAND_MAX) * (maxx - minx);
        float y = miny + ry*((float)rand()/RAND_MAX) * (maxy - miny);
        float z = minz + rz*((float)rand()/RAND_MAX) * (maxz - minz);
        glm::vec3 p = glm::vec3(x, y, z);
        points.push_back(p);
    }
    fluidSim.addFluidParticles(points);

    points = utils::createPointPanel(2, 8, 0.7*radius, 3,
                                     glm::vec3(0.0, 0.0, 1.0), glm::vec3(0.0, 1.0, 0.0), true);
    //fluidSim.addObstacleParticles(points);
    fluidSim.translateObstacle(0, glm::vec3(maxx, 0.5*(maxy-miny), 0.5*(maxz-minz)));

    double damp = t["initialDampingConstant"].cast<double>();
    fluidSim.setDampingConstant(damp);

    fluidGradient = Gradients::getSkyblueGradient();
    minColorDensity = t["minColorDensity"].cast<double>();
    maxColorDensity = t["maxColorDensity"].cast<double>();
}

void GLWidget::activateSimulation() {
    using namespace luabridge;

    lua_State* L = luaL_newstate();
    luaL_openlibs(L);

    if (luaL_dofile(L, "scripts/fluid_config.lua") != 0) {
        qDebug() << "Error loading script";
        exit(1);
    }
    LuaRef t = getGlobal(L, "settings");
    LuaRef bounds = t["finalBounds"];
    double minx = bounds["minx"].cast<double>();
    double maxx = bounds["maxx"].cast<double>();
    double miny = bounds["miny"].cast<double>();
    double maxy = bounds["maxy"].cast<double>();
    double minz = bounds["minz"].cast<double>();
    double maxz = bounds["maxz"].cast<double>();
    fluidSim.setBounds(minx, maxx, miny, maxy, minz, maxz);

    double damp = t["finalDampingConstant"].cast<double>();
    fluidSim.setDampingConstant(damp);
    isRendering = true;
}

void GLWidget::stopSimulation() {
}

void GLWidget::updateSimulationSettings() {
    using namespace luabridge;

    lua_State* L = luaL_newstate();
    luaL_openlibs(L);

    if (luaL_dofile(L, "scripts/fluid_config.lua") != 0) {
        qDebug() << "Error loading script";
        exit(1);
    }
    LuaRef t = getGlobal(L, "settings");

    minColorDensity = t["minColorDensity"].cast<double>();
    maxColorDensity = t["maxColorDensity"].cast<double>();
    simulationFPS = t["fps"].cast<double>();
    isSimulationPaused = t["isSimulationPaused"].cast<bool>();
}

void GLWidget::updateSimulation() {
    // find delta time
    float dt = (float) deltaTimer->elapsed() / 1000;
    deltaTimer->restart();
    updateCameraMovement(dt);
    fluidSim.setCameraPosition(camera.position);

    dt = 1.0 / simulationFPS;
    if (isSimulationPaused) {
        dt = 0.0;
    }
    runningTime = runningTime + dt;


    // fluidSimulation
    using namespace luabridge;

    updateSimulationSettings();
    lua_State* L = luaL_newstate();
    luaL_openlibs(L);

    if (luaL_dofile(L, "scripts/fluid_config.lua") != 0) {
        qDebug() << "Error loading script";
        exit(1);
    }
    LuaRef t = getGlobal(L, "settings");
    bool isRenderingEnabled = t["isRenderingEnabled"].cast<bool>();

    if (!isSimulationPaused) {
        fluidSim.rotateObstacle(0, Quaternion(glm::vec3(0.0, 1.0, 0.0), 0.05));
        fluidSim.update(dt);
    }


    if (isRendering && isRenderingEnabled && !isSimulationPaused) {
        writeFrame();
    } else {
        updateGL();
    }
}

void GLWidget::writeFrame() {
    updateGL();

    std::string s = std::to_string(currentFrame);
    s.insert(s.begin(), 6 - s.size(), '0');
    s = "test_render/" + s + ".png";
    bool r = saveFrameToFile(QString::fromStdString(s));

    qDebug() << r << QString::fromStdString(s);

    currentFrame += 1;
}


// drawn when animation is running
void GLWidget::drawAnimation() {

}

bool GLWidget::saveFrameToFile(QString fileName) {
    GLubyte *data = (GLubyte*)malloc(4*(int)screenWidth*(int)screenHeight);
    if( data ) {
        glReadPixels(0, 0, screenWidth, screenHeight,
                     GL_RGBA, GL_UNSIGNED_BYTE, data);
    }

    QImage image((int)screenWidth, (int)screenHeight, QImage::Format_RGB32);
    for (int j=0; j<screenHeight; j++) {
        for (int i=0; i<screenWidth; i++) {
            int idx = 4*(j*screenWidth + i);
            char r = data[idx+0];
            char g = data[idx+1];
            char b = data[idx+2];

            // sets 32 bit pixel at (x,y) to yellow.
            //uchar *p = image.scanLine(j) + i;
            //*p = qRgb(255, 0, 0);
            QRgb value = qRgb(r, g, b);
            image.setPixel(i, screenHeight-j-1, value);
        }
    }
    bool r = image.save(fileName);
    free(data);

    return r;
}

void GLWidget::drawBillboard(GLuint *tex, glm::vec3 p, float width) {
    float hw = 0.5*width;
    glm::vec3 look = glm::normalize(camera.position - p);
    glm::vec3 right = glm::normalize(glm::cross(camera.up, look));
    glm::vec3 up = glm::cross(look, right);

    glm::mat4 mat = glm::transpose(glm::mat4(right.x, up.x, look.x, p.x,
                                             right.y, up.y, look.y, p.y,
                                             right.z, up.z, look.z, p.z,
                                             0.0, 0.0, 0.0, 1.0));
    glPushMatrix();
    glMultMatrixf((GLfloat*)&mat);

    glEnable(GL_TEXTURE_2D);
    glDisable(GL_LIGHTING);
    glBindTexture(GL_TEXTURE_2D, tex[0]);
    glBegin(GL_QUADS);
    glTexCoord2f(0.0f, 0.0f); glVertex3f(-hw, -hw,  0.0f);
    glTexCoord2f(1.0f, 0.0f); glVertex3f( hw, -hw,  0.0f);
    glTexCoord2f(1.0f, 1.0f); glVertex3f( hw,  hw,  0.0f);
    glTexCoord2f(0.0f, 1.0f); glVertex3f(-hw,  hw,  0.0f);
    glEnd();
    glEnable(GL_LIGHTING);
    glDisable(GL_TEXTURE_2D);

    glPopMatrix();
}

bool compareByZDistance(const SPHParticle *p1, const SPHParticle *p2) {
    return p1->zdistance > p2->zdistance;
}

void GLWidget::paintGL()
{
    camera.set();
    utils::drawGrid();

    float scale = 2.0;
    glm::mat4 scaleMat = glm::transpose(glm::mat4(scale, 0.0, 0.0, 0.0,
                                                  0.0, scale, 0.0, 0.0,
                                                  0.0, 0.0, scale, 0.0,
                                                  0.0, 0.0, 0.0, 1.0));
    glPushMatrix();
    glMultMatrixf((GLfloat*)&scaleMat);
    fluidSim.draw();

    // draw fluid particles. Color by density value.
    std::vector<SPHParticle*> particles = fluidSim.getAllParticles();
    std::sort(particles.begin(), particles.end(), compareByZDistance);

    SPHParticle *sp;
    float size = 0.5*fluidSim.getParticleSize();
    for (uint i=0; i<particles.size(); i++) {
        sp = particles[i];
        glm::vec3 p = particles[i]->position;

        if (particles[i]->isObstacle) {
            glColor3f(0.5, 0.5, 0.5);
        } else {
            glColor3f(sp->color.x, sp->color.y, sp->color.z);
        }

        if (particles[i]->isVisible) {
            drawBillboard(&texture[0], p, size);
        }
    }

    // draw obstacle particles
    particles = fluidSim.getObstacleParticles();
    glColor3f(0.5, 0.5, 0.5);
    for (uint i=0; i<particles.size(); i++) {
        if (particles[i]->isVisible) {
            glm::vec3 p = particles[i]->position;
            drawBillboard(&texture[0], p, size);
        }
    }

    glColor3f(0.0, 0.0, 0.0);
    glPointSize(6.0);

    glPopMatrix();



    camera.unset();
}

void GLWidget::resizeGL(int width, int height)
{
    screenWidth = width;
    screenHeight = height;
    camera.resize((float)width, (float)height);
    updateGL();
}

void GLWidget::keyPressEvent(QKeyEvent *event)
{
    // turn on camera movement
    isMovingForward = event->key()  == Qt::Key_W;
    isMovingBackward = event->key() == Qt::Key_S;
    isMovingRight = event->key()    == Qt::Key_D;
    isMovingLeft = event->key()     == Qt::Key_A;
    isMovingUp = event->key()       == Qt::Key_T;
    isMovingDown = event->key()     == Qt::Key_G;

    isRotatingRight = event->key()  == Qt::Key_E;
    isRotatingLeft = event->key()   == Qt::Key_Q;
    isRotatingUp = event->key()     == Qt::Key_F;
    isRotatingDown = event->key()   == Qt::Key_R;
    isRollingRight = event->key()   == Qt::Key_X;
    isRollingLeft = event->key()    == Qt::Key_Z;

    // slow down simulation
    if (event->key() == Qt::Key_C) {
        deltaTimeModifier = minDeltaTimeModifier;
    }
}

void GLWidget::keyReleaseEvent(QKeyEvent *event)
{
    // turn off camera movement
    isMovingForward = event->key()  == Qt::Key_W ? false : isMovingForward;
    isMovingBackward = event->key() == Qt::Key_S ? false : isMovingBackward;
    isMovingRight = event->key()    == Qt::Key_D ? false : isMovingRight;
    isMovingLeft = event->key()     == Qt::Key_A ? false : isMovingLeft;
    isMovingUp = event->key()       == Qt::Key_T ? false : isMovingUp;
    isMovingDown = event->key()     == Qt::Key_G ? false : isMovingDown;

    isRotatingRight = event->key()  == Qt::Key_E ? false : isRotatingRight;
    isRotatingLeft = event->key()   == Qt::Key_Q ? false : isRotatingLeft;
    isRotatingDown = event->key()   == Qt::Key_R ? false : isRotatingDown;
    isRotatingUp = event->key()     == Qt::Key_F ? false : isRotatingUp;
    isRollingRight = event->key()   == Qt::Key_X ? false : isRollingRight;
    isRollingLeft = event->key()    == Qt::Key_Z ? false : isRollingLeft;

    if (event->key() == Qt::Key_C) {
        deltaTimeModifier = maxDeltaTimeModifier;
    }

    if (event->key() == Qt::Key_H) {
        activateSimulation();
    }

}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
    (void)event;
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{
    (void)event;
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
    (void)event;
}








