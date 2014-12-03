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
#include "glwidget.h"
#include "glm/glm.hpp"
#include "quaternion.h"

#ifndef GL_MULTISAMPLE
#define GL_MULTISAMPLE  0x809D
#endif

//! [0]
GLWidget::GLWidget(QWidget *parent)
    : QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
{

    screenWidth = 1200;
    screenHeight = 600;

    //initializealize update/draw timers
    float updatesPerSecond = 30;
    float drawsPerSecond = 30;

    drawTimer = new QTimer(this);
    connect(drawTimer, SIGNAL(timeout()), this, SLOT(updateGL()));
    drawTimer->start(1000.0/drawsPerSecond);

    updateTimer = new QTimer(this);
    connect(updateTimer, SIGNAL(timeout()), this, SLOT(updateSimulation()));
    updateTimer->start(1000.0/updatesPerSecond);

    deltaTimer = new QTime();
    deltaTimer->start();

    // Initialize camera
    glm::vec3 pos = glm::vec3(16.88, 11.82, 11.41);
    glm::vec3 dir = glm::vec3(-0.828, -0.472, -0.302);
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

    grid = SpatialGrid(1.0);

    srand(time(NULL));
    float max = 10.0;

    int n = 100000;
    for (int i=0; i<n; i++) {
        float r = (float)rand() / (float)RAND_MAX;
        float x = r * max;
        r = (float)rand() / (float)RAND_MAX;
        float y = r * max * cos(x);
        r = (float)rand() / (float)RAND_MAX;
        float z = r * max * sin(y)*cos(x);

        float vmax = 0.1;
        float vx = (2*((float)rand() / (float)RAND_MAX) - 1.0) * vmax;
        float vy = (2*((float)rand() / (float)RAND_MAX) - 1.0) * vmax;
        float vz = (2*((float)rand() / (float)RAND_MAX) - 1.0) * vmax;

        glm::vec3 p = glm::vec3(x, y, z);
        glm::vec3 v = glm::vec3(vx, vy, vz);
        int ref = grid.insertPoint(p);
        gridTest.push_back(std::tuple<int,glm::vec3, glm::vec3>(ref, p, v));
    }

    //grid.movePoint(ref, glm::vec3(5.0, 1.0, 1.0));

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
}

void GLWidget::stopSimulation() {
}

void GLWidget::updateSimulation() {
    // find delta time
    float dt = (float) deltaTimer->elapsed() / 1000;
    runningTime = runningTime + dt;
    deltaTimer->restart();
    updateCameraMovement(dt);

    dt *= deltaTimeModifier;  // speed of simulation

    dt = 1.0/30.0;
    //grid test
    for (int i=0; i<(int)gridTest.size(); i++) {
        int ref = std::get<0>(gridTest[i]);
        glm::vec3 p = std::get<1>(gridTest[i]);
        glm::vec3 v = std::get<2>(gridTest[i]);

        p = p + dt*v;
        std::get<1>(gridTest[i]) = p;

        if ((float)rand() / (float)RAND_MAX > 0.98) {
            float vmax = 0.5;
            float vx = (2*((float)rand() / (float)RAND_MAX) - 1.0) * vmax;
            float vy = (2*((float)rand() / (float)RAND_MAX) - 1.0) * vmax;
            float vz = (2*((float)rand() / (float)RAND_MAX) - 1.0) * vmax;
            std::get<2>(gridTest[i]) = glm::vec3(vx, vy, vz);
        }

        grid.movePoint(ref, p);
    }
}

// Draws coordinate axis' and floor grid
void GLWidget::drawGrid() {
    // draw axis'
    float len = 10.0;
    glLineWidth(3.0);
    glBegin(GL_LINES);
    glColor3f(1.0, 0.0, 0.0);   // x
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(len, 0.0, 0.0);
    glColor3f(0.0, 1.0, 0.0);   // y
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, len, 0.0);
    glColor3f(0.0, 0.0, 1.0);   // z
    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, 0.0, len);
    glEnd();

    // draw outline around xy, zy planes
    glLineWidth(2.0);
    glColor4f(0.0, 0.0, 0.0, 0.3);
    glBegin(GL_LINES);
    glVertex3f(0.0, len, 0.0);
    glVertex3f(len, len, 0.0);
    glVertex3f(len, len, 0.0);
    glVertex3f(len, 0.0, 0.0);
    glVertex3f(0.0, len, 0.0);
    glVertex3f(0.0, len, len);
    glVertex3f(0.0, len, len);
    glVertex3f(0.0, 0.0, len);
    glEnd();


    // draw xz plane grid
    float spacing = 0.25;
    int yLines = 120;
    int zLines = 120;
    float height = (float)yLines * spacing;
    float width = (float)zLines * spacing;

    float z = spacing;
    glLineWidth(1.0);
    glColor4f(0.0, 0.0, 0.0, 0.2);
    glBegin(GL_LINES);
    for (int i=0; i < yLines; i++) {
        glVertex3f(0.0, 0.0, z);
        glVertex3f(width, 0.0, z);
        z += spacing;
    }

    float x = spacing;
    for (int i=0; i < zLines; i++) {
        glVertex3f(x, 0.0, 0.0);
        glVertex3f(x, 0.0, height);
        x += spacing;
    }
    glEnd();

}


// drawn when animation is running
void GLWidget::drawAnimation() {

}

void GLWidget::paintGL()
{
    camera.set();
    drawGrid();

    float scale = 2.0;
    glm::mat4 scaleMat = glm::transpose(glm::mat4(scale, 0.0, 0.0, 0.0,
                                                  0.0, scale, 0.0, 0.0,
                                                  0.0, 0.0, scale, 0.0,
                                                  0.0, 0.0, 0.0, 1.0));
    glPushMatrix();
    glMultMatrixf((GLfloat*)&scaleMat);
    grid.draw();
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








