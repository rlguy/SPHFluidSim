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

#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <GL/glu.h>
#include <QGLWidget>
#include <vector>
#include <tuple>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <QImage>
#include "glm/glm.hpp"
#include "camera3d.h"
#include "utils.h"
#include "spatialgrid.h"
#include "sphfluidsimulation.h"

class QtLogo;

/**
  * @class  GLWidget
  * @brief  our OpenGL view derived from QGLWidget.
  * We have to override several functions for our
  * application-specific OpenGL drawing functionality
  */
class GLWidget : public QGLWidget
{
    Q_OBJECT

public:
    GLWidget(QWidget *parent = 0);
    ~GLWidget();

    void keyPressEvent(QKeyEvent *event);
    void keyReleaseEvent(QKeyEvent *event);

    QSize minimumSizeHint() const;
    QSize sizeHint() const;

protected:
    void initializeGL();
    void resizeGL(int width, int height);
    void mousePressEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);

public slots:
    void paintGL();
    void updateSimulation();

private:
    float screenWidth;
    float screenHeight;

    void updatePerspective();
    void updateCameraMovement(float dt);
    void drawGrid();
    void drawAnimation();
    void initializeSimulation();
    void stopSimulation();
    bool saveFrameToFile(QString fileName);

    // update/draw tiemrs
    QTimer *drawTimer;
    QTimer *updateTimer;
    QTime *deltaTimer;

    // camera and movement
    camera3d camera;
    bool isMovingForward;
    bool isMovingBackward;
    bool isMovingRight;
    bool isMovingLeft;
    bool isMovingUp;
    bool isMovingDown;
    bool isRotatingRight;
    bool isRotatingLeft;
    bool isRotatingUp;
    bool isRotatingDown;
    bool isRollingRight;
    bool isRollingLeft;

    // simulation system
    float minDeltaTimeModifier;
    float maxDeltaTimeModifier;
    float deltaTimeModifier;
    float runningTime;
    int currentFrame;

    SPHFluidSimulation fluidSim;
};

#endif







