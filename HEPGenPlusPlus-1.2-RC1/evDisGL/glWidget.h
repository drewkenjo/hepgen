#ifndef GLWIDG_H
#define GLWIDG_H
#define GL_GLEXT_PROTOTYPES 1
#define GL3_PROTOTYPES 1

/*! physics stuff */
#include "vStructs.h"
#include "physicsEngine.h"

/*! qt stuff */
#include <qt4/QtOpenGL/QGLWidget>
#include <qt4/QtGui/QWidget>
#include <qt4/QtCore/QDebug>
#include <qt4/QtOpenGL/QGLBuffer>
#include <GL/glu.h>
#include <GL/gl.h>
#include <GLES3/gl3.h>
#include <QKeyEvent>
#include <qt4/QtOpenGL/QGLShaderProgram>
#include <qt4/QtOpenGL/QGLShader>
#include <QEvent>
#include <cassert>
#include <QApplication>
#include <QTimer>
#include <QVector>
#include <QMatrix4x4>

#include "physicsEngine.h"

class myWindow;

class glWidget : public QGLWidget
{
    Q_OBJECT
public:
    glWidget(myWindow* _windowPointer,QWidget* _parent,physicsEngine* _inEngine);

    QVector2D screenDimensions;
public slots:
  void rotateView(){roty += 0.2;paintGL(true);}
  void nextEventView();
    
private:
  
    physicsEngine* theEngine;
    QVector<QMatrix4x4> modelStack;
    QVector<QMatrix4x4> viewStack;
    QVector<QMatrix4x4> projectionStack;
  
    QTimer* rotTimer;  
    QTimer* eventTimer;
    GLuint positionBuffer;
    GLuint gridBuffer;
    GLuint colourBuffer;
    GLuint vao;
    GLuint axisBuffer;
    
    void initializeGL();
    void resizeGL(int w, int h);

    void paintGL(bool _quiet = false);
    void drawEvent();
    void drawGrid();
    void drawCoordinates();
    myWindow* windowPointer;
    std::vector<line> myList;
    float rotx,roty,zZoom;

    void keyPressEvent( QKeyEvent *e );

    void mousePressEvent ( QMouseEvent * event );
    void wheelEvent(QWheelEvent *event);
    
    
    
    QGLShaderProgram* program;
    int vertexLocation;
    int matrixLocation;
    int colorLocation;
    int worldMatrixLocation;
    int screenDimensionLocation;
    

};



#endif


