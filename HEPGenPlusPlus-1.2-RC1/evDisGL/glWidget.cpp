#include "glWidget.h"
#include "shaders.h"
glWidget::glWidget(myWindow* _windowPointer, QWidget* _parent,physicsEngine* _inEngine) :  QGLWidget(QGLFormat(QGL::SampleBuffers), _parent)
{
    theEngine = _inEngine;
    windowPointer = _windowPointer;
}


void glWidget::initializeGL()
{
    qDebug() << "init" << endl;
    glEnable( GL_DEPTH_TEST);
    glEnable( GL_CULL_FACE);
    glEnable( GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);

    myList = theEngine->getNextEvent();
    
    glLineWidth(2.f);

    rotx = 20.0f;
    roty = 40;
    zZoom = -3.6f;

    rotTimer = new QTimer(this);
    connect (rotTimer, SIGNAL(timeout()),this, SLOT(rotateView()));

    eventTimer = new QTimer(this);
    connect (eventTimer, SIGNAL(timeout()),this, SLOT(nextEventView()));

    printf("/******---- Loading and compiling shaders\n");
    program = new QGLShaderProgram(context());
    program->addShaderFromSourceCode(QGLShader::Vertex,vertexShaderSource);
    program->addShaderFromSourceCode(QGLShader::Fragment,fragmentShaderSource);
    printf("/******---- Linking and binding shaders\n");
    program->link();
    program->bind();
    printf("/******---- Getting vertex,matrix and color locations from shaders\n");

    vertexLocation = program->attributeLocation("vertex");
    matrixLocation = program->uniformLocation("matrix");
    screenDimensionLocation = program->uniformLocation("screenDim");
    colorLocation = program->attributeLocation("color");
    worldMatrixLocation = program->attributeLocation("worldMatrix");
    
    

}

void glWidget::drawCoordinates()
{
//     glPushMatrix();
//     glTranslatef(-.0002f,0.0f,0.0f);

    QMatrix4x4 model = modelStack.front();
    model.translate(-0.0002,0.0,0.0);
    
    program->setUniformValue(worldMatrixLocation,model);
    static const float axisArray[18] = {-4.0,0.0,0.0,4.0,0.0,0.0,
                                        0.0,-4.0,0.0,0.0,4.0,0.0,
                                        0.0,0.0,-4.0,0.0,0.0,4.0,
                                       };

    static const float colourArray[18] = {1.0,0.5,0.0,0.0,1.0,1.0,
                                          1.0,0.5,0.0,0.0,1.0,1.0,
                                          1.0,0.5,0.0,0.0,1.0,1.0
                                         };
    QColor color(0, 255, 0, 255);


    program->setAttributeArray(vertexLocation,axisArray,3);
    QMatrix4x4 pmvMatrix = projectionStack.front() * viewStack.front() * model;
    program->setUniformValue(matrixLocation,pmvMatrix);
    program->setAttributeArray(colorLocation,colourArray,3);
    for (int i = 0; i < 3; i++)
        glDrawArrays(GL_LINES,2*i,2);
    
    program->setUniformValue(worldMatrixLocation,modelStack.front());
}

void glWidget::drawGrid()
{
    const int nLines = 6;
    float axisArray[(4*(nLines)-2)*6*2];
    float colourArray[(4*(nLines)-2)*6*2];
    fill(colourArray,colourArray+(4*(nLines)-2)*6*2,192./255.);
//     fill(colourArray,colourArray+(4*(nLines)-2)*6*2,0.f);
    
    
    memset(axisArray,0,sizeof(axisArray));
    int lineNum =0;
    for (int i = -nLines+1; i < nLines; i ++) {
        if (i == 0)
            continue;
        float varyComp = i*(4.f/nLines);
        //lines for x=const and varying z
        axisArray[lineNum*12]=-4.f;
        axisArray[lineNum*12+1]=0.f;
        axisArray[lineNum*12+2]=varyComp;
        axisArray[lineNum*12+3]=+4.f;
        axisArray[lineNum*12+4]=0.f;
        axisArray[lineNum*12+5]=varyComp;

        axisArray[lineNum*12+6]=varyComp;
        axisArray[lineNum*12+7]=0.f;
        axisArray[lineNum*12+8]=-4.f;
        axisArray[lineNum*12+9]=varyComp;
        axisArray[lineNum*12+10]=0.f;
        axisArray[lineNum*12+11]=+4.f;
        lineNum++;
    }
    //no translation here
    QMatrix4x4 pmvMatrix = projectionStack.front() * viewStack.front() * modelStack.front();
    program->setAttributeArray(vertexLocation,axisArray,3);
    program->setAttributeArray(colorLocation,colourArray,3);
    for (int i =0; i < 4*(nLines)-2; i ++)
        glDrawArrays(GL_LINES,i*2,2);

}


void glWidget::drawEvent()
{
    float vertizes[myList.size()*6];
    float vertexColours[myList.size()*6];

    memset(vertizes,0,sizeof(float)*myList.size()*6);
    memset(vertexColours,0,sizeof(float)*myList.size()*6);
    //rotate z to x
    for (int i =0; i < myList.size(); i++)   {
        vertizes[i*6]= myList.at(i).posVertex1[0];
        vertizes[i*6+1]= myList.at(i).posVertex1[1];
        vertizes[i*6+2]= myList.at(i).posVertex1[2];
        vertizes[i*6+3+0]= myList.at(i).posVertex2[0];
        vertizes[i*6+3+1]= myList.at(i).posVertex2[1];
        vertizes[i*6+3+2]= myList.at(i).posVertex2[2];
        memcpy(&vertexColours[i*6],&myList.at(i).colorVertex[0],sizeof(float)*3);
        memcpy(&vertexColours[i*6+3],&myList.at(i).colorVertex[0],sizeof(float)*3);
    }
    program->setAttributeArray(vertexLocation,vertizes,3);
    program->setAttributeArray(colorLocation,vertexColours,3);

    for (int i =0; i < myList.size(); i ++)
        glDrawArrays(GL_LINES,i*2,2);
}



void glWidget::paintGL(bool _quiet)
{

    program->enableAttributeArray(vertexLocation);
    program->enableAttributeArray(colorLocation);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

      glClearColor(0.0,0.0,0.0,1.0);
//      glClearColor(1.0,1.0,1.0,1.0);

    //make the modelmatrix ready
    QMatrix4x4 modelTransform;
    modelTransform.translate(0.0,0.0,zZoom);
    modelTransform.rotate(rotx,QVector3D(1.0,0.0,0.0));
    modelTransform.rotate(roty,QVector3D(0.0,1.0,0.0));
    modelStack.push_front(modelTransform);
    program->setUniformValue(worldMatrixLocation,modelStack.front());
    drawCoordinates();
    drawGrid();
    drawEvent();

    if (!_quiet)
        printf("rotX: %.3f rotY: %.3f zZoom:%.3f\n",rotx,roty,zZoom);

    //throw it out again - actually we are wasting time here since there is only this one model
    //but it does not really matter, since we are not using much gpu time anyway
    //and its good practice to clean up the matrix-stack after using it
    modelStack.pop_front();
    swapBuffers();
    program->disableAttributeArray(vertexLocation);
    program->disableAttributeArray(colorLocation);

}

void glWidget::wheelEvent(QWheelEvent* event)
{
    int direction = event->delta() / abs(event->delta());
    zZoom += 1.2*direction;
    paintGL();
    QWidget::wheelEvent(event);
}


void glWidget::mousePressEvent(QMouseEvent* event)
{
    setFocus();
    if (event->button() == Qt::LeftButton) {
        if (rotTimer->isActive())
            rotTimer->stop();
        else
            rotTimer->start(5);
    }
    else if (event->button() == Qt::RightButton)
    {
        if (eventTimer->isActive())
            eventTimer->stop();
        else
            eventTimer->start(750);
    }
    QWidget::mousePressEvent(event);
}



void glWidget::keyPressEvent(QKeyEvent* e)
{

    double angleStep = 1.2;
    double zStep = 0.6;
    switch( e->key() ) {
    case Qt::Key_Left:
        roty += angleStep;
        break;
    case Qt::Key_Right:
        roty -= angleStep;
        break;
    case Qt::Key_Up:
        rotx += angleStep;
        break;
    case Qt::Key_Down:
        rotx -= angleStep;
        break;
    case Qt::Key_PageUp:
        zZoom += zStep;
        break;
    case Qt::Key_PageDown:
        zZoom -= zStep;
        break;
    case Qt::Key_Space:
        myList = theEngine->getNextEvent();
        break;
    case Qt::Key_Escape:
        QApplication::quit();
        break;
    case Qt::Key_S:
        theEngine->triggerTheSave();
        break;
    default:
        QWidget::keyPressEvent(e);
    }
    paintGL();
}

void glWidget::nextEventView()
{
    myList = theEngine->getNextEvent();
    paintGL();
}


void glWidget::resizeGL(int w, int h)
{
    printf("resizeGL\n");
    if (h <=0)
        h=1;
    
    screenDimensions.setX(w);
    screenDimensions.setY(h);
    
    program->setUniformValue(screenDimensionLocation,screenDimensions);
    
    QMatrix4x4 viewPort;
    viewPort.lookAt(QVector3D(0.0,0.0,0.0),QVector3D(0.0,0.0,-6.0),QVector3D(0.0,1.0,0.0));
    if (viewStack.size() > 0)
        viewStack.pop_front();
    viewStack.push_front(viewPort);
    

//     // Establish the clipping volume by setting up an perspective projection
    double aspectRatio = double( w ) / double( h );
    glViewport(0,0,w,h);
    double verticalViewingAngle = 45.0f;
    QMatrix4x4 projection;
    projection.perspective(verticalViewingAngle,aspectRatio,1.0,100.0);
    if (projectionStack.size() > 0)
        projectionStack.pop_front();
    projectionStack.push_front(projection);
    paintGL();
}







