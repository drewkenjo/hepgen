
#ifndef GLWIN_H
#define GLWIN_H

#include <qt4/QtGui/QWidget>
#include <qt4/QtGui/QHBoxLayout>

#include "glWidget.h"
#include "physicsEngine.h"

class myWindow: public QWidget
{
    Q_OBJECT
public:
    myWindow(physicsEngine* _inEngine);



private:
    physicsEngine* theEngine;
    glWidget* myGLWidget;


};


#endif

