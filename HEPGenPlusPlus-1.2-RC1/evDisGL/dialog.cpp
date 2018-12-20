#include "dialog.h"

myWindow::myWindow(physicsEngine* _inEngine): QWidget()
{
    theEngine = _inEngine;
    myGLWidget = new glWidget(this,this,_inEngine);
    QHBoxLayout* layout = new QHBoxLayout;
    layout->addWidget(myGLWidget);
    setLayout(layout);
    resize(1024,768);


    setWindowTitle("HEPGen++ event display OpenGL");
    myGLWidget->setFocus();
}
