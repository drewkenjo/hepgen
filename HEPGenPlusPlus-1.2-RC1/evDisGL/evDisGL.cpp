/*! std c includes */
#include <iostream>
#include <ctime>
#include <string>

/*! qt stuff */
#include <qt4/QtGui/QApplication>
#include "dialog.h"

/*! physics stuff */
#include "physicsEngineBase.h"
#include "physicsEngine.h"
#include "physicsEngineLepto.h"
#include "hparammanager.h"

const string currentDateTime()
{
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
    return buf;
}

using namespace std;

int main( int argc, char** argv)
{
    if (argc < 4) {
        printf("For the evDisGL please use the following parameters:\n");
	printf("1.) Libhepgen.so use: \n");
	printf("  ./%s -hepgen [DATACARD] [SEED]\n",argv[0]);
	printf("2.)LEPTOFILE use: \n");
	printf("  ./%s -leptofile [FILE_IN] [FILE_OUT]\n",argv[0]);
	exit(0);
    }
    physicsEngine* myEngine;
    if (strcmp(argv[1],"-hepgen")==0)
      myEngine = new physicsEngineHEPGen;
    else if (strcmp(argv[1],"-leptofile")==0)
      myEngine = new physicsEngineLEPTO;
    else{
	printf("Please read the allowed operation modes!\n");
	return -1;
      }
	
    myEngine->setSettingsAndSeed(argc,argv);	
    QApplication myApp(argc,argv);
    myApp.setApplicationName("evDisGL");
    myWindow* window = new myWindow(myEngine);
    window->show();
    int r = myApp.exec();
    printf("Closing output files and saving histograms!\n");
    myEngine->closeFiles();
    return r;
}
