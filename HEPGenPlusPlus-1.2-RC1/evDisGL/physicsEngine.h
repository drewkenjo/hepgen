#ifndef PHYSICS_ENGINE
#define PHYSICS_ENGINE

#include "physicsEngineBase.h"

/* libhepgen  */
#include "hvector.h"
#include "config.h"
#include "hlorentzvector.h"
#include "hcardparser.h"
#include "hparammanager.h"
#include "hpionicdata.h"
#include "hbooker.h"
#include "hbookbackendASCII.h"
#include "hgenmanager.h"

/* libhepgen ROOT addons */
#ifdef USE_ROOT
#include "hbookbackendROOT.h"
#endif



/* std includes */
#include <iostream>
#include <string>
#include <cstring>
#include <string.h>
#include <time.h>

/* qt includes */
#include <qt4/QtCore/QString>
#include <qt4/QtCore/QDebug>





class physicsEngineHEPGen : public physicsEngine
{

public:
    std::vector<line> getNextEvent();

    void triggerTheSave(){tmp->saveLastEvent();printf("Event saved!\n");}
    void closeFiles(){tmp->closeOutPut();}

    void setSettingsAndSeed(int argc, char** argv);



    physicsEngineHEPGen();

    ~physicsEngineHEPGen();

private:
    HGenManager* tmp;
    const string currentDateTime()
    {
        time_t     now = time(0);
        struct tm  tstruct;
        char       buf[80];
        tstruct = *localtime(&now);
        strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
        return buf;
    }
};


#endif
