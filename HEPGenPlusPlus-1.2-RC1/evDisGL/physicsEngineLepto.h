#ifndef PHYSICS_ENGINE_LEPTO
#define PHYSICS_ENGINE_LEPTO

#include "physicsEngineBase.h"

/* std includes */
#include <iostream>
#include <string>
#include <cstring>
#include <string.h>
#include <time.h>
#include <cassert>


/* lepto file reading */
#include "lfread.h"
#include "lfwrite.h"


class physicsEngineLEPTO : public physicsEngine
{
public:
  std::vector<line> getNextEvent();
  void triggerTheSave();
  void closeFiles();
  void setSettingsAndSeed(int argc, char** argv);

private:
  
    lfr::lEvent myEvent;
    lfr::lHeader myHeader;
    fstream inputFile;
    fstream outputFile;
    bool largeHeaders;
  
  
};


#endif

