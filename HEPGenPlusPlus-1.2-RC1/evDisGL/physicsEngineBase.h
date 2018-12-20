#ifndef PHYSICS_ENGINE_BASE
#define PHYSICS_ENGINE_BASE


/* interchange structs */
#include "vStructs.h"

class physicsEngine{
public:
  virtual std::vector<line> getNextEvent() = 0;
  virtual void triggerTheSave() = 0;
  virtual void closeFiles() = 0;
  virtual void setSettingsAndSeed(int argc, char** argv) = 0;
};

#endif

