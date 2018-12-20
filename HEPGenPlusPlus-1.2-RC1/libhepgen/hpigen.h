/*!
 *  \file hpigen.h
 *  \date Created on: Mar 7, 2013
 *  \author Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2013 All Right Reserved
 */

#define _USE_MATH_DEFINES

#include "hphysicsgen.h"
#include "hconstants.h"
#include "hhelper.h"
#include "hevent.h"
#include "hpionicdata.h"
#include "hWeightingInterface.h"

#ifndef HPIGEN_H_
#define HPIGEN_H_

/*!
 * \brief Implementation of a Pi0-Generator
 *
 */
class HPhysicsGenPI : public HPhysicsGen
{


public:
    HPhysicsGenPI() : HPhysicsGen() {};
    HPhysicsGenPI(CLHEP::HepRandom* _randGen, HParamManager* _paramMan, HEvent* _event, HBooker* _booker, HBooker* _ddBooker);
    ~HPhysicsGenPI();

    void generateListOfParticles();

    void generateEvent();


    void generateDecay();

    bool generateMesonMass();

    void calcWeights();
    static void calcWeights(hWeightInterface _in, double &WEIGHT,HPionicData* dataTable = NULL);
    

    void printAuxVars();

    void addHistograms();

    void setUserVars();
    void setParameters();


    double weight;
    double thetapi;
    double phipi;
    double costhetapi;
private:
  double newWeightUnRounded,newWeightRounded,newWeightRoundedAbs;
  HPionicData* myPionic;
};



#endif

