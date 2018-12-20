/*!
 *  \file hphigen.h
 *  \date Created on: Feb 26, 2014
 *  \author Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2014 All Right Reserved
 */


#define _USE_MATH_DEFINES

#include "hphysicsgen.h"
#include "hconstants.h"
#include "hhelper.h"
#include "hevent.h"


#ifndef HPHIGEN_H_
#define HPHIGEN_H_




/*!
 * \brief Implementation of a Phi-Generator
 *
 */
class HPhysicsGenPHI : public HPhysicsGen
{


public:
    HPhysicsGenPHI() {};
    HPhysicsGenPHI(CLHEP::HepRandom* _randGen, HParamManager* _paramMan, HEvent* _event, HBooker* _booker, HBooker* _ddBooker);
    ~HPhysicsGenPHI();

    void generateListOfParticles();

    void generateEvent();


    void addHistograms();

    void generateDecay();

    bool generateMesonMass();

    void calcWeights();
    static void calcWeights(hWeightInterface* myInt,double& WEIGHTRET);

    void generatePolarisation();

    void setUserVars();
    void setParameters();

    int polarization;

    double weight;
    double thetapi;
    double phipi;
    double x,y,z;

    int eventcount;
    int resetcount;


};



























#endif





