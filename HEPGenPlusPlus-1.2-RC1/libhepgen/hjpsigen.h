/*!
 *  \file hrhogen.h
 *  \date Created on: Mar 2, 2013
 *  \author Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2013 All Right Reserved
 */


#define _USE_MATH_DEFINES

#include "hphysicsgen.h"
#include "hconstants.h"
#include "hhelper.h"
#include "hevent.h"


#ifndef HJPSIGEN_H_
#define HJPSIGEN_H_




/*!
 * \brief Implementation of a J/Psi-Generator
 *
 */
class HPhysicsGenJPsi : public HPhysicsGen
{


public:
    HPhysicsGenJPsi() {};
    HPhysicsGenJPsi(CLHEP::HepRandom* _randGen, HParamManager* _paramMan, HEvent* _event, HBooker* _booker, HBooker* _ddBooker);
    ~HPhysicsGenJPsi();

    void generateListOfParticles();

    void generateEvent();


    void addHistograms();

    void generateDecay();
    
    bool generateMesonMass();

    void calcWeights();
    static void calcWeights(hWeightInterface* myInt,double& WEIGHTRET);

    void generatePolarisation();
    
    void printAuxVars();

    void setUserVars();
    void setParameters();

    int polarization;

    double weight;
    double thetapi;
    double phipi;
    double x,y,z;
    
    int auxFlag;

    int eventcount;
    int resetcount;

    bool dec_muons;

};



























#endif





