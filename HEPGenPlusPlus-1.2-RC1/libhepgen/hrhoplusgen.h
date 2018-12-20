/*!
 *  \file hrhoplusgen.h
 *  \date Created on: Jul 21, 2014
 *  \author Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2014 All Right Reserved
 */


#define _USE_MATH_DEFINES

#include "hphysicsgen.h"
#include "hconstants.h"
#include "hhelper.h"
#include "hevent.h"


#ifndef HRHOGENPLUS_H_
#define HRHOGENPLUS_H_




/*!
 * \brief Implementation of a Rho0-Generator
 *
 */
class HPhysicsGenRHOPlus : public HPhysicsGen
{


public:
    HPhysicsGenRHOPlus() {};
    HPhysicsGenRHOPlus(CLHEP::HepRandom* _randGen, HParamManager* _paramMan, HEvent* _event, HBooker* _booker, HBooker* _ddBooker);
    ~HPhysicsGenRHOPlus();

    void generateListOfParticles();

    void generateEvent();


    void addHistograms();

    void generateDecay();

    bool generateMesonMass();

    void calcWeights();
    static void calcWeights(hWeightInterface _in, double &WEIGHT,HPionicData* dataTable = NULL);
    

    void generatePolarisation();

    void setUserVars();
    void setParameters();

    int polarization;

    void printAuxVars();

    double weight;
    double thetapi;
    double phipi;
    double x,y,z;

    int eventcount;
    int resetcount;

private:
    HParticle outPhoton1;
    HParticle outPhoton2;
    HParticle outPhoton1_Lab;
    HParticle outPhoton2_Lab;



};



























#endif





