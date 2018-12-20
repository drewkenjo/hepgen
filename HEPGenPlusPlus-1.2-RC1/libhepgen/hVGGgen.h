/*!
 *  \file hVGGgen.h
 *  \date Created on: Aug 1, 2014
 *  \author: Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2014 All Right Reserved
 */


#define _USE_MATH_DEFINES

#include "hphysicsgen.h"
#include "hconstants.h"
#include "hhelper.h"
#include "hevent.h"
#include "compassInterFace.h"

#ifndef HDVCSGENB_H_
#define HDVCSGENB_H_


/*!
 * \brief Implementation of a DVCS-Generator with the weight calculation from H. Moutarde
 *
 */
class HPhysicsGenDVCSVGG : public HPhysicsGen
{
public:
    HPhysicsGenDVCSVGG() {};
    HPhysicsGenDVCSVGG(CLHEP::HepRandom* _randGen, HParamManager* _paramMan, HEvent* _event, HBooker* _booker, HBooker* _ddBooker);
    ~HPhysicsGenDVCSVGG();
    void generateListOfParticles();

    static void calcWeights(hWeightInterface _in, double* BH, double* DVCS, double* INT,bool _wantLO=false);
    
    
    void generateEvent();

    void calcWeights();
    bool generateMesonMass();

    void addHistograms();

    void setUserVars();
    void setParameters();


private:


    double phir_trento;
    double weightdvcs, weightbh, weightint;
    double weightresult0,weightresult1,weightresult2,weight ;
    double none;

};


#endif


