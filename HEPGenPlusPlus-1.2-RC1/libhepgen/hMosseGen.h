/*!
 *  \file hMosseGen.h
 *  \date Created on: Mar 16, 2015
 *  \author: Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2015 All Right Reserved
 */


#define _USE_MATH_DEFINES

#include "hphysicsgen.h"
#include "hconstants.h"
#include "hhelper.h"
#include "hevent.h"
#include "myTHEO.hh"
#include "GPDQ.hh"



#ifndef HDVCSGENMOSSE_H_
#define HDVCSGENMOSSE_H_




/*!
 * \brief Implementation of a DVCS-Generator with the weight calculation from L. Mosse
 *
 */
class HPhysicsGenDVCSMosse : public HPhysicsGen
{
public:
    HPhysicsGenDVCSMosse() {};
    HPhysicsGenDVCSMosse(CLHEP::HepRandom* _randGen, HParamManager* _paramMan, HEvent* _event, HBooker* _booker, HBooker* _ddBooker);
    ~HPhysicsGenDVCSMosse();
    void generateListOfParticles();



    void generateEvent();

    void calcWeights();
    static void calcWeights(hWeightInterface _in, double& BH, double &DVCS, double &INT);
    bool generateMesonMass();

    void addHistograms();

    void setUserVars();
    void setParameters();

    static void cinematique(C_Lorentz_Vector*,C_Lorentz_Vector*,C_Lorentz_Vector*,C_Lorentz_Vector*,C_Lorentz_Vector*,C_Lorentz_Vector*,C_Lorentz_Vector*,Complexe ,Complexe ,Complexe ,Complexe ,Complexe ,Complexe );
    static void cinematiqueOZ(C_Lorentz_Vector*,C_Lorentz_Vector*,C_Lorentz_Vector*,C_Lorentz_Vector*,C_Lorentz_Vector*,C_Lorentz_Vector*,C_Lorentz_Vector*,Complexe ,Complexe ,Complexe ,Complexe ,Complexe ,Complexe );
    static void cinematiqueOZT(C_Lorentz_Vector*,C_Lorentz_Vector*,C_Lorentz_Vector*,C_Lorentz_Vector*,C_Lorentz_Vector*,C_Lorentz_Vector*,C_Lorentz_Vector*,Complexe ,Complexe ,Complexe ,Complexe ,Complexe ,Complexe );


private:



    double phir_trento;
    double weightdvcs, weightbh, weightint;
    double weightresult0,weightresult1,weightresult2,weight ;
    double none;

};


#endif



