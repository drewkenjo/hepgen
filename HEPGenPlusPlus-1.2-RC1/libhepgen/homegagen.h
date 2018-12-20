/*!
 *  \file homegagen.h
 *  \date Created on: July 27, 2014
 *  \author Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2014 All Right Reserved
 */


#define _USE_MATH_DEFINES

#include "hphysicsgen.h"
#include "hconstants.h"
#include "hhelper.h"
#include "hevent.h"
#include "hrotmat.h"

#ifndef HPhysicsGenOMEGA_H_
#define HPhysicsGenOMEGA_H_




/*!
 * \brief Implementation of a omega-Generator
 *
 */
class HPhysicsGenOMEGA : public HPhysicsGen
{


public:
    HPhysicsGenOMEGA() {};
    HPhysicsGenOMEGA(CLHEP::HepRandom* _randGen, HParamManager* _paramMan, HEvent* _event, HBooker* _booker, HBooker* _ddBooker);
    ~HPhysicsGenOMEGA();

    void generateListOfParticles();

    void generateEvent();


    void addHistograms();

    void generateDecay();

    bool generateMesonMass();

    void calcWeights();
    /*! \brief static version of calcweights for usage in the hepgen_in_phast or other features
     * Please note: This function calls hreweightKine::getFluxCompensator and therefore needs myInt to be set with:
     * Qsq,tprime, y,nu,slpin -!! it will not error if they are not set correctly!
     */
    static void calcWeights(hWeightInterface* myInt,double& WEIGHTRET);

    void generatePolarisation();

    void setUserVars();
    void setParameters();

    int polarization;

    double weight;
    double thetapi;
    double phipi;
    double x,y,z;
    
    double xD,yD;




    int eventcount;
    int resetcount;
private:
    /*! \brief changes the momentum components cyclish style */
    void changeVector(HVector3& _in);
    HRotMat rotationMatrix;
    vector<HParticle*> decayParticles;

};















#endif





