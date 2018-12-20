/*!
 *  \file hphysicsgen.h
 *  \date Created on: Feb 6, 2013
 *  \author: Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2013 All Right Reserved
 */


#define _USE_MATH_DEFINES


#ifndef HPHYSGENa_H_
#define HPHYSGENa_H_

#include <cmath>
#include <CLHEP_EMBEDDED/Random/Random/Random.h>
#include <CLHEP_EMBEDDED/Random/Random/RandomEngine.h>
#include <CLHEP_EMBEDDED/Random/Random/RandGauss.h>
#include <CLHEP_EMBEDDED/Random/Random/RanluxEngine.h>
#include "hparammanager.h"
#include "hbeamfile.h"
#include "hevent.h"
#include "hconstants.h"
#include "hbooker.h"
#include "hvector.h"
#include "hlorentzvector.h"
#include "hWeightingInterface.h"

/*!
 * \brief Base-Class for physics-generators: all allowed generators will be derived from this
 *
 */

class HPhysicsGen
{

public:
    HPhysicsGen() {
        enabledBeamFile = false;
    };
    HPhysicsGen(CLHEP::HepRandom* _randGen, HParamManager* _paramMan, HEvent* _event, HBooker* _booker, HBooker* _ddBookMan);
    virtual ~HPhysicsGen();

    void setNextBeam(HBeamEntry _newBeam) {
        beamEntry = _newBeam;
    };
    bool useBeamFile() {
        enabledBeamFile = true;
        return true;
    };
    virtual void generateEvent();
    virtual void setParameters() {};
    virtual void addHistograms();
    
    inline double getWeightSum(){return weightCounter;};

protected:

    HBooker* bookMan;
    HBooker* ddBookMan;
    HEvent* event;

    HBeamEntry beamEntry;


    CLHEP::HepRandom* myRandom;
    CLHEP::RandGauss* myRandomGauss;
    CLHEP::RanluxEngine* GaussEngine;

    HParamManager* paramMan;

    /*! \brief generates a nu */
    bool generateNu();
    /*! \brief generates a qsq */
    bool generateQSQ();
    /*! \brief checks if we need to generate an elastic or inelastic event and rolls the dice for target nuclei type*/
    bool generateElastic();
    /*! \brief sets the beamparticle initially */
    bool setBeam();
    /*! \brief generates a phi_gammavirt and calcs the kinematics */
    bool generatePhiGamma();
    /*! \brief smears the missing mass e_miss */
    bool generateSmearing();
    /*! \brief does the generation of angles of the outgoing particle (meson or gamma) */
    bool generateOutgoingParticle();

    /*! \brief calculates the angle of the gamma* in the lepton plane*/
    bool calculatePhir();

    virtual bool generateMesonMass(); /*! this actually needs to be implemented by each specialized generator. */
    virtual void generateDiffractiveDissociation() {}; /*! this needs to be implemented by each specialized generator. */
    virtual void setUserVars();
    virtual void generateListOfParticles() {}; /*! generates a list of particles */

    virtual void printAuxVars();


    /*! \brief generates the mandelstam-t */
    bool generatet();



    bool enabledBeamFile;

    void ddDecayParticle(int _index, int _startIndex);
    void ddGenerateMultiplicities();
    void ddGenerateKinematics();
    void ddGenerateKinematics2Body();
    void ddGenerateParticles();
    void ddGeneratePhaseSpace();
    void ddGenerateLongitudinalPhaseSpace();
    void ddBalanceMultiParticles();
    HLorentzVector ddSumKinematics();
    void ddGotoLab();
    void ddAddParticleList();

    int ddGetNextFreeIndex(int _startIndex);



    bool coherent;
    bool incoherent;

    double amx2;
    double m_meson;

    bool doDD;
    bool hitNeutron;
    bool   isNeutronAfterHit;

    
    double weightCounter;



    double pft_, pftqsq_, pfnu_, pfoldqsq;


    double phir, iphir;
    double phi_out,theta_out;

    double theta_gamma_gammavirt;

    double theta_proton, beta_proton;

    //for dd
    double mx;
    int num_neutral, num_charged, num_particles;
    HParticle ddParticles[30];
    int ddIterations;

    //number of particles without DD adding in
    int baseNumberOfParticles;


    hepconst::nucType afterDecay;

};


#endif





