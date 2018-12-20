/*!
 *  \file hdvcsgen.h
 *  \date Created on: Feb 9, 2013
 *  \author: Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2013 All Right Reserved
 */


#define _USE_MATH_DEFINES

#include "hphysicsgen.h"
#include "hconstants.h"
#include "hhelper.h"
#include "hevent.h"


#ifndef HDVCSGENa_H_
#define HDVCSGENa_H_



/*!
 * \brief Implementation of a DVCS-Generator
 *
 */
class HPhysicsGenDVCS : public HPhysicsGen
{
public:
    HPhysicsGenDVCS() {};
    HPhysicsGenDVCS(CLHEP::HepRandom* _randGen, HParamManager* _paramMan, HEvent* _event, HBooker* _booker, HBooker* _ddBooker);
    ~HPhysicsGenDVCS();
    void generateListOfParticles();



    void generateEvent();

    void calcWeights();
    bool generateMesonMass();

    void addHistograms();

    void setUserVars();
    void setParameters();






    /*! \brief this function does the calculation of the weight of dvcs - its the static generic version that can be called without the generator  -- interfaced struct version*/
    static double fundvcs(hWeightInterface& _data);
    /*! \brief this function does the calculation of the weight of interference terms - its the static generic version that can be called without the generator -- interfaced via struct version */
    static double funint(hWeightInterface& _data);
    /*! \brief this function does the calculation of the weight of Bethe-Heitler - static variant -- interfaced via struct version */
    static double funbh(hWeightInterface& _data);


    /*! \brief this function calculates the bethe heitler with the old propagators approximation by BMK */
    static double funbhOldProp(hWeightInterface& _data);



    /*! \brief this function does the calculation of the weight of dvcs - its the static generic version that can be called without the generator */
    static double fundvcs(double _xbj,double _qsq, double _t,double _y,double _s,double _nu,double _beamE, double _b0, double _xbj0, double _alphap);
    /*! \brief this function does the calculation of the weight of interference terms - its the static generic version that can be called without the generator */
    static double funint(double _t, double _y, double _qsq, double _xbj, HLorentzVector& _MuIn, HLorentzVector& _MuOut, HLorentzVector& _gammaOut, double _phir, double _S, double _nu, double _clept, double _slept, double _b0, double _xbj0, double _alphap);
    /*! \brief this function does the calculation of the weight of Bethe-Heitler - static variant */
    static double funbh(double _t, double _y, double _qsq, double _xbj, HLorentzVector& _MuIn, HLorentzVector& _MuOut, HLorentzVector& _gammaOut, double _S, double _nu, double _phir);



private:
    /*! \brief this function are actually kept for legacy reasons -- do not use! they will be removed in the future! use the static versions instead!*/
    double fundvcs(void);
    double funbh(void);
    double funint(void);


    

    double phir_trento;
    double weightdvcs, weightbh, weightint;
    double weightresult0,weightresult1,weightresult2,weight;
    double none;
    double HepgenInPhastRoundedBH,HepgenInPhastUnRoundedBH,HepgenInPhastRoundedBHAbs;

};


#endif


