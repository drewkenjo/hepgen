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


#ifndef HPAMGENa_H_
#define HPAMGENa_H_



/*!
 * \brief Implementation of a PAM-Guichon-BH-Generator
 *
 */
class HPhysicsGenPAMBH : public HPhysicsGen
{
public:
    HPhysicsGenPAMBH() {};
    HPhysicsGenPAMBH(CLHEP::HepRandom* _randGen, HParamManager* _paramMan, HEvent* _event, HBooker* _booker, HBooker* _ddBooker);
    ~HPhysicsGenPAMBH();
    void generateListOfParticles();



    void generateEvent();

    void calcWeights();
    bool generateMesonMass();

    void addHistograms();

    void setUserVars();
    void setParameters();






    /*! \brief this function does the calculation of the weight of Bethe-Heitler - static variant -- interfaced via struct version */
    static double funbh(hWeightInterface& _data);




    /*! \brief this function does the calculation of the weight of Bethe-Heitler - static variant */
    static double funbh(double _t, double _qsq, double _xbj,  double _elept, double _phir);



private:
  double funbh();
    double phir_trento;
    double weightdvcs, weightbh, weightint;
    double weightresult0,weightresult1,weightresult2,weight;
    double none;
    double HepgenInPhastRoundedBH,HepgenInPhastUnRoundedBH,HepgenInPhastRoundedBHAbs;

};


#endif


