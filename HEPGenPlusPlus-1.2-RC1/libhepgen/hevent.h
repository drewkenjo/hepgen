/*!
 *  \file hphysicsgen.h
 *  \date Created on: Feb 7, 2013
 *  \author Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2013 All Right Reserved
 */


#define _USE_MATH_DEFINES


#ifndef HEVENTa_H_
#define HEVENTa_H_

#include <cmath>
#include "hparammanager.h"
#include "hlorentzvector.h"
#include "hparticle.h"
#include "hconstants.h"
#include "hrotmat.h"



/*!
 *
 *
 * \brief dataholder for the events
 */

typedef struct
{
    double nu;
    double PFnu;

    double qsq;
    double PFqsq;
    double emiss;

    double s;

    double wsq;
    double t;



    double tprim;

    double PFt;

    double z;
    double y;
    double xbj;
    double ptprim;

    vector<double> PARL;
    vector<double> USERVAR;


    double amx2,m_meson;
    double ddNum_Particles, ddNum_Charged, ddNum_Neutral;
    double ddMX;
    bool ddActive;
    double ddPMax;


    double gam2,epsilon,delta,flux;


    HParticle incBeamParticle;
    HParticle scatBeamParticle;
    HParticle targetParticle;

    HParticle gammaVirt;

    HParticle outPart1, outPart2, outPart3, outPart4, outPart5, outPart6;

    HParticle outPart1_Lab;
    HParticle outPart2_Lab;
    HParticle outPart3_Lab;
    HParticle outPart4_Lab;
    HParticle outPart5_Lab;
    HParticle outPart6_Lab;









    hepconst::eventType type;


    HLorentzVector CMS;
    HParticle recoil;

    vector<HParticle*> listOfParticles;

    int epp_polarized_longitudinal;

    double dummyW;


} HEventData;



/*!
 * \brief this class stores the data of the event and does the kinematic calculations.
 *
 */

class HEvent
{
public:
    /*!  \brief Constructor */
    HEvent(HParamManager* _paramMan);
    /*! \brief Destructor */
    ~HEvent() {}
    /*! \brief Resets the event to not happened */
    void reset();
    /*! \brief calculates the muon-vertex-kinematics */
    bool calculateMuonKinematics(double _flatrandphi);
    /*! \brief Sets a nu and recalcs the muon vectors */
    bool setNu(double _nu, double _phaseFactor);
    /*! \brief Sets a Q^2 and recalcs the muon vectors */
    bool setQSQ(double _qsq, double _phaseFactor);
    /*! \brief Sets a t */
    bool setT(double _t, double _phaseFactor, double _mesonMass, double _amx2);
    /*! \brief sets the missing energy */
    bool setEMiss(double _emiss) {
        dataStruct.emiss = _emiss;
        return true;
    };
    /*! \brief sets the parameters of the primary outgoing particle (meson, gamma) */
    bool setOutgoingParticle(double _phi, double _theta, double _pout);
    /*! \brief calcs the p_T and does therefore transform the primary particle to lab system */
    void	 calcPT();
    /*! \brief Transform the particle from the CMS system to the lab system */
    HLorentzVector goToLabSystem(HLorentzVector& _particle, HLorentzVector& _rot);

    /*! \brief rotates all the vectors from x-axis to z-axis for COMPASS-standard */
    void rotXToZ();

    /*! \brief returns the nu */
    double getNu() const {
        return dataStruct.nu;
    };
    /*! \brief returns the Y */
    double getY() const {
        return dataStruct.y;
    }
    /*! \brief returns the W^2 */
    double getWsq() const {
        return dataStruct.wsq;
    }
    /*! \brief returns the Xbj */
    double getXbj() const {
        return dataStruct.xbj;
    }
    /*! \brief returns the mandelstam t */
    double getT() const {
        return dataStruct.t;
    };
    /*! \brief returns the Q^2 */
    double getQsq() const {
        return dataStruct.qsq;
    };
    /*! \brief returns the energy s */
    double getS() const {
        return dataStruct.s;
    };  /*! \brief returns the missing energy */
    double getEmiss() const {
        return dataStruct.emiss;
    };

    /*! \brief cout Debug-information about the event */
    void printDebug();

    /*! \brief returns the virtual gamma */
    HParticle& getGammaVirt() {
        return dataStruct.gammaVirt;
    };
    HParticle& getBeam() {
        return dataStruct.incBeamParticle;
    };
    HLorentzVector& getCMS() {
        return dataStruct.CMS;
    };
    HParticle& getScat() {
        return dataStruct.scatBeamParticle;
    };
    HParticle& getOutPart1() {
        return dataStruct.outPart1;
    };
    HParticle& getOutPart1_Lab() {
        return dataStruct.outPart1_Lab;
    };


    HParticle& getOutPart2() {
        return dataStruct.outPart2;
    };
    HParticle& getOutPart2_Lab() {
        return dataStruct.outPart2_Lab;
    };

    HParticle& getOutPart3() {
        return dataStruct.outPart3;
    };
    HParticle& getOutPart3_Lab() {
        return dataStruct.outPart3_Lab;
    };


    HParticle& getRecoil() {
        return dataStruct.recoil;
    };


    HEventData* getStruct() {
        return &dataStruct;
    };

    double getTotalPhaseFactor() const {
        return 1./(dataStruct.PFnu*dataStruct.PFqsq*dataStruct.PFt);

    };/*! returns the total phase factor */
    void copyFrom(HEvent* _myEvent);





    void rotateEventToBeam(HVector3 _beam );



private:
    HParamManager* paramMan;
    HEventData dataStruct;

};


#endif




