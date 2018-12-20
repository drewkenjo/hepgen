/*!
 *  \file hparammanager.h
 *  \date Created on: Jan 16, 2013
 *  \author Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2013 All Right Reserved
 */

#ifndef HPARAMM_H_
#define HPARAMM_H_
#include <iostream>
#include <string>
#include <vector>
#include "hcardparser.h"
#include "hhelper.h"
#include "hpionicdata.h"

using namespace std;


/*!
 * \brief This Struct contains all the Parameters needed by Lepto (setpar.F)
 *
 * It is managed by the class HParamManager
 */

typedef struct
{
    //run-parameters
    long NEVENT;
    double ELEPT;
    double EHADR;
    int IPART;
    int I_BEAMREAD;
    int IPROC;

    //lepto-parameters for target
    vector<double> PARL;

    //common block for vector meson decay
    vector<int> IVECM;


    //lepto-switches
    vector<int> LST;

    //lepto-soft-cuts
    vector<double> CUT;



    //aux-parameters for generators
    vector<string> aux1;
    vector<string> aux2;
    vector<string> aux3;
    vector<string> aux4;
    vector<string> aux5;



    //other parameters
    double THMAX;		// Scattered muon acceptance - theta_max
    double    ivecm;		// vector meson flag
    double    kcdec;		// absolute value of the decay partiles ID
    double    alf;		// value of the parameter for A dependence
    double    atomas;		// atomic mass (average) for the target
    double    probc;		// fraction of coherent events
    double     bcoh;		// slope of the nuclear formfactor
    double    bin;		// slope for the production on a nucleon
    double    rdiss;		// ratio of incoherent events with nucFleon dissociation
    // to incoherent elastic events
    double    clept;		// lepton beam charge - default electron
    double    slept;		// lepton beam polarisation - default no polarisation
    double    B0;		// parameter for x,t correlation |
    double    xbj0;		// parameter for x,t correlation |=> for DVCS only
    double    alphap;		// parameter for x,t correlation |

    //schildknecht_parameters
    double aksi2;
    double amt2_rho;
    double amt2_phi;
    double aml2_rho;
    double aml2_phi;
    double slopeR;
    std::string outFile;


    //target type in enum
    hepconst::nucType targetType;

} HParams;



/*!
 * \brief This class manages all the parameters needed for the main execution. (setpar.F)
 *
 * As a minimal set of dependencies it is purely relying on libstdc++ which comes with every modern OS an even slc5!
 * Most parameters here are hardcoded in the original - i will continue with that practice for now, even if this is no good style for C++
 *
 */
class HParamManager
{
public:
    /*! \brief  standard-constructor */
    HParamManager(HPionicData* _data);

    /*! \brief constructor for fiddling with the parameters */
    HParamManager(HParams* _paramStruct);

    /*! \brief standard-constructor from filename */
    HParamManager(string _filename);

    /*! \brief  standard-destructor */
    ~HParamManager() {};
    /*! \brief resets the Parameters according to  setpar.F*/
    void resetData();
    /*! \brief reads the parameters from the .data-file via HCardParser and sets the parameters in the HParams accordingly. Keywords are hardcoded here. */
    void readFile ( string _filename );
    /*! \brief returns the number of events to generate */
    int getNumEvents() const {
        return paramStruct.NEVENT;
    };
    /*! \brief returns the physics programme to use */
    int getPhysics() const {
        return paramStruct.IPROC;
    };
    /*! \brief returns the chosen meson programme */
    vector<int> getIVECM() const {
        return paramStruct.IVECM;
    };
    /*! \brief returns the read-in minimal value for nu */
    double getNuMin() const {
        return paramStruct.CUT.at(8);
    };
    /*! \brief returns the read-in maximal value for nu */
    double getNuMax() const {
        return paramStruct.CUT.at(9);
    };
    /*! \brief returns the read-in minimal value for Q^2 */
    double getQsqMin() const {
        return paramStruct.CUT.at(4);
    };
    /*! \brief returns the read-in maximal value for Q^2 */
    double getQsqMax() const {
        return paramStruct.CUT.at(5);
    };
    /*! \brief returns the read-in minimal value for t mandelstam */
    double gettMin() const {
        return paramStruct.PARL.at(15);
    };
    /*! \brief returns the read-in maximal value for t mandelstam */
    double gettMax() const {
        return paramStruct.PARL.at(16);
    };
    /*! \brief returns the read-in minimal value for y */
    double getyMin() const {
        return paramStruct.CUT.at(2);
    };
    /*! \brief returns the read-in maximal value for y */
    double getyMax() const {
        return paramStruct.CUT.at(3);
    };
    
    /*! \brief returns the read-in minimal value for y */
    double getXbjMin() const {
        return paramStruct.CUT.at(0);
    };
    /*! \brief returns the read-in maximal value for y */
    double getXbjMax() const {
        return paramStruct.CUT.at(1);
    };


    /*! \brief returns the lepton charge */
    double getclept() const {
        return paramStruct.clept;
    };
    /*! \brief returns the lepton polarisation */
    double getslept() const {
        return paramStruct.slept;
    };

    /*! \brief  returns the slope b for incoherent */
    double getSlopeIncoherent() const {
        return paramStruct.bin;
    };
    /*! \brief returns the slope b for coherent */
    double getSlopeCoherent() const {
        return paramStruct.bcoh;
    };

    /*! \brief returns the particle typenum for beam particle */
    int getBeamPart() const {
        return paramStruct.IPART;
    };
    /*! \brief returns the xbj0-parameters for x,t correlation */
    double getXbj0() const {
        return paramStruct.xbj0;
    };
    /*! \brief returns the B0-parameters for x,t correlation */
    double getB0() const {
        return paramStruct.B0;
    };
    /*! \brief returns the alphap-parameters for x,t correlation */
    double getalphap() const {
        return paramStruct.alphap;
    };

    std::string getOutName() const {
        return paramStruct.outFile;
    }


    /*! \brief returns the whole parameter struct */
    HParams* getStruct() {
        return &paramStruct;
    };





    /*! \brief returns the read-in count for events to be generated */
    long   getNEvent() const {
        return paramStruct.NEVENT;
    };
    /*! \brief returns the lepton energy */
    double getBeamE() const {
        return paramStruct.ELEPT;
    };

    /*! \brief returns a pointer to the pionic-data class*/
    HPionicData* getPionicData() {
        return pionicData;
    };

    /*! \brief Sets the pointer to the pionic data class */
    void setPionicData(HPionicData* _data) {
        pionicData = _data;
    }


    static void resetParamStruct(HParams* toReset);


    /*! \brief  returns the parameter-list from a key */
    vector<string> getKeyContents(string _keyName) {
        return paramParser.getKeyContents(_keyName);
    };




private:
    HPionicData* pionicData;
    string redFile;
    HParams paramStruct;
    HCardParser paramParser;

};






#endif

