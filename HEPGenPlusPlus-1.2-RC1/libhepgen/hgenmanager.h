/*!
 *  \file hgenmanager.h
 *  \date Created on: Feb 6, 2013
 *  \author Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2013 All Right Reserved
 */





#ifndef HGENMANAGER_H_
#define HGENMANAGER_H_

#include "config.h"
#include "hbooker.h"
#include "hhelper.h"
#include "hparammanager.h"
#include "hpionicdata.h"
#include "hevent.h"
#include "hphysicsgen.h"
#include "hdvcsgen.h"
#include "hrhogen.h"
#include "hpigen.h"
#include "hphigen.h"
#include "hrhoplusgen.h"
#include "homegagen.h"
#include "hpamgen.h"
#include "hjpsigen.h"
#include "homegagenpigamma.h"


#ifdef USE_EXPERIMENTAL
// #include "hVGGgen.h"
#include "hMosseGen.h"
#endif


#include "houtput.h"
#include "hbeamfile.h"

#include "CLHEP_EMBEDDED/Random/Random/RandomEngine.h"
#include "CLHEP_EMBEDDED/Random/Random/RanluxEngine.h"
#include "CLHEP_EMBEDDED/Random/Random/Random.h"



/*!
 * \brief This class manages the runs by calling a physics-plugin repeatedly and saving everything to a backend.
 *
 */
class HGenManager
{
public:

    /*! \brief gets the actual Instance of this singleton */

    static  HGenManager* getInstance() {
        if (instance == NULL) {
            instance = new HGenManager();
        }
        return instance;
    };

    ~HGenManager();
    /*! \brief starts the whole run */
    double startRun(const void (*_callBack)(int));
    
    /*! \brief forces the write of the last generated event to the file
     *  also fills the histograms
     */
    void saveLastEvent();
    
    /*! \brief saves the output file and writes the histograms to the file*/
    void closeOutPut();
    
    /*! \brief creates some standard histograms for the generator kinematics */
    void addStdHistos();
    
    /*! \brief creates generatorspecific histograms */
    void addGenHistos(){if (generator != NULL)generator->addHistograms();}

    /*! \brief enables beamfile-functionality */
    void enableBeamFile(string _fileName, HBeamType _meGusta);
    /*! \brief Enables ROOT-Backend for histograms */
    void enableROOTBook();
    /*! \brief Enables ASCII-Backend for histograms */
    void enableASCIIBook();
    /*! \brief Enabled Output in General */
    void enableOutPut(string _filename);
    /*! \brief Enabled Output in LEPTO-Format either with or without gfortran-compat-headers */
    void enableLeptoOutput(bool _gfortran);

    /*! \brief returns the event-pointer for a generated event */
    HEvent* getEvent() {
        return bookingEvent;
    };
    /*! \brief sets the seed */
    void setSeed(int _seed);

    /*! \brief generates a single event - no beamfile involved! you have to rotate it yourself! */
    void oneShot();


    /*! \brief enabled debug-output */
    void enableDebug() {
        debug = true;
    };

    /*! \brief disables debug-output */
    void disableDebug() {
        debug = false;
    };

    /*! \brief prints debug if debug was enabled */
    void printDebug();

    /*! \brief loads the parameters from a datacard */
    void loadSettings (string _fileName);



    /*! \brief setup the generator for oneshot-library-use */
    void setupGeneratorOneShot();

    /*! \brief setup the generator */
    void setupGenerator();

    /*! \brief setup an external param manager in case you want to do some fiddling by hand */
    void setParamManager(HParamManager* _params);


    /*! \brief gets the paramManager actually in use, so you can do some fiddling with it */
    HParamManager* getParamManager(void) {
        return paramMan;
    };



private:

    static HGenManager* instance;

    /*! \brief  */
    HGenManager();

    bool enabledBeamFile;
    bool enabledROOTBook;
    bool enabledASCIIBook;
    int seed;

    HBeamType wantedBeamType;


    CLHEP::HepRandom* randMan;
    CLHEP::HepRandomEngine* randEng;

    HBeamEntry thisBeam;
    double aux[10];




    bool debug;



    HBeamFile* beamFile;
    HPhysicsGen* generator;
    HParamManager* paramMan;
    HBooker* hbookMan;
    HBooker* ddBookMan;
    HEvent* bookingEvent;
    HOutPut* outputMan;
    HEvent* myEvent;

};







#endif

