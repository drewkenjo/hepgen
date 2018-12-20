/*!
 *  \file houtput.h
 *
 *  \date created on: Feb 13, 2013
 *  \author Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2013 All Right Reserved
 */


#ifndef HOUTPUT_H_
#define HOUTPUT_H_

#include <vector>
#include "houtputbackend.h"
#include "config.h"

#ifdef USE_ROOT
#include "houtputbackendROOT.h"
#endif

#include "houtputbackendLEPTO.h"
#include "houtputbackendLEPTOGF.h"

#include "hparammanager.h"
/*! \brief This class manages the outputs
 *  It keeps track of all the running output backends and messages them.
 *
 */
class HOutPut
{
public:
    /*! \brief constructor */
    HOutPut(HParamManager* _params);
    ~HOutPut();
    /*! \brief enables ascii backend */
    void enableASCII();
    /*! \brief enables LEPTO backend */
    void enableLEPTO(bool _gfortran);
    /*! \brief enables ROOT backend (NOT WORKING!!!) */
    void enableROOT();
    /*! \brief sets the filename */
    void setFileName(std::string _fileName);
    /*! \brief writes the file-headers in all the outputs */
    void writeHeader();
    /*! \brief opens the files and makes them ready for dumping data */
    void initFile(HEvent* _event);
    /*! \brief dumps an event in all active backends */
    void dumpEvent(HEvent* _event);
    /*! \brief closes file after use */
    void closeFile();


private:


    std::vector<HOutPutBackEnd*> backEndList;
    std::string fileName;
    HParamManager* paramMan;

    bool enabledROOT, enabledASCII, enabledLEPTO;



};








#endif
