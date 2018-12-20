/*!
 *  \file houtputbackendROOT.h
 *  \date Created on: Feb 13, 2013
 *  \author Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2013 All Right Reserved
 */



#ifndef HOUTPUTBACKENDROOT_H_
#define HOUTPUTBACKENDROOT_H_


#include "houtputbackend.h"
#include "hevent.h"
#include <TFile.h>
#include <TTree.h>

/*! \brief Output-Backend for events in ROOT-Format currently NOT WORKING */
class HOutPutBackEndROOT: public HOutPutBackEnd
{
public:
    HOutPutBackEndROOT() {};
    ~HOutPutBackEndROOT() {
        delete eventTree;
    };

    void writeHeader();
    void dumpEvent();
    void initFile(std::string _fileName, HEvent* _eventPointer);
    void closeFile();
    void setParams(HParamManager* _params);


private:
    std::string fileName;
    TFile* eventFile;
    TTree* eventTree;
    HEvent* eventPointer;
    int eventBranch;

};



#endif


