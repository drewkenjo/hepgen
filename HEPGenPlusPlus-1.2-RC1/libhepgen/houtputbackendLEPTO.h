/*!
 *  \file houtputbackendROOT.h
 *
 *  \date Created on: Feb 13, 2013
 *  \author Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2013 All Right Reserved
 */


#ifndef HOUTPUTBACKENDLEPTO_H_
#define HOUTPUTBACKENDLEPTO_H_


#include "houtputbackend.h"
#include "hparammanager.h"
#include "hevent.h"
#include <iostream>
#include <fstream>
#include <inttypes.h>


/*! \brief Output-Backend for events in binary LEPTO-format with extended gfortran headers */
class HOutPutBackEndLEPTO: public HOutPutBackEnd
{
public:
    HOutPutBackEndLEPTO();
    ~HOutPutBackEndLEPTO() {
        delete[] floatBuffer, delete[] intBuffer;
    };

    void writeHeader();
    void dumpEvent(HEvent* _event);
    void initFile(std::string _fileName, HEvent* _eventPointer);
    void closeFile();
    void setParams(HParamManager* _params);


private:
    void putFloat(float _f);
    void putInt(int32_t _i);
    void putInt64(int64_t _i);
    char* floatBuffer;
    char* intBuffer;
    int floatSize;
    int intSize;

    std::string fileName;
    ofstream outFile;
    int eventBranch;

};



#endif

