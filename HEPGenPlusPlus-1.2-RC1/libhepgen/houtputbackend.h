/*!
 *  \file houtputbackend.h
 *
 *  \date Created on: Feb 13, 2013
 *  \author Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2013 All Right Reserved
 */


#ifndef HOUTPUTBACKEND_H_
#define HOUTPUTBACKEND_H_

#include <string>
#include "hevent.h"
/*! \brief The virtual base-class for all the file-output-backends
 * all file-outputs are derived from this class
 *
 */
class HOutPutBackEnd
{
public:
    /*! \brief constructor */
    HOutPutBackEnd() {};

    virtual ~HOutPutBackEnd() {};
    /*! \brief opens the file  or creates if needed*/
    virtual void initFile(std::string _fileName, HEvent* _eventPointer) {};
    /*! \brief writes possible headers */
    virtual void writeHeader() {};
    /*! \brief dumps a data-event */
    virtual void dumpEvent(HEvent* _event) {};
    /*! \brief close the file */
    virtual void closeFile() {};
    /*! \brief set additional parameters */
    virtual void setParams(HParamManager* _params);

protected:
    HParamManager* paramMan;
    HEvent* eventPointer;

};




#endif

