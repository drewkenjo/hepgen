/*!
 *  \brief hbooker.h
 *  \date Created on: Jan 29, 2013
 *  \author: Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2013 All Right Reserved
 */


#ifndef HBOOKER_H_
#define HBOOKER_H_


#include <iostream>
#include <string>
#include <vector>
#include "config.h"
#include "hbookbackend.h"
#include "hbookbackendASCII.h"
#ifdef USE_ROOT
#include "hbookbackendROOT.h"
#endif

/*!
 * \brief Manages the HBOOKs, routes all calls to the listed backends
 *
 *
 */
class HBooker
{
public:
    HBooker(std::string _nameSuffix) {
        nameSuffix = _nameSuffix;
    }
    ~HBooker();
    /*! \brief resets the backendlist */
    void reset();
    /*! \brief dumps all data via all the backends */
    void dumpToFile(std::string _fileName);
    /*! \brief does use the internal fileName variable */
    void dumpToFile();
    /*! \brief adds a double histogram to all backends */
    void addHBook1D(double* _var, double* _weight, int _binCount, double _lowerBound, double _upperBound, std::string _title);
    /*! \brief adds a double histogram to all backends */
    void addHBook2D(double* _varX, double* _varY , double* _weight, int _binCountX,int _binCountY , double _lowerBoundX,double _lowerBoundY, double _upperBoundX,double _upperBoundY, std::string _title);
    
    /*! \brief adds a float histogram to all backends */
    void addHBook1F(float* _var, float* _weight, int _binCount, float _lowerBound, float _upperBound, std::string _title);
    /*! \brief adds a float histogram to all backends */
    void addHBook2F(float* _varX, float* _varY , float* _weight, int _binCountX,int _binCountY , float _lowerBoundX,float _lowerBoundY, float _upperBoundX,float _upperBoundY, std::string _title);
    
    
    /*! \brief adds an integer histogram to all backends */
    void addHBook1I(int* _var, double* _weight, int _binCount, int _lowerBound, int _upperBound, std::string _title);
    /*! \brief adds an ASCII-Backend to the backendlist */
    void addASCII();
    /*! \brief adds a ROOT-Backend to the backendlist */
    void addROOT();
    /*! \brief fills all the histograms in all the backends */
    void fill();
    /*! \brief fills only the histogram with the number _num */
    void fill(int _num);

private:
    /*! \brief adds a new backendpointer to the backendlist */
    void setBackEnd(HBookBackEnd* _newBackEnd);
    std::vector<HBookBackEnd*> backEnd;
    std::string fileName;
    std::string nameSuffix;
};







#endif
