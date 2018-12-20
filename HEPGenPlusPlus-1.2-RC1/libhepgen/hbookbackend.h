/*!
 *  \file hbookbackend.cc
 *  \date Created on: Jan 29, 2013
 *  \author: Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2013 All Right Reserved
 */

#include <iostream>
#include <string>



#ifndef HBOOKBACKEND_H_
#define HBOOKBACKEND_H_

/*!
 * \class HBookBackEnd
 * \brief this class is just for inheritance (purely virtual except constructor). Its a basic BackEnd class for histograms that is used in hbookbackendASCII and hbookbackendROOT
 *
 */
class HBookBackEnd
{
public:
    HBookBackEnd(std::string _nameSuffix) {
        nameSuffix = _nameSuffix;
    }
    virtual ~HBookBackEnd() {};
    virtual int addHBook1D(double* _var, double* _weight, int _binCount, double _lowerBound, double _upperBound, std::string _title) {
        return 0;
    };
    virtual int addHBook2D(double* _varX,double* _varY , double* _weight, int _binCountX,int _binCountY , double _lowerBoundX,double _lowerBoundY, double _upperBoundX,double _upperBoundY, std::string _title) {
        return 0;
    };


    virtual int addHBook1I(int* _var, double* _weight, int _binCount, int _lowerBound, int _upperBound, std::string _title) {
        return 0;
    };
    
    virtual int addHBook1F(float* _var, float* _weight, int _binCount, float _lowerBound, float _upperBound, std::string _title) {
        return 0;
    };
    virtual int addHBook2F(float* _varX,float* _varY , float* _weight, int _binCountX,int _binCountY , float _lowerBoundX,float _lowerBoundY, float _upperBoundX,float _upperBoundY, std::string _title) {
        return 0;
    };
    


    virtual std::string getTitle(int _num) {
        return "";
    };

    virtual void dumpToFile(std::string _fileName) {};
    virtual void fill() {};
    virtual int fill(int _histNum) {
        return 0;
    };

protected:
    std::string nameSuffix;
};





#endif
