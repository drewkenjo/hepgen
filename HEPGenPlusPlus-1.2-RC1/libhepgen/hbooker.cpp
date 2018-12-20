/*
 * hbooker.h
 *
 *  Created on: Jan 29, 2013
 *      Author: Christopher Regali
 *
 *
 */

#include "hbooker.h"


void HBooker::dumpToFile()
{
    for (unsigned int i =0; i < backEnd.size(); i++)
        backEnd.at(i)->dumpToFile(fileName);
}
void HBooker::addASCII()
{
    HBookBackEndASCII* tmp = new HBookBackEndASCII(nameSuffix);
    setBackEnd(tmp);
}

HBooker::~HBooker()
{
    for (unsigned int i =0; i < backEnd.size(); i++)
        delete backEnd.at(i);
}


void HBooker::addROOT()
{
#ifdef USE_ROOT
    HBookBackEndROOT* tmp = new HBookBackEndROOT(nameSuffix);
    setBackEnd(tmp);
#else
    cout << "ROOT-HISTOGRAMMING DISABLED" << endl;
#endif
}



void HBooker::dumpToFile(std::string _fileName)
{
    fileName = _fileName;
    for (unsigned int i =0; i < backEnd.size(); i++)
        backEnd.at(i)->dumpToFile(_fileName);

}

void HBooker::fill()
{
    for (unsigned int i =0; i < backEnd.size(); i++)
        backEnd.at(i)->fill();

}

void HBooker::reset()
{
    backEnd.clear();
}

void HBooker::setBackEnd(HBookBackEnd* _newBackEnd)
{
    backEnd.push_back(_newBackEnd);

}

void HBooker::addHBook1D(double* _var, double* _weight, int _binCount, double _lowerBound, double _upperBound, std::string _title)
{
    //cout << "adding book " << _title << " to backends " << backEnd.size() << endl;
    for (unsigned int i =0; i < backEnd.size(); i++)
        backEnd.at(i)->addHBook1D(_var, _weight, _binCount, _lowerBound, _upperBound, _title);
}


void HBooker::addHBook2D(double* _varX, double* _varY, double* _weight, int _binCountX, int _binCountY, double _lowerBoundX, double _lowerBoundY, double _upperBoundX, double _upperBoundY, string _title)
{
    for (unsigned int i =0; i < backEnd.size(); i++)
        backEnd.at(i)->addHBook2D( _varX,  _varY,  _weight,  _binCountX,  _binCountY,  _lowerBoundX,  _lowerBoundY,  _upperBoundX,  _upperBoundY,  _title);
}

void HBooker::addHBook1F(float* _var, float* _weight, int _binCount, float _lowerBound, float _upperBound, std::string _title)
{
    //cout << "adding book " << _title << " to backends " << backEnd.size() << endl;
    for (unsigned int i =0; i < backEnd.size(); i++)
        backEnd.at(i)->addHBook1F(_var, _weight, _binCount, _lowerBound, _upperBound, _title);
}


void HBooker::addHBook2F(float* _varX, float* _varY, float* _weight, int _binCountX, int _binCountY, float _lowerBoundX, float _lowerBoundY, float _upperBoundX, float _upperBoundY, string _title)
{
    for (unsigned int i =0; i < backEnd.size(); i++)
        backEnd.at(i)->addHBook2F( _varX,  _varY,  _weight,  _binCountX,  _binCountY,  _lowerBoundX,  _lowerBoundY,  _upperBoundX,  _upperBoundY,  _title);
}





void HBooker::addHBook1I(int* _var, double* _weight, int _binCount, int _lowerBound, int _upperBound, std::string _title)
{

    for (unsigned int i =0; i < backEnd.size(); i++)
        backEnd.at(i)->addHBook1I(_var, _weight, _binCount, _lowerBound, _upperBound, _title);

}


void HBooker::fill(int _num)
{
    for (unsigned int i =0; i < backEnd.size(); i++)
        backEnd.at(i)->fill(_num);

}

