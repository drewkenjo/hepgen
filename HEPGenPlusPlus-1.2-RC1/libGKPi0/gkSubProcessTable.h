/*!
 *  \file gkSubProcessTable.h
 *  \date Created on: 8.4.2016
 *  \author: Christopher Regali <christopher.regali@cern.ch>
 */


#ifndef GKSUBTABLE_HH
#define GKSUBTABLE_HH

#include "libGKPi0.h"

#include <libgen.h>
#include <map>

using namespace std;


struct subProcessAmplitudeLine{
  double xlow,xhigh,weight,xpseudo;
  TComplex xlowT2,xlowT3,xhighT2,xhighT3;
  inline subProcessAmplitudeLine operator+(const subProcessAmplitudeLine& other) const {
         subProcessAmplitudeLine res;
    res.xlow =xlow+other.xlow;
    res.xhigh =xhigh+other.xhigh;
    res.weight =weight;
    res.xlowT2 =xlowT2+other.xlowT2;
    res.xlowT3 =xlowT3+other.xlowT3;
    res.xhighT2 =xhighT2+other.xhighT2;
    res.xhighT3 =xhighT3+other.xhighT3;
    return res;
  }
  inline subProcessAmplitudeLine operator-(const subProcessAmplitudeLine& other) const {
    subProcessAmplitudeLine res;
    res.xlow =xlow-other.xlow;
    res.xhigh =xhigh-other.xhigh;
    res.weight =weight;
    res.xlowT2 =xlowT2-other.xlowT2;
    res.xlowT3 =xlowT3-other.xlowT3;
    res.xhighT2 =xhighT2-other.xhighT2;
    res.xhighT3 =xhighT3-other.xhighT3;
        return res;
  }
  inline subProcessAmplitudeLine operator*(const double scale)const {
    subProcessAmplitudeLine res;
    res.xlow =xlow*scale;
    res.xhigh =xhigh*scale;
    res.weight =weight;
    res.xlowT2 =xlowT2*scale;
    res.xlowT3 =xlowT3*scale;
    res.xhighT2 =xhighT2*scale;
    res.xhighT3 =xhighT3*scale;
    return res;
  }
    
};

inline double StrToDouble(string _arg)
{
    stringstream sstr;
    double i;
    sstr << _arg;
    sstr >> i;
    return i;
}


inline vector<string> explodeStringWhiteSpace(string _input)
{
    vector<string> array;

    int inputlength = _input.length();
    int whitespaceLength = 1;
    int i = 0;
    int k = 0;
    while (i < inputlength) {
        int j = 0;
        while (i + j < inputlength && j < whitespaceLength
                && (_input[i + j] == ' ' || _input[i + j] == '\t'))
            j++;
        if (j == whitespaceLength) { //found whitespace
            if (_input.substr(k, i - k) != "")
                array.push_back(_input.substr(k, i - k));
            i += whitespaceLength;
            k = i;
        } else {
            i++;
        }
    }
    if (_input.substr(k, i - k) != "")
        array.push_back(_input.substr(k, i - k));
    return array;
};


class gkSubProcessTableCache{
  public:
    /*! \brief initializes the cache from the tableIndex file in _fileName
     *  \return the code of error, 0 is fine, -1 $HEPGEN not found, -2 tableIndex not found,
     *          -3 could not open one of the cache-files */
    int loadCache(string _fileName);
    /*! \brief loads the default installation cache 
     */
    int loadCache();
    
    /*! \brief interpolates subprocessAmplitudes in the qsq and xbj range. Then calculates the full amplitudes using libGKPi0
     */
     GKPI0::amplitude getAmpsForKine(double _qsq, double _xbj, double _t);
    
     
    /*! \brief interpolates subprocessAmplitudes in the qsq and xbj range, using forward and backwards extrapolation. Then calculates the full amplitudes using libGKPi0
     */
    GKPI0::amplitude  getAmpsForKineDoubleWay(double _qsq, double _xbj, double _t);

     
    subProcessAmplitudeLine getLineForKine(double _qsq, double _xbj,int i);
  private:
    subProcessAmplitudeLine loadLineFromFile(string _fileName);
    vector<double> qbins;
    vector<double> wbins;
    vector<double> xpseudo;
    vector<double> weights;
    //map [qsq_bin][w_bin]
    map<int, map<int, vector< subProcessAmplitudeLine > > > cacheMap;
    
};

#endif
