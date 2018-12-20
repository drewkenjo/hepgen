/*!
 *  \file hpionicdata.h
 *  \date Created on: Jan 29, 2013
 *  \author Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2013 All Right Reserved
 */


#include <iostream>
#include <fstream>
#include <string>
#include "hhelper.h"


#ifndef HPIONICDATA_H_
#define HPIONICDATA_H_




using namespace std;









/*!
 * \brief loads the pionic cross section data from goloskokov and kroll via the pi0_input.dat
 *
 */

class HPionicData
{
public:
    HPionicData();
    HPionicData(string _fileName);

    ~HPionicData() {};
    void resetData();
    /*! \brief loads the data from the file _fileName */
    void loadFile(string _fileName);
    
    /*! \brief loads the data from the ng-format-file _fileName
     * This file-format features embedded information about the binning
     */
    void loadFileNG(string _fileName);


    /*! \brief returns the vector to the pionic data */
    vector< vector<double> > getData() const {
        return pi0data;
    };


    vector<double>* getBinningW() {
        return &pi0w;
    }
    vector<double>* getBinningQ2() {
        return &pi0qsq;
    }
    vector<double>* getBinningTPR() {
        return &pi0tpr;
    }





    double getPi0w(int _i) const {
        return pi0w[_i];
    };
    double getPi0Qsq(int _j) const {
        return pi0qsq[_j];
    };
    double getPi0tpr(int _k) const {
        return pi0tpr[_k];
    };
    double getPi0sigl(int _i, int _j, int _k) const {
        return pi0sigl.at(_i).at(_j).at(_k);
    }
    double getPi0sigt(int _i, int _j, int _k) const {
        return pi0sigt.at(_i).at(_j).at(_k);
    }

    /*! \brief these function gets the nearest bin in w*/
    int getBinW(double _w);
    /*! \brief these function gets the nearest bin in Q^2*/
    int getBinQsq(double _qsq);
    /*! \brief these function gets the nearest bin in tprime*/
    int getBinTPrim(double _tprim);




private:
    /*! \brief fills the legacy arrays with the same indexing as in the original hepgen */
    void fillLegacyArrays(void);

    string fileName;
    vector< vector<double> > pi0data;


    /*! \brief one dimensional arrays that index the others */
    vector<double> pi0w;
    /*! \brief one dimensional arrays that index the others */
    vector<double> pi0qsq;
    /*! \brief one dimensional arrays that index the others */
    vector<double> pi0tpr;



//     static const double pi0w[11];
//     static const double pi0qsq[9];
//     static const double pi0tpr[16];
//

    /*! \brief these are the three dimensinal arrays in form of crappy std::vectors but since we dont want to use libboost we have to do it this way.  */
    vector< vector< vector<double> > > pi0sigl, pi0sigt;
};



#endif


