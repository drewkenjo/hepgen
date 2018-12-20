/*!
 *  \file hrotmat.h
 *  \date Created on: July 30, 2014
 *  \author Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2014 All Right Reserved
 */

#include "hvector.h"
#include <cstdio>
#include <cstring>


#ifndef HROTMAT_H
#define HROTMAT_H


/*! \brief a simple rotation matrix class */
class HRotMat
{
public:
    /*! \brief constructor setting matrix to all identity */
    HRotMat() {
        resetMatrix();
    };



    /*! \brief constructor setting matrix according to angles and axis */
    HRotMat(double _theta, double _phi, const HVector3& _axis) {
        setFromThetaPhiVector(_theta,_phi,_axis);
    };


    /*! \brief destructor */
    ~HRotMat() {};

    HVector3 rotateVector(HVector3 _in);



    void resetMatrix();

    /*! \brief sets the rotation matrix from two angles and an axis */
    void setFromThetaPhiVector(double _theta, double _phi, const HVector3& _axis);

    /*! \brief prints the matrix */
    void printDebug();

private:
    double components[3][3];
};










#endif
