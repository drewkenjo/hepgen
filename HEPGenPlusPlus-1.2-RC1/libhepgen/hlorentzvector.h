/*!
 *  \file hlorentzvector.h
 *  \date Created on: Jan 16, 2013
 *  \author Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2013 All Right Reserved
 */

#ifndef HLVECTOR_H_
#define HLVECTOR_H_


#include "hvector.h"
#include <iostream>
#include <cstdio>


/*!
 * \brief a simple Lorentz-4-Vector
 *
 * This is somewhat similar to the TLorentzVector of root but it has only some basic features needed for the hepgen-project.
 * Not taking ROOT's TLorentzVector is a choice of because of the design-guideline: NO DEPENDENCIES
 *
 */
class HLorentzVector
{
public:
    /*! \brief raw constructor from 4 doubles, note that e could be e or t, depending on the interpretation of the other 3 doubles */
    HLorentzVector(double _x=0.0, double _y=0.0, double _z=0.0, double _e=0.0);
    /*! \brief  copy constructor */
    HLorentzVector(const HLorentzVector& _copy);
    /*! \brief  constructor from 3-vector of momentum and energy */
    HLorentzVector(HVector3 _vec, double _e);
    ~HLorentzVector() {};

    /*! \brief returns 3-vector of momentum CONST */
    const HVector3 getVector() const;

    HVector3& getVectorRef() {
        return threeVector;
    };

    /*! \brief returns energy */
    double getEnergy() const {
        return energy;
    };
    
    /*! \brief returns the invariant mass */
    double getMass() const{
        return sqrt(getQuare());
    }
    
    /*! \brief rounds this vector to float precision only */
    void roundToFloat();

    /*! \brief gets the square of the vector */
    double getQuare() const;

    /*! \brief rotates the vector from x to z axis */
    void rotXToZ() {
        threeVector.rotXToZ();
    };
    
    
    /*! \brief makes a dot product of two lorentz vectors */
    double dotProduct(const HLorentzVector& _in) const{
      return (-getEnergy()*_in.getEnergy() + getVector().dotProduct(_in.getVector()));
    }


    /*! \brief prints vector to cout */
    void print();


    //setter
    /*! \brief  sets 3-vector */
    void setVector(const HVector3& _in);
    /*! \brief sets the energy */
    void setEnergy(double _e) {
        energy = _e;
    };

    /*! \brief  sets the vector-properties from spherical coordinates of the momentum and the energy */
    void setLVectorAngular(double _P, double _theta, double _phi,double _E);


    /*! \brief gets transverse momentum wrs z-axis*/
    double getPtrans();

    /*! \brief Boosts the vector */
    void boost(double _u, const HLorentzVector& _ps);
    /*! \brief lorenf reimpls from fortran */
    void lorenf(double _u, const HLorentzVector& _ps);


    //overloaded operators

    bool operator== (const HLorentzVector& _in);
    HLorentzVector& operator= (const HLorentzVector& _in);

    HLorentzVector operator- (const HLorentzVector & _in) const;
    HLorentzVector operator+ (const HLorentzVector & _in) const;
    






private:
    HVector3 threeVector;
    double energy;

};


#endif /* HLVECTOR_H_ */

