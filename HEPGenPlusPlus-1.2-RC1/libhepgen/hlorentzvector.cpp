/*
 * hlorentzvector.cpp
 *
 *  Created on: Jan 16, 2013
 *      Author: Christopher Regali
 *
 *
 */

#include "hlorentzvector.h"


void HLorentzVector::print()
{
    double mass = getQuare();
    if (mass < 0 )
        mass = - sqrt(-mass);
    else
        mass = sqrt (mass);

    printf("VEC: %e\t%e\t%e\t%e\n",threeVector.X(),threeVector.Y(),threeVector.Z(),getEnergy());


}

void HLorentzVector::roundToFloat()
{
  threeVector.roundToFloat();
  float en = static_cast<float>(energy);
  energy = en;

}




HLorentzVector::HLorentzVector(double _x, double _y, double _z, double _e)
{
    threeVector.setXYZ(_x,_y,_z);
    energy=_e;
}




double HLorentzVector::getPtrans()
{
    double pTransSquare = threeVector.Z()*threeVector.Z() + threeVector.Y()*threeVector.Y();
    return sqrt(pTransSquare);
}




double HLorentzVector::getQuare() const
{
    return (getEnergy()*getEnergy())-threeVector.dotProduct(threeVector);
}


HVector3 const HLorentzVector::getVector() const
{

    return HVector3(threeVector.X(), threeVector.Y(), threeVector.Z());
}


void HLorentzVector::lorenf(double _u, const HLorentzVector& _ps)
{
    if (_u == _ps.getEnergy())
    {
        setVector(_ps.getVector());
        setEnergy(_ps.getEnergy());
    }

    double PF4 = (_ps.getEnergy()*getEnergy() - _ps.getVector().dotProduct(getVector()))    /  _u;
    double FN = (PF4 + getEnergy()) / (_ps.getEnergy() + _u );

    //set the new values
    setEnergy(PF4);
    setVector(HVector3 (getVector().X() - _ps.getVector().X() * FN ,
                        getVector().Y() - _ps.getVector().Y() * FN ,
                        getVector().Z() - _ps.getVector().Z() * FN ) );


}


void HLorentzVector::boost(double _u,const HLorentzVector& _ps)
{
    if (_u == _ps.getEnergy())
    {
        setVector(_ps.getVector());
        setEnergy(_ps.getEnergy());
    }

    double PF4 = (_ps.getEnergy()*getEnergy() + _ps.getVector().dotProduct(getVector()))/_u;
    double FN = (PF4 + getEnergy()) / (_ps.getEnergy() + _u );

    //set the new values
    setEnergy(PF4);
    setVector(HVector3 (getVector().X() + FN * _ps.getVector().X(),
                        getVector().Y() + FN * _ps.getVector().Y(),
                        getVector().Z() + FN * _ps.getVector().Z()) );

}



void HLorentzVector::setVector(const HVector3& _in)
{
    threeVector = _in;
}

HLorentzVector::HLorentzVector(HVector3 _vec, double _e)
{
    threeVector = _vec;
    energy = _e;
}


HLorentzVector::HLorentzVector(const HLorentzVector& _copy)
{
    threeVector = _copy.getVector();
    energy = _copy.getEnergy();
}


HLorentzVector HLorentzVector::operator-(const HLorentzVector& _in) const
{
    return HLorentzVector(getVector() - _in.getVector(), getEnergy() - _in.getEnergy());
}

HLorentzVector HLorentzVector::operator+(const HLorentzVector& _in) const
{
    return HLorentzVector(getVector() + _in.getVector(), getEnergy() + _in.getEnergy());
}

bool HLorentzVector::operator== (const HLorentzVector& _in)
{
    return (  (threeVector == _in.getVector()) && (energy == _in.getEnergy()) );
}

HLorentzVector& HLorentzVector::operator= (const HLorentzVector& _in)
{
    if (this != &_in)
    {
        threeVector = _in.getVector();
        energy = _in.getEnergy();
    }
    return *this;
}


void HLorentzVector::setLVectorAngular(double _P, double _theta, double _phi, double _E)
{
    threeVector.setXYZ(_P*cos(_theta),_P*sin(_theta)*cos(_phi),_P*sin(_theta)*sin(_phi));
    energy = _E;
}

