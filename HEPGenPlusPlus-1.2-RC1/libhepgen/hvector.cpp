/*
 * hvector.cpp
 *
 *  Created on: Jan 16, 2013
 *      Author: Christopher Regali
 */

#include "hvector.h"

HVector3::HVector3(double _x, double _y, double _z)
{
    x=_x;
    y=_y;
    z=_z;
}

void HVector3::roundToFloat()
{
  float xn,yn,zn;
  xn = static_cast<float>(x);
  yn = static_cast<float>(y);
  zn = static_cast<float>(z);
  x=xn;
  y=yn;
  z=zn;
}

void HVector3::rotateAxisAngle(const HVector3& _axis, const double _phi)
{
  HVector3 copyVec = *this;
  HVector3 axis = _axis;
  axis.normalize(1.0);
  HVector3 result = copyVec*cos(_phi) + (copyVec.crossProduct(axis))*sin(_phi) +axis*(1-cos(_phi))*(axis.dotProduct(copyVec));
  x = result.X();
  y = result.Y();
  z = result.Z();
}



HVector3::HVector3(const HVector3& _copy)
{
    x = _copy.X();
    y = _copy.Y();
    z = _copy.Z();
}


void HVector3::rotXToZ()
{
    double old[3] = {Y(),Z(),X()};
    setXYZ(old[0],old[1],old[2]);
}

double HVector3::at(unsigned int _index)
{
    switch(_index)
    {
    case 0:
        return X();
        break;
    case 1:
        return Y();
        break;
    case 2:
        return Z();
        break;
    default:
        return -1;
    }

}



void HVector3::rotPhi(double _angle)
{
    double xnew =  X()*cos(_angle)+Y()*sin(_angle);
    double ynew = -X()*sin(_angle)+Y()*cos(_angle);
    setX(xnew);
    setY(ynew);
}


double& HVector3::operator[](const int& _index)
{
    switch(_index)
    {
    case 0:
        return x;
        break;
    case 1:
        return y;
        break;
    case 2:
        return z;
        break;
    }
}

const double& HVector3::operator[](const int& _index) const
{
    switch(_index)
    {
    case 0:
        return x;
        break;
    case 1:
        return y;
        break;
    case 2:
        return z;
        break;
    }
}




void HVector3::print()
{

    std::cout << "HVector3: (" << X() << ", " << Y() << ", " << Z()<<")" << std::endl;

}


void HVector3::setXYZ(double _x, double _y, double _z)
{
    x=_x;
    y=_y;
    z=_z;
}

void HVector3::normalize(double _newlength)
{
    double len = length();
    x = _newlength / len * x;
    y = _newlength / len * y;
    z = _newlength / len * z;
}


void HVector3::fromStdVector(const std::vector< double >& _in)
{
    x = _in.at(0);
    y = _in.at(1);
    z = _in.at(2);
}



const std::vector<double> HVector3::toStdVector() const
{
    std::vector<double> tmp;
    tmp.push_back(x);
    tmp.push_back(y);
    tmp.push_back(z);
    return tmp;
}



double HVector3::dotProduct(const HVector3& _in) const
{
    return x*_in.X()+y*_in.Y()+z*_in.Z();
}

HVector3 HVector3::operator-(const HVector3& _in) const
{
    return HVector3(X() - _in.X(), Y() - _in.Y(), Z() - _in.Z());
}


HVector3 HVector3::operator+(const HVector3& _in) const
{
    return HVector3(X() + _in.X(), Y() + _in.Y(), Z() + _in.Z());
}

HVector3 HVector3::operator*(double _in)
{
  return HVector3(X() * _in, Y() * _in, Z() * _in);

}


HVector3 HVector3::crossProduct(const HVector3& _in) const
{

    double newx=y*_in.Z() - z*_in.Y();
    double newy=z*_in.X() - x*_in.Z();
    double newz=x*_in.Y() - y*_in.X();
    return HVector3(newx,newy,newz);
}

const double HVector3::length() const
{
    return sqrt(x*x+y*y+z*z);
}



void HVector3::extrapNeutral(const HVector3& _startPoint, HVector3& _endPoint, double distanceZ)
{
    //get the slopes
    double dx = X()/Z();
    double dy = Y()/Z();


    _endPoint.setX(_startPoint.X()+distanceZ*dx);
    _endPoint.setY(_startPoint.Y()+distanceZ*dy);
    _endPoint.setZ(_startPoint.Z()+distanceZ);
}





bool HVector3::operator== (const HVector3& _in)
{
    return (  (x==_in.X()) && (y == _in.Y()) && (z == _in.Z())  );
}

HVector3& HVector3::operator= (const HVector3& _in)
{
    if (this != &_in)
    {
        x= _in.X();
        y= _in.Y();
        z= _in.Z();
    }
    return *this;
}
