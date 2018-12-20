/*!
 *  \file hvector.h
 *  \date Created on: Jan 16, 2013
 *  \author Christopher Regali <christopher.regali@cern.ch
 *  Copyright (c) 2013 All Right Reserved
 */

#ifndef HVECTOR_H_
#define HVECTOR_H_
#include <cmath>
#include <vector>
#include <iostream>

/*!
 * \brief a simple 3-value-double-precision vector.
 *
 * This is somewhat similar to the TVector3 of root but it has only some basic features needed for the hepgen-project.
 * Not taking ROOT's TVector3 is a choice of because of the design-guideline: NO DEPENDENCIES
 *
 */
class HVector3
{
public:
    /*! \brief raw constructor from 3 doubles */
    HVector3(double _x=0.0, double _y=0.0, double _z=0.0);
    //todo: maybe add a constructor for reading std::vector or other stuff.
    /*! \brief copy constructor */
    HVector3(const HVector3& _copy);
    /*! \brief constructor from std::vector */
    HVector3(const std::vector< double >& _in) {
        fromStdVector(_in);
    };
    
    /*! \brief breaks down double precision to float precision - for test use only! */
    void roundToFloat();

    /*! \brief get the component at position _index (has to be 0,1,2 or else it returns nan */
    double at(unsigned int _index);


    /*! \brief prints the vector to cout */
    void print();

    /*! \brief destructor */
    ~HVector3() {};



    //static getters
    /*! \brief Read-Only access to X */
    double X() const {
        return x;
    };

    /*! \brief Read-Only access to Y */
    double Y() const {
        return y;
    };

    /*! \brief Read-Only access to Z */
    double Z() const {
        return z;
    };
    //setters


    /*! \brief Sets a new X */
    void setX(double _x) {
        x=_x;
    };

    /*! \brief Sets a new Y */
    void setY(double _y) {
        y=_y;
    };

    /*! \brief Sets a new Z */
    void setZ(double _z) {
        z=_z;
    };

    /*! \brief Sets a new X,Y,Z */
    void setXYZ(double _x, double _y, double _z);  /*! gives write access to components */




    /*! \brief Linearly extrapolates along this vector from start vector to end vector */
    void extrapNeutral(const HVector3& _startPoint, HVector3& _endPoint, double distanceZ);


    /*! \brief rotates vector around axis with angle phi in radians */
    void rotateAxisAngle(const HVector3& _axis, const double _phi);


    //just add some more math stuff here as needed
    /*! \brief returns a simple dot-product between the 2 vectors*/
    double dotProduct(const HVector3& _in) const;
    /*! \brief returns the vectorial-product (AxB) between the 2 vectors*/
    HVector3 crossProduct(const HVector3& _in) const;
    /*! \brief returns the length */
    const double length() const;

    /*! \brief normalizes the vector to a new length */
    void normalize(double _newlength);

    /*! \brief rotates the X to the Z axis (COMPASS compat) */
    void rotXToZ();


    /*! \brief rotates the x/y plane with angle _angle */
    void rotPhi(double _angle);


    /*! \brief returns a std::vector with the entries */
    const std::vector<double> toStdVector() const;

    /*! \brief initializes from an std-vector */
    void fromStdVector(const std::vector<double>& _in);


    //overloaded operators

    bool operator== (const HVector3& _in);
    HVector3& operator= (const HVector3& _in);
    HVector3 operator- (const HVector3& _in)const;
    HVector3 operator+ (const HVector3& _in)const;
    HVector3 operator* (double _in);
    const double& operator[] (int const& _index)const;
    double& operator[] (int const& _index);


private:
    double x,y,z;


};



#endif /* HVECTOR_H_ */


