/*!
 *  \file hparticle.h
 *  \date Created on: Feb 11, 2013
 *  \author Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2013 All Right Reserved
 */

#ifndef HPARTICLE_H_
#define HPARTICLE_H_


#include "hlorentzvector.h"
#include "hhelper.h"

#include <cstdio>
/*!
 * \brief a class representing a particle with a pdg-group-particle-code and a lorentz-vector
 *
 */

class HParticle
{
public:
    /*! \brief Constructor*/
    HParticle() {};

    /*! \brief Constructor from HLorentzVector and Particle Type*/
    HParticle(HLorentzVector _vec, int _particleType) {
        vector = _vec;
        partType = _particleType;
    };

    /*! \brief Prints Debug-Table-Header to cout*/
    void printDebugHeader();

    /*! \brief prints particle to cout*/
    void printDebug();



    /*! \brief Destructor*/
    ~HParticle() {};

    /*! \brief sets particle type - types are defined in hconstants.h */
    void setParticleType(int _type) {
        partType = _type;
    };
    /*! \brief gets the particle type */
    int getParticleType(void) const {
        return partType;
    };
    /*! \brief sets the particle origin */
    void setParticleOrigin(int _origin) {
        partOrigin = _origin;
    };
    /*! \brief gets the particle origin */
    int getParticleOrigin(void) const {
        return partOrigin;
    };
    /*! \brief sets the particle daughter particle start-line */
    void setParticleDaughter1(int _daughter) {
        partDaughter1 = _daughter;
    };

    /*! \brief gets the particle daughter particle start-line */
    int getParticleDaughter1(void) const {
        return partDaughter1;
    };

    /*! \brief sets the particle daughter particle stop-line */
    void setParticleDaughter2(int _daughter) {
        partDaughter2 = _daughter;
    };

    /*! \brief gets the particle daughter particle stop-line */
    int getParticleDaughter2(void) const {
        return partDaughter2;
    };


    /*! \brief Sets if the particle has to have a sharp mass from PDG or its mass is calculated from the HLorentzVector */
    void setSharpMass(bool _flag) {
        useSharpMass = _flag;
    };


    /*! \brief Gets if the particle has to have a sharp mass from PDG or its mass is calculated from the HLorentzVector */
    bool getSharpMass()const {
        return useSharpMass;
    };


    /*! \brief Gets the Mass (according to wether it has to have a sharp mass or not) */
    double getMass();

    /*! \brief Returns the Mass-Square of the particle - uses getMass() function for setting sharp masses right */
    double getMassSq() {
        return pow(getMass(),2.);
    };

    /*! \brief Scales the Momentum by a factor */
    void scaleMomentum(double _factor);

    /*! \brief gets the Momentum */
    inline double getMomentum() {
        return getTVector().length();
    };


    /*! \brief is set 1 for particles that are alive, and 21 for particles that are decayed */
    void setParticleAuxFlag(int _flag) {
        partAuxFlag = _flag;
    };
    /*! \brief gets particle aux flag */
    int getParticleAuxFlag(void) const {
        return partAuxFlag;
    };



    /*! \brief sets the HLorentzVector */
    void setVector(const HLorentzVector& _newVec) {
        vector = _newVec;
    };

    /*! \brief gets the HLorentzVector */
    HLorentzVector& getVector(void) {
        return vector;
    };

    /*! \brief gets directly the HVector3 */
    void setTVector(HVector3 _newVec) {
        vector.setVector(_newVec);
    };


    /*! \brief gets directly the HVector3 */
    HVector3 getTVector(void) const {
        return vector.getVector();
    };
    /*! \brief gets the Energy */
    double getEnergy(void) const {
        return vector.getEnergy();
    };
    /*! \brief rotates the particle momentum from x to z axis for compass */
    void RotXToZ() {
        vector.rotXToZ();
    };





private:
    HLorentzVector vector;
    bool useSharpMass;
    int		 partType;
    int 		 partOrigin;
    int     	 partAuxFlag;
    int 		 partDaughter1; //index of first daughter particle
    int 		 partDaughter2; //index of last daughter particle





};

#endif


