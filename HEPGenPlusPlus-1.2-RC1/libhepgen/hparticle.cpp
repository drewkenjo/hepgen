#include "hparticle.h"

double HParticle::getMass()
{
    double mass;
    if (!getSharpMass())
    {
        mass = getVector().getQuare();
        if (mass < 0)
            mass = -sqrt(-mass);
        else
            mass = sqrt(mass);
    }
    else
        mass = hephelp::getMassByID(getParticleType());

    return mass;
}


void HParticle::printDebugHeader()
{
    std::cout << std:: endl << " Type \t (K) \t Or \t d1 \t d2 \t Vec x \t y \t z \t p \t m " << std::endl;
}



void HParticle::printDebug()
{
    printf("PART:%+i\t%i\t%i\t%i\t%i\t(% 6.3e\t% 6.3e\t% 6.3e\t% 6.3e)\n",getParticleType(),getParticleAuxFlag(),getParticleOrigin(),getParticleDaughter1(),getParticleDaughter2(),getVector().getVector().X(),getVector().getVector().Y(),getVector().getVector().Z(),getEnergy());
}


void HParticle::scaleMomentum(double _factor)
{
    HVector3 vec = HVector3(getTVector().X()*_factor,
                            getTVector().Y()*_factor,
                            getTVector().Z()*_factor
                           );

    vector.setVector(vec);
    double mass = getMass();
    vector.setEnergy(sqrt(pow(mass,2) + pow(getTVector().length(),2) ));

}
