#include "houtputbackendLEPTOGF.h"


void HOutPutBackEndLEPTOGF::closeFile()
{
    outFile.flush();
    outFile.close();
}

HOutPutBackEndLEPTOGF::HOutPutBackEndLEPTOGF()
{
    //set up our buffers for writing once
    floatSize = sizeof(float);
    floatBuffer = new char[floatSize];


    intSize = sizeof(int32_t);
    intBuffer = new char[intSize];
}


void HOutPutBackEndLEPTOGF::dumpEvent(HEvent* _event)
{
    //write event start-flag
//     if (eventPointer->getStruct()->type == hepconst::DVCS)
//         putInt64(432);
//     else if (eventPointer->getStruct()->type == hepconst::RHO0 || eventPointer->getStruct()->type == hepconst::PI0)
//         putInt64(512);
    int startflag = 192 + eventPointer->getStruct()->listOfParticles.size()*40;
    putInt64(startflag);

    //first we write the number of outgoing particles

    putInt(eventPointer->getStruct()->listOfParticles.size());
    //cout << "put size: " << eventPointer->getStruct()->listOfParticles.size() << endl;


    //now we write the LST-variables again, but just 20 to 30
    for (int i=20; i < 31; i ++)
        putInt(paramMan->getStruct()->LST.at(i));



    //now we write 10 parls again
    for (int i=19; i < 30; i ++)
    {
        //cout << "PARL AT " << i << " with value " << eventPointer->getStruct()->PARL.at(i)<< endl;
        putFloat(eventPointer->getStruct()->PARL.at(i));
    }
    //then we put the kinematic variables
    putFloat(eventPointer->getXbj());
    putFloat(eventPointer->getY());
    putFloat(eventPointer->getWsq());
    putFloat(eventPointer->getQsq());
    putFloat(eventPointer->getNu());

    //now we put the uservars in
    for (int i = 0; i < 20; i++)
        putFloat(eventPointer->getStruct()->USERVAR.at(i));



    //then we can start with the dump of the vectors
    //first we write the KS-Values 21 for incoming, 1 for outgoing

    //we always have three incoming, so its fixed
    for (unsigned int i=0; i < eventPointer->getStruct()->listOfParticles.size(); i++)
        putInt(eventPointer->getStruct()->listOfParticles.at(i)->getParticleAuxFlag());


    //now we put the particle ids

    for (unsigned int i=0; i < eventPointer->getStruct()->listOfParticles.size(); i++)
        putInt(eventPointer->getStruct()->listOfParticles.at(i)->getParticleType());

    //then we put the origins
    for (unsigned int i=0; i < eventPointer->getStruct()->listOfParticles.size(); i++)
        putInt(eventPointer->getStruct()->listOfParticles.at(i)->getParticleOrigin());
    //then we put the Daughters1
    for (unsigned int i=0; i < eventPointer->getStruct()->listOfParticles.size(); i++)
        putInt(eventPointer->getStruct()->listOfParticles.at(i)->getParticleDaughter1());

    //then we put the Daughters2
    for (unsigned int i=0; i < eventPointer->getStruct()->listOfParticles.size(); i++)
        putInt(eventPointer->getStruct()->listOfParticles.at(i)->getParticleDaughter2());



    //now put px, py, pz,E,m
    for (unsigned int i=0; i < eventPointer->getStruct()->listOfParticles.size(); i++)
        putFloat(eventPointer->getStruct()->listOfParticles.at(i)->getTVector().X());
    for (unsigned int i=0; i < eventPointer->getStruct()->listOfParticles.size(); i++)
        putFloat(eventPointer->getStruct()->listOfParticles.at(i)->getTVector().Y());
    for (unsigned int i=0; i < eventPointer->getStruct()->listOfParticles.size(); i++)
        putFloat(eventPointer->getStruct()->listOfParticles.at(i)->getTVector().Z());
    for (unsigned int i=0; i < eventPointer->getStruct()->listOfParticles.size(); i++)
        putFloat(eventPointer->getStruct()->listOfParticles.at(i)->getEnergy());
    for (unsigned int i=0; i < eventPointer->getStruct()->listOfParticles.size(); i++)
    {
        double mass;
        //everything except for rho0, gamma and pi0 should be sharp so we can lookup the mass
        if (!eventPointer->getStruct()->listOfParticles.at(i)->getSharpMass())
        {
            mass = eventPointer->getStruct()->listOfParticles.at(i)->getVector().getQuare();
            if (mass < 0)
                mass = -sqrt(-mass);
            else
                mass = sqrt(mass);
        }
        else
        {
            mass = hephelp::getMassByID(eventPointer->getStruct()->listOfParticles.at(i)->getParticleType());
        }
        putFloat(mass);
    }







    //finally we write event stop-flag
//     if (eventPointer->getStruct()->type == hepconst::DVCS)
//         putInt64(432);
//     else if (eventPointer->getStruct()->type == hepconst::RHO0 || eventPointer->getStruct()->type == hepconst::PI0)
//         putInt64(512);
    putInt64(startflag);
    outFile.flush();
}

void HOutPutBackEndLEPTOGF::putFloat(float _f)
{
    outFile.write( (char*) &_f, floatSize);
}

void HOutPutBackEndLEPTOGF::putInt(int32_t _i)
{
    outFile.write( (char*) &_i, intSize);
}

void HOutPutBackEndLEPTOGF::putInt64(int64_t _i)
{
    outFile.write( (char*) &_i, sizeof(int64_t));
}



void HOutPutBackEndLEPTOGF::initFile(string _fileName, HEvent* _eventPointer)
{
    fileName = _fileName+"_binary.bin";
    eventPointer = _eventPointer;
    outFile.open(fileName.c_str(), std::ios::out | std::ios::binary);



}



void HOutPutBackEndLEPTOGF::setParams(HParamManager* _params)
{
    paramMan = _params;
}



void HOutPutBackEndLEPTOGF::writeHeader()
{
    //here we write the header of the lepto file - its just a binary dump of the params
    //i will describe the file-format here in the comments

    //write the start-flag
    //this is the bytesize of the header
    putInt64(228);


    //first we write a 2 for the IDMPU- Variable in hepgen
    putInt(2);

    //now we put the lepto parameters which are 14 doubles, but we convert them to float for length-reasons
    for (int i = 0; i < 14; i++)
        putFloat(paramMan->getStruct()->CUT.at(i));

    //then we put the lepto switches which are 20 ints
    for (int i = 0; i < 20; i++)
        putInt(paramMan->getStruct()->LST.at(i));



    //then we put 2 more ints there for whatever reason:
    putInt(2);
    putInt(1);



    //then we put the PARL-vector which contains all the parameters
    //i actually dont know how many entries there have to be i hope it just fits

    for (int i=0; i < 19; i ++)
        putFloat((float)paramMan->getStruct()->PARL.at(i));

    //again putting a float for no idea why!
    putFloat(0.1);

    //write the stop-flag
    putInt64(228);








    outFile.flush();

}


