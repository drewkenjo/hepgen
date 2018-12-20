#include "lfwrite.h"



int  lfr::lOpenFile(string _fileName, fstream& fileHandle, bool _useLongHeaders)
{
    fileHandle.open(_fileName.c_str(), std::ios::out | std::ios::binary);
    if (fileHandle.good())
        return LFWR_SUCCESS;
    return LFWR_ERROR;
}

int  lfr::lWriteHeader(fstream& fileHandle, lfr::lHeader& writeMe, bool _largeHeaders)
{
    //here we write the header of the lepto file - its just a binary dump of the params
    //i will describe the file-format here in the comments

    //write the start-flag
    lfr::putInt(fileHandle,228);

    //first we write a 2 for the IDMPU- Variable in hepgen
    lfr::putInt(fileHandle,2);

    //now we put the lepto parameters which are 14 doubles, but we convert them to float for length-reasons
    for (int i = 0; i < 14; i++)
        lfr::putFloat(fileHandle,writeMe.cutl[i]);

    //then we put the lepto switches which are 20 ints
    for (int i = 0; i < 20; i++)
        lfr::putInt(fileHandle,writeMe.lst[i]);
    //then we put 2 more unidentified
    fileHandle.write((char*) &writeMe.unidentified[0],sizeof(int));
    fileHandle.write((char*) &writeMe.unidentified[1],sizeof(int));
    //then we put the PARL-vector which contains all the parameters
    //i actually dont know how many entries there have to be i hope it just fits
    for (int i=0; i < 19; i ++)
        lfr::putFloat(fileHandle,writeMe.parl[i]);
    //again putting a float for no idea why!
    fileHandle.write((char*) &writeMe.unidentified[2],sizeof(float));
    //write the stop-flag
    lfr::putInt(fileHandle,228);
    fileHandle.flush();
    return 0;
}



int  lfr::lWriteNextEvent(fstream& fileHandle, lfr::lEvent* event, bool _largeHeaders)
{
    int byteInt;
    float byteFloat;
    int size = sizeof(int);
    if (_largeHeaders)
        size = sizeof(int64_t);


    //eventsize
    lfr::putInt(fileHandle,event->size);
    //num of particles:
    lfr::putInt(fileHandle,event->nBeamParticle);
    //next write lst and parl and kinematics

    for (unsigned int i = 0; i < 11; i++)
        lfr::putInt(fileHandle,event->lst[i + 20]);

    for (unsigned int i = 0; i < 11; i++)
        lfr::putFloat(fileHandle,event->parl[i + 19]);


    lfr::putFloat(fileHandle,event->x_bj);
    lfr::putFloat(fileHandle,event->y);
    lfr::putFloat(fileHandle,event->w2);
    lfr::putFloat(fileHandle,event->q2);
    lfr::putFloat(fileHandle,event->nu);

    for (unsigned int i = 0; i < 20; i++)
        lfr::putFloat(fileHandle,event->uservar[i]);


    for (unsigned int i = 0; i < 5; i++)
        for (unsigned int j = 0; j < event->nBeamParticle; j++)
            lfr::putInt(fileHandle,event->beamParts[j].k[i]);


    for (unsigned int i = 0; i < 5; i++)
        for (unsigned int j = 0; j < event->nBeamParticle; j++)
            lfr::putFloat(fileHandle,event->beamParts[j].p[i]);


    lfr::putInt(fileHandle,event->size);
    return 0;

}

void lfr::putFloat(fstream& fileHandle, float _f)
{
    fileHandle.write( (char*) &_f, sizeof(float));
}

void lfr::putInt(fstream& fileHandle, int32_t _i)
{
    fileHandle.write( (char*) &_i, sizeof(int32_t));
}

void lfr::putInt64(fstream& fileHandle, int64_t _i)
{
    fileHandle.write( (char*) &_i, sizeof(int64_t));
}




