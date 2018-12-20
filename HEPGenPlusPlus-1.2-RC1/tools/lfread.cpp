#include "lfread.h"



void lfr::lPrintLogo(void)
{
    cout << "  _    ___ ___ _____ ___      ___             _         "<<endl;
    cout << " | |  | __| _ \\_   _/ _ \\ ___| _ \\___ __ _ __| |___ _ _ "<<endl;
    cout << " | |__| _||  _/ | || (_) |___|   / -_) _` / _` / -_) '_|"<<endl;
    cout << " |____|___|_|   |_| \\___/    |_|_\\___\\__,_\\__,_\\___|_|  "<<endl;






    cout << "  Easy-Access library to LEPTOv2 binary format in C and C++  "<<endl;
    cout << "  C. Regali  "<<endl;

}




int lfr::lLoadFile(string _fileName, fstream& fileHandle)
{
    fileHandle.open(_fileName.c_str(), ios::in | ios::binary);

    //check if file is good at all
    if (!fileHandle.good())
        return -1;

    //check for the headers now
    bool largeHeaders;
    int retVal = lIsGFortranFile(fileHandle,largeHeaders);
    //could not be identified!
    if (retVal < 0)
        return -2;
    else if (largeHeaders)
        return 0;
    else if (!largeHeaders)
        return 1;
}


int lfr::lFreeBuffer(lfr::lEvent& event)
{
    delete[] event.beamParts;
}


int lfr::lReadHeader(fstream& fileHandle, lfr::lHeader& readResult, bool _largeHeaders, bool _beQuiet)
{
    //first we rewing completely!
    fileHandle.seekg(0);
    //read the header-start-bit
    int byteInt;
    float byteFloat;
    char byteChar;
    int size = sizeof(int);
    if (_largeHeaders)
        size = sizeof(int64_t);
    fileHandle.read((char*) &byteInt,size);
    //safety check
    if (byteInt != 228)
    {
      if (!_beQuiet)
        printf("startbit failed! read: %i \n",byteInt);
        return -1;
    }

    fileHandle.read((char*) &byteInt, sizeof(int)); // lepto binary file version = 2

    for (unsigned int i = 0; i < 14; i++) {
        fileHandle.read((char*) &byteFloat, sizeof(float));
        readResult.cutl[i] = byteFloat;
    }
    for (unsigned int i = 0; i < 20; i++) {
        fileHandle.read((char*) &byteInt, sizeof(int));
        readResult.lst[i] = byteInt;
    }

    fileHandle.read((char*)&byteFloat, sizeof(int)); // unidentified
    readResult.unidentified[0] = byteFloat;
    fileHandle.read((char*)&byteFloat, sizeof(int)); // unidentified
    readResult.unidentified[1] = byteFloat;

    for (unsigned int i = 0; i < 19; i++) {
        fileHandle.read((char*) &byteFloat, sizeof(float));
        readResult.parl[i] = byteFloat;
    }

    fileHandle.read((char*) &byteFloat, sizeof(float)); // unidentified = 0.1
    readResult.unidentified[2] = byteFloat;
    fileHandle.read((char*) &byteInt, size); // end flag of header = 228

    if (byteInt != 228) {
      if (!_beQuiet)
        printf("stopbit failed! read: %i \n",byteInt);
        return -2;
    }
    return 0;
}



int lfr::lNextEvent(fstream& fileHandle, lfr::lEvent& event, bool _largeHeaders)
{
    int byteInt;
    float byteFloat;
    int size = sizeof(int);
    if (_largeHeaders)
        size = sizeof(int64_t);



    fileHandle.read((char*) &byteInt, size); // start flag (432 for DVCS, 512 for rho0)
    int start = byteInt;
    //the startflag is the bytesize of the event
    event.size = start;
    if (fileHandle.eof())
        return -2;

    // num particles
    fileHandle.read((char*) &byteInt, sizeof(int));
    event.nBeamParticle = byteInt;

    for (unsigned int i = 0; i < 11; i++) {
        fileHandle.read((char*) &byteInt, sizeof(int));
        event.lst[i + 20] = byteInt;
    }

    for (unsigned int i = 0; i < 11; i++) {
        fileHandle.read((char*) &byteFloat, sizeof(float));
        event.parl[i + 19] = byteFloat;
    }

    fileHandle.read((char*) &byteFloat, sizeof(float));
    event.x_bj = byteFloat;
    fileHandle.read((char*) &byteFloat, sizeof(float));
    event.y = byteFloat;
    fileHandle.read((char*) &byteFloat, sizeof(float));
    event.w2 = byteFloat;
    fileHandle.read((char*) &byteFloat, sizeof(float));
    event.q2 = byteFloat;
    fileHandle.read((char*) &byteFloat, sizeof(float));
    event.nu = byteFloat;

    for (unsigned int i = 0; i < 20; i++) {
        fileHandle.read((char*) &byteFloat, sizeof(float));
        event.uservar[i] = byteFloat;
    }


    event.beamParts = new lBeamParticle[event.nBeamParticle];

    for (unsigned int i = 0; i < 5; i++) {
        for (unsigned int j = 0; j < event.nBeamParticle; j++) {
            fileHandle.read((char*) &byteInt, sizeof(int));
            event.beamParts[j].k[i] = byteInt;
        }
    }

    for (unsigned int i = 0; i < 5; i++) {
        for (unsigned int j = 0; j < event.nBeamParticle; j++) {
            fileHandle.read((char*) &byteFloat, sizeof(float));
            event.beamParts[j].p[i] = byteFloat;
        }
    }


    fileHandle.read((char*) &byteInt, size); // stop flag
    if (byteInt != start) {
        return -1;
    }
    return 0;
}


int lfr::lIsGFortranFile(fstream& fileHandle,bool &largeHeaders)
{
    int byteInt;
    float byteFloat;
    lHeader tmp; 
    //try short headers
    if (lfr::lReadHeader(fileHandle,tmp,false) == 0){
      largeHeaders = false;
      fileHandle.seekg(0);
      return 0;
    }
    //rewind and try long headers
    fileHandle.seekg(0);
    if (lfr::lReadHeader(fileHandle,tmp,true) == 0){
      largeHeaders = true;
      fileHandle.seekg(0);
      return 0;
    }
    //give up
    fileHandle.seekg(0);
    return -1;

}









//






