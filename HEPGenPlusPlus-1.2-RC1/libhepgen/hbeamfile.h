/*!
 *  \file hbeamfile.h
 *  \date Created on: Sept 23, 2013
 *  \author: Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2013 All Right Reserved
 */


#ifndef HBEAMFILE
#define HBEAMFILE

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <inttypes.h>
#include "hlorentzvector.h"
#include "hvector.h"
#include "hconstants.h"


using namespace std;


enum HBeamType
{
    Beam=1,
    Halo=2,
    Both=3
};

/*! \brief Dataholder struct for a single Beam-Entry with momentum and energy */
struct  HBeamEntry
{
    HLorentzVector momentum;
    HVector3 position;
    HBeamType type;
    double energy;
    float slopex,slopey;
};


/*!
 * \class HBeamFile
 * \brief This Class implements the beamfile-reading support for HEPGen++
 * This file loads the beamfile - automatically checks for f77 or gfortran headersize and makes the beamentries available to the generator
 */
class HBeamFile
{
public:
    /*! \brief Constructor */
    HBeamFile();
    /*! \brief Constructor with direct file-initialization */
    HBeamFile(string _fileName);
    /*! \brief Destructor */
    ~HBeamFile();
    /*! \brief Loads a beamfile */
    void loadFile(string _fileName);

    /*! \brief writes the beamfile to a given filename */
    void saveFile(string _fileName,bool _header64 = false);

    /*! \brief Gets the next beam-entry from the file */
    HBeamEntry& getNextEntry(HBeamType _typeWanted);
    /*! \brief Returns the number of entries in the beamfile */
    int getNEntries() {
        return beamList.size();
    };
    /*! \brief check if file is ok to use */
    bool good() {
        return isInitialized;
    }

    /*! \brief add new entry to list */
    inline void addEntry (HBeamEntry& newEntry){
      beamList.push_back(newEntry);
    }
    
    /*! \brief clears the list of beamentries */
    inline void clearList(){beamList.clear();};



private:

    HBeamEntry& rollOver(HBeamType _typeWanted);
    bool is4BitFile;
    bool is4BitHeader();
    bool isInitialized;
    ifstream inFile;
    ofstream oFile;
    string fileName;
    vector<HBeamEntry>::iterator index;
    vector<HBeamEntry> beamList;
    float tmpDouble;
    int tmpInt;
    int64_t tmpInt64;

    inline int getLong() {
        inFile.read((char*) &tmpInt64, sizeof(int64_t));
        return tmpInt64;
    };
    inline int getShort() {
        inFile.read((char*) &tmpInt, sizeof(int));
        return tmpInt;
    };

    inline int getInt() {
        if (is4BitFile) return getShort();
        else return getLong();
    };
    inline float getDouble() {
        inFile.read((char*) &tmpDouble, sizeof(float));
        return tmpDouble;
    };
    inline void putFloat(float _f)
    {
        oFile.write( (char*) &_f, sizeof(float));
    };
    inline void putInt(int32_t _i)
    {
        oFile.write( (char*) &_i, sizeof(int));
    };
    inline void putInt64(int64_t _i)
    {
        oFile.write( (char*) &_i, sizeof(int64_t));
    };



};


















#endif


