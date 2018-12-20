#include "hbeamfile.h"


HBeamFile::HBeamFile()
{
    isInitialized = false;
    index = beamList.begin();

}


HBeamFile::HBeamFile(std::string _fileName)
{
    cout << "HBeamFile reading file: " << _fileName << endl;
    isInitialized = false;
    loadFile(_fileName);
    index = beamList.begin();
}

/*this function checks for f77 header-size (4bit) if it fails, it should be 8 bit */
bool HBeamFile::is4BitHeader()
{

    //first we try to find a 4bit int with 24 content.
    while (!inFile.eof()) {
        if (getShort() != 24)
            continue;
        //expect we have one here
        //next we expect a type - it has to be 1 or 2, else our headerlength-assumption was wrong!
        int type = getShort();
        if (type == 1 || type ==2)
            return true;
    }

    //try 8 bit

    inFile.seekg(0);
    while (!inFile.eof()) {
        if (getLong() != 24)
            continue;
        //expect we have one here
        //next we expect a type - it has to be 1 or 2, else our headerlength-assumption was wrong!
        int type = getShort();
        if (type == 1 || type ==2)
            return false;
    }

    //so if we are still here its not a valid beamfile
    isInitialized = false;
    cout << "HBeamFile: Error! Neither 4 Bit-Headers nor 8 Bit-Headers found! Beamfile is corrupted!" << endl;

    return false;


}



HBeamFile::~HBeamFile()
{


}

void HBeamFile::saveFile(string _fileName, bool _header64)
{
    oFile.open(_fileName.c_str(),ios::out | ios::binary);
    if (!inFile.good()) {
        printf("Error, could not open outfile %s for writing! Aborting save!\n",_fileName.c_str());
        return;
    }

    for (unsigned int myIndex = 0; myIndex < beamList.size(); myIndex++) {
        //write header
        if (_header64)
            putInt64(24);
        else
            putInt(24);
        //type
        if (beamList.at(myIndex).type == Beam)
            putInt(1);
        else
            putInt(2);
        //position
        putFloat(beamList.at(myIndex).position.X());
        putFloat(beamList.at(myIndex).position.Y());
        //slopes of momentum
        putFloat(beamList.at(myIndex).slopex);
        putFloat(beamList.at(myIndex).slopey);
        //energy
        putFloat(beamList.at(myIndex).energy);
        //stopflag
        if (_header64)
            putInt64(24);
        else
            putInt(24);
    }
    oFile.close();
    printf("Successfully wrote beamfile %s!\n",_fileName.c_str());
}



void HBeamFile::loadFile(string _fileName)
{
    fileName = _fileName;
    inFile.open(_fileName.c_str(), ios::in | ios::binary);
    if (inFile.good())
        isInitialized = true;
    else
        isInitialized = false;



    //now we check if this is a f77 or gfortran binary fileName
    //difference is the different size in start and stop-flag
    //to check this, we read the header in 4 bit and see if the second (first typeflag) is valid (i.e. int 1 or  int 2)
    //else if it does not validate - retry with 8 bit int length



    is4BitFile =  is4BitHeader();
    //reset the filepointer
    inFile.seekg(0);

    if (is4BitFile && isInitialized)
        cout << "HBeamFile: 4 Bit-Headers detected!" << endl;
    else
        cout << "HBeamFile: 8 Bit-Headers detected!" << endl;




    if (isInitialized)
    {
        while (!inFile.eof()) {

            //wait for startflag
            if (getInt() != 24)
                continue;
            HBeamEntry newEntry;
            int type = getShort();
            switch (type)
            {
            case 1:
                newEntry.type = Beam;
                break;
            case 2:
                newEntry.type = Halo;
                break;
            default:
                continue;
                break;
            }
            // cout << type << endl;

            newEntry.position.setX(getDouble());

            newEntry.position.setY(getDouble());

            double slopeX = getDouble();
            newEntry.slopex = slopeX;
            double slopeY = getDouble();
            newEntry.slopey = slopeY;
            double energy = getDouble();
            int flagEnd = getInt();
            //endflag here or event corrupted
            if (flagEnd != 24)
                continue;
            //now do the kinematics of the beam and make the beam ready
            double sinx,siny,px;
            sinx = sin(slopeX*1e-3);

            siny = sin(slopeY*1e-3);


            //here called energy is just absolute momentum
            px = energy/(sqrt(1.+pow(sinx,2)+ pow(siny,2)));

            HVector3 momThree = HVector3(px, px*sinx, px*siny);

            newEntry.momentum.setVector(momThree);

            newEntry.momentum.setEnergy(sqrt(momThree.dotProduct(momThree)+hepconst::w2mu));
            newEntry.energy = sqrt(momThree.dotProduct(momThree)+hepconst::w2mu);

            beamList.push_back(newEntry);
        }
    }

    cout << "HBeamFile: read " << beamList.size() << " entries." << endl;

}



HBeamEntry& HBeamFile::getNextEntry(HBeamType _typeWanted)
{
    if (index != beamList.end()) {
        while (((*index).type != _typeWanted) && (_typeWanted != Both))
            if (index + 1 != beamList.end())
                index++;
            else
                return rollOver(_typeWanted);
        return (*index++);
    } else
        return rollOver(_typeWanted);
}


HBeamEntry& HBeamFile::rollOver(HBeamType _typeWanted)
{
    cout
            << "T4BeamFile: Warning! Beamfile is too small! Starting over with the same file again."
            << std::endl;
    index = beamList.begin();
    return getNextEntry(_typeWanted);
}


