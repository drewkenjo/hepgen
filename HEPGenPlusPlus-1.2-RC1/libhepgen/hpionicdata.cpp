#include "hpionicdata.h"




HPionicData::HPionicData()
{
    //resetData();
}

HPionicData::HPionicData(string _fileName)
{
    resetData();
    cout << "loading Pi0-File: " << _fileName << endl;
    loadFile(_fileName);
}

void HPionicData::loadFileNG(string _fileName)
{
    fileName = _fileName;
    ifstream inputFile;
    inputFile.open(fileName.c_str(),ios::in);
    if (!inputFile.is_open() || !inputFile.good()) {
        printf("Cannot open requested table %s \n",_fileName.c_str());
        exit(-1);
    }
    string line;
    bool inData = false;
    while (getline(inputFile,line)) {
        if (line.length() < 5)
            continue;
        if (inData) {
//             printf("%s\n",line.c_str());
            //end data if requested
            if (line == "#END_DATA") {
                inData=false;
                //finalize
                fillLegacyArrays();
                continue;
            }
            else {
                vector<string>boomResult = hephelp::explodeStringWhiteSpace(line);
                for (int i = 0; i < 5; i++) {
//                     printf("%i %f\n",i,hephelp::StrToDouble(boomResult.at(i)));
                    pi0data.at(i).push_back(hephelp::StrToDouble(boomResult.at(i)));
                }
            }
        }
        //special line
        else if (line[0] == '#') {
            vector<string> boomResult = hephelp::explodeStringWhiteSpace(line);
            //resize the arrays
            if (boomResult.front() == "#NUM_BINS_W")
                pi0w.resize(hephelp::StrToInt(boomResult.back()));
            else if (boomResult.front() == "#NUM_BINS_TPRIME")
                pi0tpr.resize(hephelp::StrToInt(boomResult.back()));
            else if (boomResult.front() == "#NUM_BINS_QSQ")
                pi0qsq.resize(hephelp::StrToInt(boomResult.back()));
            //now fill the binning values
            else if (boomResult.front() == "#BINS_W")
                for (int i = 1; i < boomResult.size(); i++)
                    pi0w.at(i-1) = hephelp::StrToDouble(boomResult.at(i));
            else if (boomResult.front() == "#BINS_TPRIME")
                for (int i = 1; i < boomResult.size(); i++)
                    pi0tpr.at(i-1) = hephelp::StrToDouble(boomResult.at(i));
            else if (boomResult.front() == "#BINS_QSQ")
                for (int i = 1; i < boomResult.size(); i++)
                    pi0qsq.at(i-1) = hephelp::StrToDouble(boomResult.at(i));
            else if (boomResult.front() == "#BEGIN_DATA") {
                inData = true;
                //this prepares the arrays from the set ranges above
                resetData();
            }
        }
        //commentary also skip
        else if (line.substr(0,2)=="//")
            continue;
    }

}





void HPionicData::loadFile(string _fileName)
{
    fileName = _fileName;
    //this replaces the format(7F10.4) as we actually just read ALL the decimals
    //and convert them to double because RAM is cheap now and we have 7 columns.
    hephelp::Sfilereader2dim<double>(fileName,pi0data,5);
    fillLegacyArrays();

}






int HPionicData::getBinW(double _w)
{
    int binl = -1;
    double min = 1000;

    for (unsigned int i =0; i < pi0w.size(); i ++)
    {
        if ( min > abs(pi0w[i] - _w)) //maybe do a while-looping here
        {
            min = abs(pi0w[i] - _w);
            binl = i;
        }
    }
    return binl;
}


int HPionicData::getBinQsq(double _qsq)
{
    int binl = -1;
    double min = 1000;
    for (unsigned int i =0; i < pi0qsq.size(); i ++)
    {
        if ( min > abs(pi0qsq[i] - _qsq)) //maybe do a while-looping here
        {
            min = abs(pi0qsq[i] - _qsq);
            binl = i;
        }
    }

    return binl;
}


int HPionicData::getBinTPrim(double _tprim)
{
    int binl = -1;
    double min = 1000;
    for (unsigned int i =0; i < pi0tpr.size(); i ++)
    {
        if ( min > abs(pi0tpr[i] - _tprim)) //maybe do a while-looping here
        {
            min = abs(pi0tpr[i] - _tprim);
            binl = i;
        }
    }
    return binl;
}












void HPionicData::resetData()
{
    pi0data.clear();
    pi0sigl.clear();
    pi0sigt.clear();
    vector<double> TMP;
    for (int i =0; i<5; i++)
        pi0data.push_back(TMP);





    //we initialize them all with these 17 elements
    vector<double> itmp;
    for (int i =0; i < pi0tpr.size(); i++)
        itmp.push_back(0.0);

//     vector<double> ktmp;
    vector< vector<double> > jtmp;

    for (int i=0; i <= pi0w.size(); i++)
    {
        for (int j=0; j <= pi0qsq.size(); j++)
            jtmp.push_back(itmp);
        pi0sigl.push_back(jtmp);
        pi0sigt.push_back(jtmp);
        jtmp.clear();
    }
}

void HPionicData::fillLegacyArrays(void)
{

    //time for some nasty looping - maybe we can think of some restructuring in the future?
    //i dont know :)
    for (int i=0; i < pi0w.size(); i++)
        for (int j=0; j < pi0qsq.size(); j++)
            for (int k=0; k < pi0tpr.size(); k++)
            {
                pi0sigl.at(i).at(j).at(k)=pi0data.at(3).at((i*pi0tpr.size()*pi0qsq.size()+ j*pi0tpr.size() + k));
//                 cout << "filled array with " << pi0data.at(3).at((i*pi0tpr.size()*pi0qsq.size()+ j*pi0tpr.size() + k)) << " with indizes: " <<i << " " << j << " " << k << endl;
                pi0sigt.at(i).at(j).at(k)=pi0data.at(4).at((i*pi0tpr.size()*pi0qsq.size()+ j*pi0tpr.size() + k));
            }
            cout << "Prepared table with binnings in w,qsq,tpr " << pi0w.size() << " " << pi0qsq.size() << " " << pi0tpr.size() << endl;
}
