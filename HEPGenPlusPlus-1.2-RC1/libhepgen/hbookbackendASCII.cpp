#include "hbookbackendASCII.h"



HBookBackEndASCII::~HBookBackEndASCII()
{
    for (unsigned int i =0; i<histList.size(); i++)
        delete histList.at(i);
}

int HBookBackEndASCII::fill(int _histNum)
{
    if (_histNum < 0 || _histNum >= static_cast<int>(histList.size()))
        return -2;

    return histList.at(_histNum)->fill();
}


void HBookBackEndASCII::fill()
{
    for (unsigned int i =0; i<histList.size(); i++)
        histList.at(i)->fill();
}

int HBookBackEndASCII::addHBook2D(double* _varX, double* _varY, double* _weight, int _binCountX, int _binCountY, double _lowerBoundX, double _lowerBoundY, double _upperBoundX, double _upperBoundY, string _title)
{
    ASCIIhist2d* tmp = new ASCIIhist2d(_varX,_varY, _weight, _binCountX, _binCountY,_lowerBoundX,_lowerBoundY,_upperBoundX,_upperBoundY,_title);
    histList.push_back(tmp);
    return histList.size()-1;
}

int HBookBackEndASCII::addHBook2F(float* _varX, float* _varY, float* _weight, int _binCountX, int _binCountY, float _lowerBoundX, float _lowerBoundY, float _upperBoundX, float _upperBoundY, string _title)
{
    ASCIIhist2f* tmp = new ASCIIhist2f(_varX,_varY, _weight, _binCountX, _binCountY,_lowerBoundX,_lowerBoundY,_upperBoundX,_upperBoundY,_title);
    histList.push_back(tmp);
    return histList.size()-1;
}



void HBookBackEndASCII::dumpToFile(string _fileName)
{
    string finalFile = _fileName+".dat";
    ofstream outfile(finalFile.c_str(), ios::out);
    outfile<<"# HEPGEN_C-DUMP" << endl;
    outfile<<"# Number of Histograms: " << histList.size() << endl;
    outfile<<"# Beginning dump! " << endl;
    for (unsigned int i = 0; i < histList.size(); i++)
    {
        outfile<<"# Histogram number: " << i << endl;
        histList.at(i)->dumpToFile(outfile);
        outfile<< endl << endl << endl << flush;

    }

}

string HBookBackEndASCII::getTitle(int _num)
{
    if (_num <0 || _num > static_cast<int>(histList.size())-1)
        return "_ERROR_OUT_OF_RANGE";
    else
        return histList.at(_num)->getTitle();

}




int HBookBackEndASCII::addHBook1F(float* _var, float* _weight, int _binCount, float _lowerBound, float _upperBound, string _title)
{
    ASCIIhistf* tmp = new ASCIIhistf(_var, _weight, _binCount,_lowerBound,_upperBound,_title);
    histList.push_back(tmp);
    return histList.size()-1;
}

int HBookBackEndASCII::addHBook1D(double* _var, double* _weight, int _binCount, double _lowerBound, double _upperBound, string _title)
{
    ASCIIhistd* tmp = new ASCIIhistd(_var, _weight, _binCount,_lowerBound,_upperBound,_title);
    histList.push_back(tmp);
    return histList.size()-1;
}

int HBookBackEndASCII::addHBook1I(int* _var, double* _weight, int _binCount, int _lowerBound, int _upperBound, string _title)
{
    cout << "STUB! HBookBackEndASCII::addHBOOK1I is not yet implemented! DO NOT USE!!!!";
    return -2;
}

ASCIIhistd::ASCIIhistd(double* _var, double* _weight, int _binCount, double _lowerBound, double _upperBound, string _title)
{
    reset(_var, _weight, _binCount, _lowerBound, _upperBound, _title);
}

int ASCIIhistd::fill()
{
    fillCount++;
    double value = *variable;
    //check for over or underflow
    if (value < lowerBound)
    {
        histData[0]+= *weight;
        return 0;
    }
    if (value > upperBound)
    {
        histData[binCount+1]+= *weight;
        return -1;
    }
    //if we are still here we are in a the range!
    int binNum = int((value - lowerBound)/stepLength);
    histData[binNum]+= *weight;

    return binNum;
}

void ASCIIhistd::reset(double* _var, double* _weight, int _binCount, double _lowerBound, double _upperBound, string _title)
{
    weight = _weight;
    variable = _var;
    binCount = _binCount;
    lowerBound = _lowerBound;
    upperBound = _upperBound;
    fillCount=0;
    stepLength = (upperBound - lowerBound)/binCount;
    title = _title;
    histData.clear();
    //double dist = (upperBound-lowerBound)/binCount;

    //we need 2 more because of over and underflow bins
    for (int i = 0; i < binCount+2; i++)
    {
        histData[i]=0;
    }

}


void ASCIIhistd::dumpToFile(ofstream& _file)
{
    _file << "HIST" << endl << "TITLE " << title << endl;
    _file << "LOW_X " << lowerBound << endl << "HIGH_X " << upperBound << endl << "NBINS_X" << binCount << endl;
    for (int i=0; i<binCount+2; i++)
        _file<<i<<" " << histData[i] << endl;
    _file << "END HIST" << endl;
}




ASCIIhist2d::ASCIIhist2d(double* _varX, double* _varY, double* _weight, int _binCountX, int _binCountY, double _lowerBoundX, double _lowerBoundY, double _upperBoundX, double _upperBoundY, string _title)
{
    weight = _weight;
    variableX = _varX;
    variableY = _varY;
    binCountX = _binCountX;
    binCountY = _binCountY;
    lowerBoundX = _lowerBoundX;
    lowerBoundY = _lowerBoundY;
    upperBoundX = _upperBoundX;
    upperBoundY = _upperBoundY;
    title = _title;
    reset();
}








void ASCIIhist2d::reset()
{
    fillCount=0;
    stepLengthX = (upperBoundX - lowerBoundX)/binCountX;
    stepLengthY = (upperBoundY - lowerBoundY)/binCountY;
    histData.clear();
    histData.resize(binCountY+2);

    for (int i = 0; i < binCountY+2; i++)
    {
        histData.at(i).resize(binCountX+2);
        memset(&histData.at(i)[0],0,sizeof(double)*(binCountX+2));
    }

}

void ASCIIhist2d::dumpToFile(ofstream& _file)
{
    _file << "HIST-2D" << endl << "TITLE " << title << endl;
    _file << "LOW_X " << lowerBoundX << endl << "HIGH_X " << upperBoundX << endl << "NBINS_X" << binCountX << endl;
    _file << "LOW_Y " << lowerBoundY << endl << "HIGH_Y " << upperBoundY << endl << "NBINS_Y" << binCountY << endl;

    for (int a=0; a<binCountY+2; a++)
    {
        _file << a << ": ";
        for (int i=0; i<binCountX+2; i++)
            _file<<" " << histData.at(a).at(i);
        _file << endl;
    }
    _file << "END HIST" << endl;
}



int ASCIIhist2d::fill()
{
    fillCount++;
    double valueX = *variableX;
    double valueY = *variableY;

    int binX = int( (valueX - lowerBoundX)/stepLengthX);
    int binY = int( (valueY - lowerBoundY)/stepLengthY);



    //check for over or underflow
    if (valueX < lowerBoundX)
    {
        binX = 0;
    }
    else if (valueX > upperBoundX)
    {
        binX = binCountX+1;
    }
    if (valueY < lowerBoundY)
    {
        binY = 0;
    }
    else if (valueY > upperBoundY)
    {
        binY = binCountY+1;
    }
    //if we are still here we are in a the range!
    histData.at(binY).at(binX) += *weight;

    return 0;
}


//------------------------------- float type 


ASCIIhistf::ASCIIhistf(float* _var, float* _weight, int _binCount, float _lowerBound, float _upperBound, string _title)
{
    reset(_var, _weight, _binCount, _lowerBound, _upperBound, _title);
}

int ASCIIhistf::fill()
{
    fillCount++;
    float value = *variable;
    //check for over or underflow
    if (value < lowerBound)
    {
        histData[0]+= *weight;
        return 0;
    }
    if (value > upperBound)
    {
        histData[binCount+1]+= *weight;
        return -1;
    }
    //if we are still here we are in a the range!
    int binNum = int((value - lowerBound)/stepLength);
    histData[binNum]+= *weight;

    return binNum;
}

void ASCIIhistf::reset(float* _var, float* _weight, int _binCount, float _lowerBound, float _upperBound, string _title)
{
    weight = _weight;
    variable = _var;
    binCount = _binCount;
    lowerBound = _lowerBound;
    upperBound = _upperBound;
    fillCount=0;
    stepLength = (upperBound - lowerBound)/binCount;
    title = _title;
    histData.clear();
    //float dist = (upperBound-lowerBound)/binCount;

    //we need 2 more because of over and underflow bins
    for (int i = 0; i < binCount+2; i++)
    {
        histData[i]=0;
    }

}


void ASCIIhistf::dumpToFile(ofstream& _file)
{
    _file << "HIST" << endl << "TITLE " << title << endl;
    _file << "LOW_X " << lowerBound << endl << "HIGH_X " << upperBound << endl << "NBINS_X" << binCount << endl;
    for (int i=0; i<binCount+2; i++)
        _file<<i<<" " << histData[i] << endl;
    _file << "END HIST" << endl;
}




ASCIIhist2f::ASCIIhist2f(float* _varX, float* _varY, float* _weight, int _binCountX, int _binCountY, float _lowerBoundX, float _lowerBoundY, float _upperBoundX, float _upperBoundY, string _title)
{
    weight = _weight;
    variableX = _varX;
    variableY = _varY;
    binCountX = _binCountX;
    binCountY = _binCountY;
    lowerBoundX = _lowerBoundX;
    lowerBoundY = _lowerBoundY;
    upperBoundX = _upperBoundX;
    upperBoundY = _upperBoundY;
    title = _title;
    reset();
}








void ASCIIhist2f::reset()
{
    fillCount=0;
    stepLengthX = (upperBoundX - lowerBoundX)/binCountX;
    stepLengthY = (upperBoundY - lowerBoundY)/binCountY;
    histData.clear();
    histData.resize(binCountY+2);

    for (int i = 0; i < binCountY+2; i++)
    {
        histData.at(i).resize(binCountX+2);
        memset(&histData.at(i)[0],0,sizeof(float)*(binCountX+2));
    }

}

void ASCIIhist2f::dumpToFile(ofstream& _file)
{
    _file << "HIST-2D" << endl << "TITLE " << title << endl;
    _file << "LOW_X " << lowerBoundX << endl << "HIGH_X " << upperBoundX << endl << "NBINS_X" << binCountX << endl;
    _file << "LOW_Y " << lowerBoundY << endl << "HIGH_Y " << upperBoundY << endl << "NBINS_Y" << binCountY << endl;

    
    for (int a=0; a<binCountY+2; a++)
    {
        _file << a << ": ";
        for (int i=0; i<binCountX+2; i++)
            _file<<" " << histData.at(a).at(i);
        _file << endl;
    }
    _file << "END HIST" << endl;
}



int ASCIIhist2f::fill()
{
    fillCount++;
    float valueX = *variableX;
    float valueY = *variableY;

    int binX = int( (valueX - lowerBoundX)/stepLengthX);
    int binY = int( (valueY - lowerBoundY)/stepLengthY);



    //check for over or underflow
    if (valueX < lowerBoundX)
    {
        binX = 0;
    }
    else if (valueX > upperBoundX)
    {
        binX = binCountX+1;
    }
    if (valueY < lowerBoundY)
    {
        binY = 0;
    }
    else if (valueY > upperBoundY)
    {
        binY = binCountY+1;
    }
    //if we are still here we are in a the range!
    histData.at(binY).at(binX) += *weight;

    return 0;
}









