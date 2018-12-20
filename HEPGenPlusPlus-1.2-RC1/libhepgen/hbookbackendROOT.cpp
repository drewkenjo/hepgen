#include "hbookbackendROOT.h"


int HBookBackEndROOT::addHBook1D(double* _var, double* _weight, int _binCount, double _lowerBound, double _upperBound, string _title)
{
    string name = "histogram_"+hephelp::IntToStr(histCounter++)+nameSuffix;
    TH1D* tmp = new TH1D(name.c_str(), _title.c_str(), _binCount,_lowerBound,_upperBound);
    histList.push_back(tmp);
    varList.push_back(_var);
    weightList.push_back(_weight);
    return histList.size()-1;
}

int HBookBackEndROOT::addHBook1F(float* _var, float* _weight, int _binCount, float _lowerBound, float _upperBound, string _title)
{
    string name = "histogram_"+hephelp::IntToStr(histCounter++)+nameSuffix;
    TH1F* tmp = new TH1F(name.c_str(), _title.c_str(), _binCount,_lowerBound,_upperBound);
    histListF.push_back(tmp);
    varListF.push_back(_var);
    weightListF.push_back(_weight);
    return histListF.size()-1;
}


int HBookBackEndROOT::addHBook1I(int* _var, double* _weight, int _binCount, int _lowerBound, int _upperBound, string _title)
{
    cout << "ERROR! HBookBackEndROOT::addHBook1I IS A STUB! NOT YET IMPLEMENTED!!!" << endl;
    return -2;
}


void HBookBackEndROOT::dumpToFile(string _fileName)
{
     string finalFile = _fileName+".root";
     TFile* outfile = new TFile(finalFile.c_str(),"RECREATE");
     for (unsigned int i =0; i < histList.size(); i++)
         histList.at(i)->Write();
     for (unsigned int i = 0; i < histList2.size(); i++)
         histList2.at(i)->Write();
     for (unsigned int i =0; i < histListF.size(); i++)
         histListF.at(i)->Write();
     for (unsigned int i = 0; i < histList2F.size(); i++)
         histList2F.at(i)->Write();
     
    
    outfile->Write();
//     outfile->Flush();
    outfile->Close();
    delete outfile;
    
}


void HBookBackEndROOT::fill()
{
    for (unsigned int i = 0; i < histList.size(); i++)
        histList.at(i)->Fill(*varList.at(i),*weightList.at(i));

    for (unsigned int i = 0; i < histList2.size(); i++)
        histList2.at(i)->Fill(*varListX.at(i),*varListY.at(i),*weightList2D.at(i));
    
    for (unsigned int i = 0; i < histListF.size(); i++)
        histListF.at(i)->Fill(*varListF.at(i),*weightListF.at(i));

    for (unsigned int i = 0; i < histList2F.size(); i++)
        histList2F.at(i)->Fill(*varListXF.at(i),*varListYF.at(i),*weightList2DF.at(i));
    
    

}

int HBookBackEndROOT::fill(int _histNum)
{
    if (_histNum < 0 || _histNum >= static_cast<int>(histList.size()+histList2.size()))
        return -2;
    if (_histNum < histList.size())
        return histList.at(_histNum)->Fill(*varList.at(_histNum),*weightList.at(_histNum));
    else
        return histList2.at(_histNum%histList.size())->Fill(*varListX.at(_histNum%histList.size()),*varListY.at(_histNum%histList.size()),*weightList2D.at(_histNum%histList.size()));
}


HBookBackEndROOT::~HBookBackEndROOT()
{
    for (unsigned int i = 0; i < histList.size(); i++)
        delete histList.at(i);
    for (unsigned int i = 0; i < histList2.size(); i++)
        delete histList2.at(i);

}


string HBookBackEndROOT::getTitle(int _num)
{
    if (_num < 0 || _num >= static_cast<int>(histList.size()))
        return "_ERROR_OUT_OF_RANGE";
    return string(histList.at(_num)->GetTitle());
}

int HBookBackEndROOT::addHBook2D(double* _varX, double* _varY, double* _weight, int _binCountX, int _binCountY, double _lowerBoundX, double _lowerBoundY, double _upperBoundX, double _upperBoundY, string _title)
{
    string name = "histogram_"+hephelp::IntToStr(histCounter++)+nameSuffix;
    TH2D* tmp = new TH2D(name.c_str(), _title.c_str(), _binCountX,_lowerBoundX,_upperBoundX,_binCountY,_lowerBoundY,_upperBoundY);
    histList2.push_back(tmp);
    varListX.push_back(_varX);
    varListY.push_back(_varY);
    weightList2D.push_back(_weight);
    return histList2.size()-1;
}

int HBookBackEndROOT::addHBook2F(float* _varX, float* _varY, float* _weight, int _binCountX, int _binCountY, float _lowerBoundX, float _lowerBoundY, float _upperBoundX, float _upperBoundY, string _title)
{
    string name = "histogram_"+hephelp::IntToStr(histCounter++)+nameSuffix;
    TH2F* tmp = new TH2F(name.c_str(), _title.c_str(), _binCountX,_lowerBoundX,_upperBoundX,_binCountY,_lowerBoundY,_upperBoundY);
    histList2F.push_back(tmp);
    varListXF.push_back(_varX);
    varListYF.push_back(_varY);
    weightList2DF.push_back(_weight);
    return histList2F.size()-1;
}
