/*!
 *  \file hbookbackendROOT.h
 *  \date Created on: Jan 31, 2013
 *  \author Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2013 All Right Reserved
 */

#ifndef HBOOKBACKENDROOT_H_
#define HBOOKBACKENDROOT_H_




#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <TH1D.h>
#include <TH2D.h>
#include <TH2F.h>
#include <TFile.h>
#include <TStyle.h>

#include "hbookbackend.h"
#include "hhelper.h"

using namespace std;
/*!
 * \brief A ROOT-based histogram backend-class
 *
 * This will only be built, if ROOT is found by CMAKE during compile-time.
 *
 */


static void ROOTTemplate(void)
{
    gStyle->SetLabelSize(0.06, "xy");
    gStyle->SetTitleSize(0.06, "xy");
    gStyle->SetTitleOffset(0.75, "xy");
    gStyle->SetLabelOffset(0., "X");
    gStyle->SetLabelOffset(0.01, "Y");

    gStyle->SetHistLineColor(1);
    gStyle->SetHistLineWidth(2);
    gStyle->SetHistFillColor(5);
    gStyle->SetHistFillStyle(1001);
}

class HBookBackEndROOT: public HBookBackEnd
{
public:
    HBookBackEndROOT(std::string _nameSuffix):HBookBackEnd(_nameSuffix) {
        ROOTTemplate();
	histCounter = 0;
    };
    ~HBookBackEndROOT();
    /*! \brief adds a double-typed histogram to the list, returns the new index else */
    int addHBook1D(double* _var, double* _weight, int _binCount, double _lowerBound, double _upperBound, string _title);
    /*! \brief NOT IMPLEMENTED YET! (dont know if it ever will be!) */
    int addHBook1I(int* _var, double* _weight, int _binCount, int _lowerBound, int _upperBound, string _title);
    /*! \brief saves the histograms to a file with name _fileName. Extension .root will be added hardcoded to make the backends use different files */
    void dumpToFile(string _fileName);

    /*! \brief add a 2-dim Double-Histogram */
    int addHBook2D(double* _varX, double* _varY , double* _weight, int _binCountX,int _binCountY , double _lowerBoundX,double _lowerBoundY, double _upperBoundX,double _upperBoundY, std::string _title);



    /*! \brief adds a float-typed histogram to the list, returns the new index else */
    int addHBook1F(float* _var, float* _weight, int _binCount, float _lowerBound, float _upperBound, string _title);
    /*! \brief add a 2-dim Double-Histogram */
    int addHBook2F(float* _varX, float* _varY , float* _weight, int _binCountX,int _binCountY , float _lowerBoundX,float _lowerBoundY, float _upperBoundX,float _upperBoundY, std::string _title);


    /*! \brief fills all histograms */
    void fill();
    /*!  \brief fills only histogram # _histNum */
    int fill(int _histNum);
    /*! \brief returns title of histogram if _num is in range, returns "_ERROR_OUT_OF_RANGE" else */
    string getTitle(int _num);

private:

    int histCounter;
    vector<TH1D*>  histList;
    vector<TH2D*> histList2;
    vector<double*> varList;
    vector<double*> weightList;

    vector<double*> varListX;
    vector<double*> varListY;
    vector<double*> weightList2D;

    vector<TH1F*>  histListF;
    vector<TH2F*> histList2F;
    vector<float*> varListF;
    vector<float*> weightListF;

    vector<float*> varListXF;
    vector<float*> varListYF;
    vector<float*> weightList2DF;
};


#endif
