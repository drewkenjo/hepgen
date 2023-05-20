#include <iostream>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <vector>

#include <string>
#include <sys/stat.h>
#include <unistd.h>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "libGKPi0.h"
#include "gkSubProcessTable.h"
#include "TString.h"
#include "TDatime.h"

struct prepvecs {
  std::vector<double> weights;
  std::vector<double> xpseudoList;
  std::vector<double> xlow;
  std::vector<double> xhigh;
  std::vector<TComplex> xlowResult;
  std::vector<TComplex> xhighResult;
  std::vector<TComplex> xlowResultTwist2;
  std::vector<TComplex> xhighResultTwist2;
};


struct kinpoint {
  double qsq;
  double xbj;
  double ww;
  double xi;
  double tt;
  double s0;
  double ds0;
  double sTT;
  double dsTT;
  double sLT;
  double dsLT;
  double E0;
};


prepvecs readFile(std::string fileName) {
  prepvecs vecs;

  std::ifstream myInFile;
  myInFile.open(fileName,std::ios::in);

  if (!myInFile.good()) {
    std::cout<<"file is not found "<<fileName<<std::endl;
    exit(111);
  }

  while (!myInFile.eof()) {
    double buf,buf2;
    myInFile >> buf;
    vecs.xlow.push_back(buf);
    myInFile >> buf;
    vecs.xlowResult.push_back(TComplex(buf,0.0));
    myInFile >> buf;
    vecs.xhigh.push_back(buf);
    myInFile >> buf;
    myInFile >> buf2;
    vecs.xhighResult.push_back(TComplex(buf,buf2));
    myInFile >> buf;
    vecs.xlowResultTwist2.push_back(TComplex(buf,0.0));
    myInFile >> buf;
    myInFile >> buf2;
    vecs.xhighResultTwist2.push_back(TComplex(buf,buf2));
    myInFile>> buf;
    vecs.weights.push_back(buf);
    myInFile>> buf;
    vecs.xpseudoList.push_back(buf);
  }

  return vecs;
}



std::vector<kinpoint> points[3];
std::vector<prepvecs> preps[3];

double chisquare(const double *xx) {
  GKPI0::SetETbarUNorm(xx[0]);
  GKPI0::SetETbarDNorm(xx[1]);

  GKPI0::SetETbarUtSlope(xx[2]);
  GKPI0::SetETbarDtSlope(xx[3]);

  GKPI0::SetAlpha0U(xx[4]);
  GKPI0::SetAlphaStrU(xx[5]);
  GKPI0::SetAlpha0D(xx[6]);
  GKPI0::SetAlphaStrD(xx[7]);

//  GKPI0::SetETbarUCoefs(xx[5], xx[6], xx[7]);
//  GKPI0::SetETbarDCoefs(xx[8], xx[9], xx[10], xx[11], xx[12]);

  double chi2 = 0;
  for(int ich=0;ich<3;ich++){
    GKPI0::SetReactionPar(1+ich);

    for(int ip=0;ip<points[ich].size();ip++) {

      GKPI0::loadPreparationFromRam(points[ich][ip].xbj, preps[ich][ip].xlowResult, preps[ich][ip].xhighResult,
         preps[ich][ip].xlowResultTwist2, preps[ich][ip].xhighResultTwist2,
         preps[ich][ip].xlow, preps[ich][ip].xhigh, preps[ich][ip].weights, preps[ich][ip].xpseudoList);

      GKPI0::amplitude myAmp =  GKPI0::getAmplitude(points[ich][ip].qsq, points[ich][ip].xi, points[ich][ip].xbj, points[ich][ip].tt);
      double sT = GKPI0::getCX(myAmp,points[ich][ip].ww);
      double s0 = sT;

      if(points[ich][ip].E0>0) {
        double eps = GKPI0::getEpsilon(points[ich][ip].qsq, points[ich][ip].xbj, points[ich][ip].E0);
        double sL = GKPI0::getCXL(myAmp,points[ich][ip].ww);
        s0 += eps*sL;
      }

      chi2 += pow((points[ich][ip].s0 - s0)/points[ich][ip].ds0, 2);
//      chi2 += pow((points[ich][ip].sLT - sLT)/points[ich][ip].dsLT, 2);
//      chi2 += pow((points[ich][ip].sTT - sTT)/points[ich][ip].dsTT, 2);
    }
  }

  return chi2;
}


int main (int argc, char** argv) {
  if (getenv("HEPGEN") == NULL) {
    printf("$HEPGEN env is not set, do not know where to load grid table values from!\n");
    return 11;
  }

  std::string hepPath = getenv("HEPGEN");
  TString status("oldmix_200");
  bool minos = false;

  TString prepdirs[] = {"pi0P/mu.eq.2", "etaP/mu.eq.2", "pi0N/mu.eq.2"};
  for(int iarg=1;iarg<argc;iarg++) {
    TString arg(argv[iarg]);
    int ich = 0;
    double E0 = 5.75;

    if(arg.EqualTo("-newmix")) {
      GKPI0::set_new_eta_mixing_angle(true);
      status.ReplaceAll("oldmix", "newmix");
      continue;

    } else if(arg.EqualTo("-minos")) {
      minos = true;
      continue;

    } else if(arg.EqualTo("-176")) {
      GKPI0::set_mu_eta(1.76);
      prepdirs[1] = "etaP/mu.eq.1.76";
      status.ReplaceAll("200", "176");
      continue;

    } else if(arg.Contains("neutron")) {
      ich = 2;

    } else if(arg.Contains("eta")) {
      ich = 1;

    } else {
      ich = 0;
    }

    if(arg.Contains("hallA.pi0.y21")) E0 = 10;
    else if(arg.Contains("compass")) E0 = 160;

    
    TString prepdir(prepdirs[ich]);
    GKPI0::SetReactionPar(1+ich);

    std::ifstream ff(argv[iarg]);
    while(ff.good()) {
      double qnom,xnom,qq,xb,tt,s0,ds0,dds0,slt,dslt,ddslt,stt,dstt,ddstt;
      ff>>qnom>>xnom>>qq>>xb>>tt>>s0>>ds0>>dds0>>slt>>dslt>>ddslt>>stt>>dstt>>ddstt;

      if(ff.good()) {
        double ww = GKPI0::getWsq(qq, xb);
        double xi = GKPI0::compassxi(qq, xb);
        double tmin = GKPI0::getTmin(qq, xb);

        TString prepname(Form("preparation/%s/preparation_%.4f_%.4f.dat",prepdir.Data(),qq,ww));
        std::ifstream fprep(prepname.Data());
        if(!fprep.good()) {
          std::cerr<<"file "<<prepname<<" does NOT exist!!!!"<<std::endl;
          return 111;
        }

        kinpoint pp = {qq, xb, ww, xi, -tt, s0, ds0, stt, dstt, slt, dslt, E0};
        points[ich].push_back(pp);
        preps[ich].push_back(readFile(prepname.Data()));
      }
    }
    ff.close();
  }

  std::cout<<status<<std::endl;

  std::cout<<"npoints: "<<points[0].size() + points[1].size() + points[2].size()<<std::endl;
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
//  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("GSLSimAn", "");
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  min->SetMaxIterations(10000);  // for GSL
  min->SetTolerance(0.05);
  min->SetPrintLevel(20);


  int npars = 8;
  ROOT::Math::Functor fchi2(&chisquare, npars);
  min->SetFunction(fchi2);


  //double cc[] = {YYYYYY};
  //double cc[] = {4.83, 3.57, 0.5, 0.5, 0.3, 0.45, 0.3, 0.45};
  double cc[] = {2.0747, 1.3451, 0.5, 0.5, 0.3, 0.45, 0.3, 0.45};

  double etNu = cc[0], etNd = cc[1];
  double etbu = cc[2], etbd = cc[3];
  double alpha0U = cc[4], alphastrU = cc[5];
  double alpha0D = cc[6], alphastrD = cc[7];

/*
  etNu = 2.17106;
  etNd = 0.721776;
  etbu = -0.599638;;
  etbd = 2.4842;
  alphastr  = 0.488117;
*/

  min->SetVariable(0, "Nu", etNu, 0.01);
  min->SetVariable(1, "Nd", etNd, 0.01);
  min->SetVariable(2, "bu", etbu, 0.01);
  min->SetVariable(3, "bd", etbd, 0.01);
  min->SetVariable(4, "alpha0U", alpha0U, 0.01);
  min->SetVariable(5, "alphastrU", alphastrU, 0.01);
  min->SetVariable(6, "alpha0D", alpha0D, 0.01);
  min->SetVariable(7, "alphastrD", alphastrD, 0.01);

//  min->SetVariableLowerLimit(2,0);
//  min->SetVariableLowerLimit(3,0);
//  min->FixVariable(4);

//  double cc[] = {1, 0, -1, 1, 0, -2, 0, 1};
//  double cc[] = {1.17336, -0.764575, -2.58361, 1.44918, 0.787773, -2.58634, -0.938394, 0.565228};
//  double cc[] = {0.789148,3.49199,-1.20457,12.3844,-337.816,2947.18,-6704.85,5311.65};
//  for(int ic=0;ic<8;ic++)
//    min->SetVariable(ic+5, Form("c%d",ic), cc[ic], 0.01);

  min->Minimize();
  const double *xs = min->X();
  const double *errs = min->Errors();
  std::cout << "Minimum: f(): " << min->MinValue()  << std::endl;

  int nsize = points[0].size() + points[1].size() + points[2].size();
  status.Append(Form(" %f %d", min->MinValue(), nsize));

  for(int ic=0;ic<npars;ic++) {
    double errLow = errs[ic], errUp = errs[ic];
    if(minos) min->GetMinosError(ic, errLow, errUp);
    std::cout<<"fitpi0_parameters: "<<xs[ic]<<" "<<errLow<<" "<<errUp<<std::endl;
    status.Append(Form(" %f %f %f", xs[ic], errLow, errUp));
  }
  std::cout<<"fitpi0_npoints: "<<nsize<<std::endl;

  TDatime out;
  TString suf(status(0,status.First(" ")));
  auto outff = fopen(Form("%d%02d%02d_%06d_%s_%d.fitoutput", out.GetYear(), out.GetMonth(), out.GetDay(), out.GetTime(), suf.Data(), nsize), "w");
  status.Puts(outff);
  fclose(outff);

  /*
  double prevq=0, prevx=0;
  for(int ip=0;ip<points.size();ip++) {
    if(prevq == points[ip].qnom && prevx == points[ip].xnom)
      continue;

    prevq = points[ip].qnom;
    prevx = points[ip].xnom;

    for(double tt=0; tt>-2; tt-=0.1) {
      if(tt<tmin) {
        GKPI0::amplitude myAmp = myCache.getAmpsForKine(points[ip].qnom, points[ip].xnom, tt);
        double sL = GKPI0::getCXL(myAmp,points[ip].wnom);
        double sT = GKPI0::getCX(myAmp,points[ip].wnom);
        double sTT = GKPI0::getCXTT(myAmp,points[ip].wnom, 0);
        double sLT = GKPI0::getCXLT(myAmp,points[ip].wnom, 0);
        double s0 = sT + eps*sL;

        printf("%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",points[ip].qnom, points[ip].xnom, tt,s0,sT,sL,sTT,sLT);
      }
    }
  }
  */

  return 0;
}
