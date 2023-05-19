#include <iostream>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fstream>

#include <sys/stat.h>
#include <unistd.h>

#include <TString.h>
#include "libGKPi0.h"
#include "gkSubProcessTable.h"

bool exists_test (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

int main (int argc, char** argv) {
  string fileNameInterpolation;
  if(exists_test(argv[argc-1]) && TString(argv[argc-1]).Contains("tableIndex")) {
    fileNameInterpolation = argv[argc-1];
    argc--;
  }

  std::vector<double> kins;
  for(int iarg=1;iarg<argc;iarg++) {
    TString arg(argv[iarg]);
    if(arg.Contains("-newmix")) {
      GKPI0::set_new_eta_mixing_angle(true);
    } else if(arg.Contains("-176")) {
      GKPI0::set_mu_eta(1.76);
    } else if(arg.Contains("-eta")) {
      GKPI0::SetReactionPar(2);
    } else if(arg.Contains("-neutron")) {
      GKPI0::SetReactionPar(3);
    } else {
      kins.push_back(atof(argv[iarg]));
    }
  }

  double qsq, w2, xbj, phi=0, tt=-2;
  if(kins.size()>=2) {
    qsq = kins[0];
    xbj = kins[1];
    w2 = kins[1];

    if(kins.size()==3) {
      tt = kins[2];
    } else if(kins.size()==4) {
      phi = kins[3];
    }
  } else {
    std::cerr<<Form("Usage: %s [-newmix] [-176] [-eta] [-neutron] qsq xbj [t] [phi] [interpolationCacheIndex]\n",argv[0]);
    std::cerr<<Form("Give phi in radians ;) \n");
    return -1;
  }

  if(w2<1) {
    w2 = GKPI0::getWsq(qsq,xbj);
  } else {
    xbj = GKPI0::getXbj(qsq,w2);
  }

  GKPI0::amplitude myAmp;
  double xi = GKPI0::compassxi(qsq, xbj);
  double tmin = GKPI0::getTmin(qsq, xbj);

  double tprime = tt-tmin;
  printf("Final parameters: Q^2 %.4f, W %.4f, xbj %.4e, xi %.4e,t %.4f t' %.4f, phi %.4f\n",qsq,w2,xbj,xi,tt,tprime,phi);
  if (fileNameInterpolation.empty()) {
    //check if file with buffered subprocess amplitudes exist
    std::string myPath="";
    if (getenv("HEPGEN_CREATESPA") != NULL)
     myPath = getenv("HEPGEN_CREATESPA");
    else
      myPath ="./";

    TString fprep(Form("%s/preparation_%.4f_%.4f.dat",myPath.c_str(),qsq,w2));

    if (exists_test(fprep.Data())) // if so, load it
      GKPI0::loadPreparationFromFile(fprep.Data(),qsq,xi);
    else { //else we calculate them new
      GKPI0::prepareConvolution(qsq,xi);
      GKPI0::savePreparationToFile(fprep.Data()); //and we save them, for the next run
    }

    //finally calculate amplitudes
    myAmp =  GKPI0::getAmplitude(qsq,xi,xbj,tt);
  } else {

    gkSubProcessTableCache myCache;
    myCache.loadCache(fileNameInterpolation);
    myAmp = myCache.getAmpsForKine(qsq,xbj,tt);
  }


  printf("Amplituden: \n");
  std::cout << "Mmmp0 Mmpp0 Mpmp0 Mppp0" << std::endl;

  std::cout << myAmp.Mmmp0 << " " << myAmp.Mmpp0 << " " << myAmp.Mpmp0 << " " << myAmp.Mppp0 << std::endl;

  std::cout << "Twist-2 contribution! \n M0pp M0mp \n " << myAmp.M0pp << " " << myAmp.M0mp << "\n";


  //and in the end, calculate the cross section
  double sigmaT = GKPI0::getCX(myAmp,w2);
  double sigmaTT = GKPI0::getCXTT(myAmp,w2,phi);
  double sigmaL = GKPI0::getCXL(myAmp,w2);
  double sigmaLT = GKPI0::getCXLT(myAmp,w2,phi);

  double E0 = 5.75;
  double eps = GKPI0::getEpsilon(qsq, xbj, E0);
  double sigma0 = sigmaT + eps*sigmaL;

  printf("Final: %.4f %.4f %.4f dSigmaTot/dt = %.4f + e * %.4f + e* %.4f + sqrt(2e(1+e)) * %.4f\n",qsq,xbj,tt,sigmaT,sigmaL,sigmaTT,sigmaLT);
  printf("Output: %.4f %.4f %.8f %.4f %.4f %.4f %.4f %.4f\n",qsq,xbj,tt,sigma0,sigmaT,sigmaL,sigmaTT,sigmaLT);

  return 0;
}
