
#include <iostream>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fstream>

#include <string>
#include <sys/stat.h>
#include <unistd.h>

#include "libGKPi0.h"
#include "gkSubProcessTable.h"

bool exists_test (const std::string& name) {
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}

int main (int argc, char** argv){   
  if (argc < 4){
    printf("Usage: %s [qsq] [xbj] [t] [phi] (-I [interpolationCacheIndex])\n",argv[0]);
    printf("Alternate Usage: %s [qsq] [W] [t] [phi] -a (-I [interpolationCacheIndex])\n",argv[0]);    
    printf("Give phi in radians ;) \n");
    return -1;
  }

  double xbj,qsq,w2,t;
  qsq=atof(argv[1]);
  bool useInterPolation = false;
  double phi = atof(argv[4]);
  string fileNameInterpolation;
  if (argc == 6 || argc == 8){
    w2=atof(argv[2]);
//     w2 = pow(w2,2.0);
     xbj=GKPI0::getXbj(qsq,w2);
//      w2 = pow(w2,2.0);  
  }
  else if ( argc == 5 || argc == 7 ){
    xbj=atof(argv[2]);
    w2=GKPI0::getWsq(qsq,xbj);
  }
  
  if (argc > 6){
    useInterPolation = true;
    fileNameInterpolation = argv[argc-1];
  }
    
  
  double m=0.93827203;

  GKPI0::amplitude myAmp;
  t=atof(argv[3]);
  double xi=GKPI0::compassxi(xbj);
  double t0 = -4*pow(m,2.0)*pow(xi,2.0)/(1.-pow(xi,2.0));
  
  double tprime = t-t0;
  printf("Final parameters: Q^2 %.4f, W %.4f, xbj %.4e, xi %.4e,t %.4f t' %.4f, phi %.4f\n",qsq,w2,xbj,xi,t,tprime,phi);
  if (!useInterPolation){

  char buffer[200];
  //check if file with buffered subprocess amplitudes exist
  std::string myPath="";
  if (getenv("HEPGEN_CREATESPA") != NULL)
   myPath = getenv("HEPGEN_CREATESPA");
  else
    myPath ="./";
  sprintf(buffer,"%s/preparation_%.4f_%.4f.dat",myPath.c_str(),qsq,w2);
  if (exists_test(buffer)) // if so, load it
    GKPI0::loadPreparationFromFile(buffer,qsq,xi);
  else{ //else we calculate them new
    GKPI0::prepareConvolution(qsq,xi);
    GKPI0::savePreparationToFile(buffer); //and we save them, for the next run
  }
  //finally calculate amplitudes
  myAmp =  GKPI0::getAmplitude(qsq,xi,xbj,t);
  }
  else{
    gkSubProcessTableCache myCache;
    myCache.loadCache(fileNameInterpolation);
    myAmp = myCache.getAmpsForKine(qsq,xbj,t);
  }
  printf("Amplituden: \n");
  std::cout << "Mmmp0 Mmpp0 Mpmp0 Mppp0" << std::endl;
  
  std::cout << myAmp.Mmmp0 << " " << myAmp.Mmpp0 << " " << myAmp.Mpmp0 << " " << myAmp.Mppp0 << std::endl;
  
  std::cout << "Twist-2 contribution! \n M0pp M0mp \n " << myAmp.M0pp << " " << myAmp.M0mp << "\n";
  
  
  
  
  //and in the end, calculate the cross section
  double sigma = GKPI0::getCX(myAmp,w2);
  double sigmaTT = GKPI0::getCXTT(myAmp,w2,phi);
  double sigmaL = GKPI0::getCXL(myAmp,w2);
  double sigmaLT = GKPI0::getCXLT(myAmp,w2,phi);
  
  
  printf("Final: %.4f %.4f %.4f dSigmaTot/dt = %.4f + e * %.4f + e* %.4f + sqrt(2e(1+e)) * %.4f\n",qsq,xbj,t,sigma,sigmaL,sigmaTT,sigmaLT);
  return 0;
}
