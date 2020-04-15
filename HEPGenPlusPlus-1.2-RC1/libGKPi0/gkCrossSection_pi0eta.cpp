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

int main (int argc, char** argv) {
  if (argc < 4) {
    printf("Usage: %s [qsq] [xbj] [t] [phi] (-I [interpolationCacheIndex])\n",argv[0]);
    printf("Alternate Usage: %s [qsq] [W] [t] [phi] -a (-I [interpolationCacheIndex])\n",argv[0]);
    printf("Give phi in radians ;) \n");
    return -1;
  }

  gkSubProcessTableCache myCache[2];
  for(int ich=0;ich<2;ich++)
    myCache[ich].loadCache(argv[argc-2+ich]);
  
  std::ifstream ff(argv[1]);
  while(ff.good()) {
    double qsq, w2, t, phi;
    ff>>qsq>>w2>>t>>phi;
    if(ff.good()) {
      double xbj;
      if (w2<1) {
        xbj = w2;
        w2 = GKPI0::getWsq(qsq,xbj);
      } else {
         xbj = GKPI0::getXbj(qsq,w2);
      }

      /*
      double xi=GKPI0::compassxi(xbj);
      printf("Final parameters: Q^2 %.4f, W %.4f, xbj %.4e, xi %.4e,t %.4f t' %.4f, phi %.4f\n",qsq,w2,xbj,xi,t,tprime,phi);
      */

      for(int ich=0;ich<2;ich++) {
        GKPI0::SetReactionPar(1+ich);

        GKPI0::amplitude myAmp = myCache[ich].getAmpsForKine(qsq,xbj,t);

        /*
        printf("Amplituden: \n");
        std::cout << "Mmmp0 Mmpp0 Mpmp0 Mppp0" << std::endl;
        std::cout << myAmp.Mmmp0 << " " << myAmp.Mmpp0 << " " << myAmp.Mpmp0 << " " << myAmp.Mppp0 << std::endl;
        std::cout << "Twist-2 contribution! \n M0pp M0mp \n " << myAmp.M0pp << " " << myAmp.M0mp << "\n";
        */

        //and in the end, calculate the cross section
        double sigma = GKPI0::getCX(myAmp,w2);
        double sigmaTT = GKPI0::getCXTT(myAmp,w2,phi);
        double sigmaL = GKPI0::getCXL(myAmp,w2);
        double sigmaLT = GKPI0::getCXLT(myAmp,w2,phi);

        //printf("Final: %.4f %.4f %.4f dSigmaTot/dt = %.4f + e * %.4f + e* %.4f + sqrt(2e(1+e)) * %.4f\n",qsq,xbj,t,sigma,sigmaL,sigmaTT,sigmaLT);
        printf("%.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",qsq,xbj,t,sigma,sigmaL,sigmaTT,sigmaLT);
      }
    }
  }

  return 0;
}
