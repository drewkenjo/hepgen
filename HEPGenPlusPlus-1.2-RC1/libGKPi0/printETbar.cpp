#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fstream>

#include <string>
#include <sys/stat.h>
#include <unistd.h>

#include "libGKPi0.h"
#include "gkSubProcessTable.h"
#include "TString.h"

int main (int argc, char** argv) {
//  GKPI0::SetETbarUCoefs(0.430284, 0.184742, 9.77666);
//  GKPI0::SetETbarDCoefs(1.83919, -38.1062, 339.345, -832.334, 668.352);

/*
  GKPI0::SetETbarUCoefs(0.321477, 1.80601, -0.556723);
  GKPI0::SetETbarDCoefs(18.1052, -571.803, 5167.55, -11884.9, 9527.68);
*/

/*
  GKPI0::SetETbarUCoefs(-1.14386, 19.5226, -17.3655);
  GKPI0::SetETbarDCoefs(-1.77821, 25.4954, -199.433, 459.157, -390.572);
  GKPI0::SetETbarUNorm(2.01537);
  GKPI0::SetETbarDNorm(-4.77236);
*/

  GKPI0::SetETbarUNorm(3.73173);
  GKPI0::SetETbarDNorm(6.69741);

  GKPI0::SetETbarUtSlope(-1.69245);
  GKPI0::SetETbarDtSlope(0.691707);

//  GKPI0::SetAlphaStr(0.933758);

  for(double x=-0.975;x<1;x+=0.05) {
    for(double xi=0.025;xi<1;xi+=0.05) {
      double etbarU = GKPI0::EBarU(x,xi,-0.1,2,0.5);
      double etbarD = GKPI0::EBarD(x,xi,-0.1,2,0.5);
      printf("%.4f %.4f %.4f %.4f\n",x,xi,etbarU,etbarD);
    }
  }

  return 0;
}
