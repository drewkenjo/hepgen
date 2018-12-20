#ifndef H_WEIGHT_INTERFACE
#define H_WEIGHT_INTERFACE

#include "hlorentzvector.h"

typedef struct {
    double xbj;
    double qsq;
    double t;
    double y;
    double s;
    double nu;
    double beamE;
    double b0;
    double xbj0;
    double alphap;
    HLorentzVector MuIn;
    HLorentzVector MuOut;
    HLorentzVector gammaOut;
    double phir;
    double clept;
    double slept;

    double tprim;

    //kinematic range used for original generation -- important for phase-space-calculation

    double numin;
    double numax;
    double qsqmin;
    double qsqmax;

    double slpin;
    double slpinn;
    double tmin;
    double tmax;


    double w2;

} hWeightInterface;


#endif
