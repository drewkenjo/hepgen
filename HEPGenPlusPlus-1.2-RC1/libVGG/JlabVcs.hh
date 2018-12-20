#include <fstream>
#include <string.h>
#include <math.h>

#include "myTHEO.hh"
#include "GPDQ.hh"

const Complexe me(0.511e-3);
const Complexe mmu(0.105);
const Complexe i(0.,1.);
const Complexe ni(0.,-1.);
const Complexe PI(3.14159265359);
const Complexe ELEC(0.3028619);

const CMatrix g0("gamma0",0);
const CMatrix g1("gamma1",1);
const CMatrix g2("gamma2",2);
const CMatrix g3("gamma3",3);
const CMatrix gamma5("gamma5",5);

const CM_Lorentz_Tensor sigma("sig",0,"rho",0,"sigma tensor");
const CM_Lorentz_Vector gamm("mug",1,g0,g1,g2,g3);


void SigmaPM(double* SP,double* SM,
             double ek, double eQ2, double exb,
             double et,double ephi,
             int phasespace,
             GpdInfo info);




// SP: positive helicty incident lepton
// SM: negative helicty incident lepton

// ek = incoming lepton energy in GeV/c
// eQ2 in GeV^2
// et = transfert in Gev^2
// ephi =  azimuthal angle in GRAD.

// phasespace: lab(0) or invariant(1) cross section?
// type: factorised(0) or xi dependant(1) SPDs ? (MRST98)"<<endl;
// bv: bv profile parameter for valence DD"<<endl;
// bs: bs profile parameter for sea DD"<<endl;
// dterm: D-term ? (yes=1)"<<endl;

void cinematiqueOZT(C_Lorentz_Vector*,C_Lorentz_Vector*,C_Lorentz_Vector*,C_Lorentz_Vector*,C_Lorentz_Vector*,C_Lorentz_Vector*,C_Lorentz_Vector*,Complexe ,Complexe ,Complexe ,Complexe ,Complexe ,Complexe );



