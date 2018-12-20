/*!
 *  \file hconstants.h
 *  \date Created on: Feb 7, 2013
 *  \author: Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2013 All Right Reserved
 */






#include <string>


#ifndef HCONST_H_
#define HCONST_H_

#define HEPGEN_INTERFACE_NEW
#include "config.h"
/*! \namespace hepconst
 *  \brief This namespace contains physics constants, masses, particle IDs and further information 
 */
namespace hepconst
{
//version in major.minor.revision
//const string version="0.0.40";

//dd-Parameters
const double ddmassMaxSq = 30.0;
const double ddmassMinSq = 1.165;
const double ddmassPeak = 2.0;
const double ddmassThreshold = 1.165;
const double sigel = 7.0;
const double par1sd = 0.68;
const double par2sd = 36.0;


//squared masses
const double w2mu = 0.011164; /*! Mass-Square of Myon */
const double w2pi = 0.018219; /*! Mass-Square of Pion neutral */
const double mpi = 0.134977776; /*! Mass of Pion neutral */
const double w2pic = 0.019488; /*! Mass-Square of Pion charged */
const double mpic = 0.139599427;  /*! Mass of Pion charged */
const double w2kaon = 0.24374; /*! Mass-Square of Kaon */

const double w2KaonCharged = 0.243716980329; /*! Mass-Square of Charged-Kaon */
const double w2proton = 0.88032; /*! Mass-Square of Proton */
const double w2electron = 0.000000261; /*! Mass-Square of Electron */
const double w2prma = 0.938256; /* Mass of Proton */
const double w2neutron = 0.882783; /* Mass Square of Neutron */
const double w2omega = 0.61254;

//constants for the phi-generator
const double mPhi = 1.019445;
const double gPhi = 0.00426;

//these constants are specifically for the rho0-generator
const double mRho0 = .770;

//constants for j/psi generators
const double mJPsi = 3.0969;

//constants for the omega generator
const double mOmega = 0.78265;
const double gOmega = 0.00849;
const double qOmega = 0.327;
const double probMax = 0.016208238;

//new omega width
const double qOmegaNew = 0.380;

// constants for the rho mass generation
// 
const double gRho0 = .153;
const double qCMS0 = .358;

//const double mPhi  = 1.01946;

const double randomValues[] = {0.945894897, 0.473478496, 0.951527894, 0.429719746, 0.0912738442, 0.310068548};

//assorted constants
const double alfa = 0.007299; /*! alfa */

//dd constants

const double ptpi = 0.330;
const double ptpr = 0.5;




//particle types
const int typeMuonPlus = -13;
const int typeMuonMinus = 13;
const int typeGamma = 22;
const int typeElectronMinus = 11;
const int typeElectronPlus = -11;
const int typeProton = 2212;
const int typeProtonExcited = 2210;
const int typeNeutron = 2112;
const int typeNeutronExcited = 2110;
const int typeJPsi = 443;


const int typeRhoPlus = 213;


const int typeRho0 = 113;
const int typePi0 = 111;
const int typePiPlus = 211;
const int typePiMinus = -211;

const int typePhi = 333;
const int typeKPlus = 321;
const int typeKMinus = -321;

const int typeOmega = 223;



enum eventType {DVCS,RHO0,PI0,PHI,RHOPLUS,OMEGA,JPSI,OMEGANEW};
enum nucType {Proton,Neutron};

#ifdef USE_ROOT

#ifdef USE_EXPERIMENTAL
const std::string version = "1.1-RC2-01-EXPERIMENTAL";
#else
const std::string version = "1.1-RC2-01";
#endif

#else

#ifdef USE_EXPERIMENTAL
const std::string version = "1.1-RC2-02-NO_ROOT-EXPERIMENTAL";
#else
const std::string version = "1.1-RC2-02-NO_ROOT";
#endif

#endif

#ifdef USE_EXPERIMENTAL
const std::string generators = "DVCS (IVECM=0), RHO0 (IVECM=2), PI0 (IVECM=1), PHI (IVECM=3), J/Psi->ee,mumu (IVECM=4), Omega (IVECM=6), Rho+(IVECM=7), Omega->Pi0Gamma (IVECM=8), DVCS-LMOSSE(IVECM=12), PAM-BH(IVECM=13)";
#else
const std::string generators = "DVCS (IVECM=0), RHO0 (IVECM=2), PI0 (IVECM=1), PHI (IVECM=3), J/Psi->ee,mumu (IVECM=4), Omega (IVECM=6), Rho+(IVECM=7), Omega->Pi0Gamma (IVECM=8), DVCS-LMOSSE(IVECM=12), PAM-BH(IVECM=13)";
#endif

const int logolength = 8;
const std::string logo[] = {
"  _    _ ______ _____   _____                        ",   
" | |  | |  ____|  __ \\ / ____|            _     _   ",
" | |__| | |__  | |__) | |  __  ___ _ __ _| |_ _| |_ ",
" |  __  |  __| |  ___/| | |_ |/ _ \\ '_ \\_   _|_   _|",
" | |  | | |____| |    | |__| |  __/ | | ||_|   |_|  ",
" |_|  |_|______|_|     \\_____|\\___|_| |_|           ",
"                                                    ",
"                                                    "
};






/*
const std::string logo[] = {
" ██░ ██ ▓█████  ██▓███    ▄████ ▓█████  ███▄    █ ",
"▓██░ ██▒▓█   ▀ ▓██░  ██▒ ██▒ ▀█▒▓█   ▀  ██ ▀█   █ ",
"▒██▀▀██░▒███   ▓██░ ██▓▒▒██░▄▄▄░▒███   ▓██  ▀█ ██▒",
"░▓█ ░██ ▒▓█  ▄ ▒██▄█▓▒ ▒░▓█  ██▓▒▓█  ▄ ▓██▒  ▐▌██▒",
"░▓█▒░██▓░▒████▒▒██▒ ░  ░░▒▓███▀▒░▒████▒▒██░   ▓██░",
" ▒ ░░▒░▒░░ ▒░ ░▒▓▒░ ░  ░ ░▒   ▒ ░░ ▒░ ░░ ▒░   ▒ ▒ ",
" ▒ ░▒░ ░ ░ ░  ░░▒ ░       ░   ░  ░ ░  ░░ ░░   ░ ▒░",
" ░  ░░ ░   ░   ░░       ░ ░   ░    ░      ░   ░ ░ ",
" ░  ░  ░   ░  ░               ░    ░  ░         ░ "
};*/




}








#endif
