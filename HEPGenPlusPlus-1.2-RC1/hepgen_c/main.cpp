/*
 * main.cpp
 *
 *  Created on: Jan 16, 2013
 *  Author: Christopher Regali
 *  christopher.regali@cern.ch
 *  Copyright (c) 2013 All Right Reserved
 */


#include "hvector.h"
#include "config.h"
#include "hlorentzvector.h"
#include "hcardparser.h"
#include "hparammanager.h"
#include "hpionicdata.h"
#include "hbooker.h"
#include "hbookbackendASCII.h"
#ifdef USE_ROOT
#include "hbookbackendROOT.h"
#endif
#include "hgenmanager.h"
#include <iostream>
#include <string>
#include <cstring>
#include <string.h>
#include <time.h>
using namespace std;

void printstringlist ( vector<string> _in ) {
    cout << "Printing Vector of String:" << endl;
    cout << "Entrynum: Value" << endl;

    for ( unsigned int i = 0; i < _in.size(); i++ ) {
        cout << i << ": " << _in.at ( i ) << endl;
    }
}

void printlogo() {
    cout << endl <<  endl;
    for ( int i = 0; i < hepconst::logolength; i++ )
        cout << hepconst::logo[i] << endl;
    cout << endl <<  endl;

}

const string currentDateTime() {
    time_t     now = time ( 0 );
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime ( &now );
    strftime ( buf, sizeof ( buf ), "%Y-%m-%d.%X", &tstruct );
    return buf;
}




void printhelp() {
    //printlogo();
    cout << "----------------------------------------------------------------" << endl;
    cout << "usage: HEPGen++ [path/to/hepgen.data] [random-seed] " << endl;
    cout << "----------------------------------------------------------------" << endl << endl;

    cout << "This is HEPGen++ version: " << hepconst::version << endl;
    cout << "Supported generators in this version: " << hepconst::generators << endl << endl << endl;
    cout << "------------------------!!!!!!!!!!!!!!!!------------------------" << endl;
    cout << "WARNING: DO NOT USE MULTILINED USERCARDS!!! - Common error is CUTL to look like: " << endl;
    cout << "'CUTL   0.0001  1.0  0.0  1.0  0.5  80.0  0.00" << endl << "  1000.0  5.0  155.0  0.0   200.0  0.00  6.28318'" << endl;
    cout << "Instead make it look like 'CUTL   0.0001  1.0  0.0  1.0  0.5  80.0  0.00  1000.0  5.0  155.0  0.0   200.0  0.00  6.28318' " << endl;
    cout << "----------------------------------------------------------------" << endl;
    cout << endl << endl << endl;
    cout << "HEPGen++ and libhepgen are a C++ rewrite of the original HEPGen by A. Sandacz." << endl;
    cout << "Written 2013-2014 by Christopher Regali (christopher.regali@cern.ch)" << endl;
    cout << "Original HEPGen-Page: http://project-gpd-full-chain-mc.web.cern.ch/project-gpd-full-chain-mc/hepgen/" << endl;
}

int evCountMax;

static char const spin_chars[] = "//-\\|";

const void printStatus ( int nEvent ) {
    if ( ( nEvent + 1 ) * 100 / evCountMax > nEvent * 100 / evCountMax ) {
        putchar ( spin_chars[ ( nEvent + 1 ) * 100 / evCountMax % sizeof ( spin_chars )] );
        putchar ( ' ' );
        fflush ( stdout );
    }

    if ( ( nEvent + 1 ) * 10 / evCountMax > nEvent * 10 / evCountMax )
        cout << ( nEvent + 1 ) * 100 / evCountMax << "% " << flush;
    putchar ( '\b' );
}



int main ( int _argc, char** _argv ) {
    printlogo();

    HCardParser parserTester;
    if ( _argc < 3 ) {
        if ( _argc >= 2 )
            if ( strcmp ( _argv[1], "-h" ) == 0 ) {
                printhelp();
                return 0;
            }
        cout << "usage: HEPGen++ [path/to/hepgen.data] [random-seed]" << endl;
        cout << "or do HEPGen++ -h for full help" << endl;
        cout << endl << endl << endl;
        return 1;
    }


    string fileName = _argv[1];
    string seed = _argv[2];

    HGenManager* tmp = HGenManager::getInstance();
    tmp->loadSettings ( fileName );
    tmp->setupGenerator();
    HParamManager ParamManager ( fileName );







    tmp->setSeed ( hephelp::StrToInt ( seed ) );



    if (ParamManager.getKeyContents("HISTOS_ASCII").at(1) == "1")
        tmp->enableASCIIBook();

#ifdef USE_ROOT
    if (ParamManager.getKeyContents("HISTOS_ROOT").at(1) == "1")
	tmp->enableROOTBook();
#endif



    bool gfortran = false;


    //TODO: Merge this into HGenManager
    if (ParamManager.getKeyContents("ENABLE_GFORTRAN").size() > 1)
        if ( ParamManager.getKeyContents ( "ENABLE_GFORTRAN" ).at ( 1 ) == "1" )
            gfortran = true;



    if (ParamManager.getKeyContents("ENABLE_DEBUG").size() > 1)
        if ( ParamManager.getKeyContents ( "ENABLE_DEBUG" ).at ( 1 ) == "1" ) {
            cout << " Debugging enabled!" << endl;
            tmp->enableDebug();
        }


    evCountMax = ParamManager.getNumEvents();


    tmp->enableLeptoOutput ( gfortran );

    string outfile;
    if ( ParamManager.getOutName() == "UNSET" ) {
        outfile = "HEPGen++-outfile-";
        outfile.append ( currentDateTime() );
        outfile.append ( "-" );
        outfile.append ( seed );
    }
    else {
        outfile = ParamManager.getOutName();
        outfile.append ( "_bin" );
    }

    tmp->enableOutPut ( outfile );

    double weightSum =  tmp->startRun ( &printStatus );
    printf("Sum of weights is: %e\n",weightSum);

    /* this is for debug reasons, if we ever need to check with the oneshot-routine thats used in TGEANT*/
//     double beamE = 165.47143;
//     for (int i =0; i < 100000; i ++)
//     {
//       cout << "----ELEPT: " << beamE << endl;
//       ParamManager.getStruct()->ELEPT = beamE;
//       ParamManager.getStruct()->PARL.at(2) = beamE;
//
//       tmp.oneShot();
//
//       tmp.getEvent()->printDebug();
//     }

    return EXIT_SUCCESS;
}


