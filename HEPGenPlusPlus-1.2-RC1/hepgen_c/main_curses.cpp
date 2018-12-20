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


#include <ncurses.h>
#include <cstdlib>


using namespace std;

void printstringlist ( vector<string> _in ) {
    cout << "Printing Vector of String:" << endl;
    cout << "Entrynum: Value" << endl;

    for ( unsigned int i = 0; i < _in.size(); i++ ) {
        cout << i << ": " << _in.at ( i ) << endl;
    }
}


//functions
void mvwprintwcenter ( WINDOW *win, int _y, string _string );
void cleanUp() {
    endwin();
}

//globals
WINDOW *titleWin, *infoContainer;
int x, y;
int evCountMax;
double ratioGlob;
int nEvent;
char generatorName[40];
char outFileFormat[10];
HParamManager* myParamManGlobal;
string outFile;
string outFileHistos;

void printLogo() {
    wbkgd ( titleWin, 1 );
    wcolor_set ( titleWin, 1, 0 );

    for ( int i = 0; i < hepconst::logolength; i++ ) {
        mvwprintwcenter ( titleWin, i + 1, hepconst::logo[i] );
    }
    refresh();
    wrefresh ( titleWin );
}


const string currentDateTime() {
    time_t     now = time ( 0 );
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime ( &now );
    strftime ( buf, sizeof ( buf ), "%Y-%m-%d.%X", &tstruct );
    return buf;
}
void updateInfos() {
    int line = 1;
    mvwprintw ( infoContainer, line++, 1, "HEPGen++ version %s running", hepconst::version.c_str() );
    mvwprintw ( infoContainer, line++, 1, "---------------------------" );

    switch ( myParamManGlobal->getIVECM().at ( 0 ) ) {
    case 0:
        sprintf ( generatorName, "DVCS - Andrzej's Model" );
        break;
    case 1:
        sprintf ( generatorName, "Pi0" );
        break;
    case 2:
        sprintf ( generatorName, "Rho0" );
        break;
    case 3:
        sprintf ( generatorName, "Phi" );
        break;
    case 4:
        sprintf ( generatorName, "Omega" );
        break;
    case 5:
        sprintf ( generatorName, "Rho+" );
        break;
    case 13:
        sprintf ( generatorName, "Pam Bethe Heitler" );
        break;
    case 11:
        sprintf ( generatorName, "H. Moutarde's GK11 model" );
        break;
    case 12:
        sprintf ( generatorName, "L. Mosse's VGG code" );
        break;
    default:
        sprintf ( generatorName, "INVALID ID FOUND %i", myParamManGlobal->getIVECM().at ( 1 ) );
    };

    mvwprintw ( infoContainer, line++, 1, "Generator: %s", generatorName );
    mvwprintw ( infoContainer, line++, 1, "Nu limits: %.2e - %.2e", myParamManGlobal->getNuMin(), myParamManGlobal->getNuMax() );
    mvwprintw ( infoContainer, line++, 1, "Q^2 limits: %.2e - %.2e", myParamManGlobal->getQsqMin(), myParamManGlobal->getQsqMax() );
    mvwprintw ( infoContainer, line++, 1, "tprime limits: %.2e - %.2e", myParamManGlobal->gettMin(), myParamManGlobal->gettMax() );
    mvwprintw ( infoContainer, line++, 1, "DD is %i", myParamManGlobal->getStruct()->LST.at ( 19 ) );
    mvwprintw ( infoContainer, line++, 1, "---------------------------" );
    mvwprintw ( infoContainer, line++, 1, "Outfile: %s\n", outFile.c_str() );

    if ( myParamManGlobal->getKeyContents ( "ENABLE_GFORTRAN" ).at ( 1 ) == "1" )
        sprintf ( outFileFormat, "gfortran" );
    else
        sprintf ( outFileFormat, "f77" );
    mvwprintw ( infoContainer, line++, 1, "Output is %s header compatible", outFileFormat );
//   mvwprintw(infoContainer,line++,1,"Debug is enabled",outFileFormat);


}


const void updateStatus ( int _nEvent ) {
    //calcpercent
    _nEvent++;
    nEvent = _nEvent;
    if (_nEvent > evCountMax)
        _nEvent = evCountMax;
    double ratio = ( double ) _nEvent / evCountMax;
    ratio *= 100.;
    int evCountPercent = (double)evCountMax / 100.;
    if ( _nEvent % evCountPercent == 0 ) {
        mvwprintw ( infoContainer, 20, 1, "Progress: %i of %i -- %.2f %\n", _nEvent, evCountMax, ratio );
        refresh();
        wrefresh ( infoContainer );
    }
}



void printhelp() {
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
    cout << "Written 2013-2016 by Christopher Regali (christopher.regali@cern.ch)" << endl;
    cout << "Original HEPGen-Page: http://project-gpd-full-chain-mc.web.cern.ch/project-gpd-full-chain-mc/hepgen/" << endl;
}



int main ( int _argc, char** _argv ) {

    HCardParser parserTester;
    if ( _argc < 3 ) {
        if ( _argc >= 2 )
            if ( strcmp ( _argv[1], "-h" ) == 0 ) {
                printhelp();
                return 0;
            }

        printf ( "usage: %s [path/to/hepgen.data] [random-seed]\n", _argv[0] );
        printf ( "or do %s -h for full help\n", _argv[0] );
        return 1;
    }


    string fileName = _argv[1];
    string seed = _argv[2];


    HGenManager* tmp = HGenManager::getInstance();
    tmp->loadSettings ( fileName );
    tmp->setupGenerator();
    HParamManager ParamManager ( fileName );
    myParamManGlobal = &ParamManager;
    tmp->setSeed ( hephelp::StrToInt ( seed ) );

    if (ParamManager.getKeyContents("HISTOS_ASCII").at(1) == "1")
        tmp->enableASCIIBook();

#ifdef USE_ROOT
    if (ParamManager.getKeyContents("HISTOS_ROOT").at(1) == "1")
      tmp->enableROOTBook();
#endif

    bool gfortran = false;
    if ( ParamManager.getKeyContents ( "ENABLE_GFORTRAN" ).at ( 1 ) == "1" )
        gfortran = true;

    if ( ParamManager.getKeyContents ( "ENABLE_DEBUG" ).at ( 1 ) == "1" ) {

        tmp->enableDebug();
    }
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
    outFile = outfile;


    noecho();
    initscr();
    clear();
    atexit ( cleanUp );
    start_color();

    cbreak();
    keypad ( stdscr, TRUE );
    getmaxyx ( stdscr, y, x );

    init_pair ( 1, COLOR_RED, COLOR_BLACK );
    init_pair ( 2, COLOR_WHITE, COLOR_BLUE );


    evCountMax = ParamManager.getNumEvents();


    titleWin = newwin ( 10, x - 2, 1, 1 );
    box ( titleWin, 0, 0 );

    infoContainer = newwin ( y - 12, x - 2, 11, 1 );

    wbkgd ( infoContainer, COLOR_PAIR ( 2 ) );
    wcolor_set ( infoContainer, 2, 0 );

    box ( infoContainer, 0, 0 );
    //print header
    printLogo();
    mvwprintwcenter ( infoContainer, 0, "HEPGen++ - High excl. event generator" );
    updateInfos();


    refresh();
    wrefresh ( titleWin );
    wrefresh ( infoContainer );

    int infoX, infoY;

    tmp->enableOutPut ( outfile );

    double weightSum = tmp->startRun ( &updateStatus );

    flash();
    char buffer[70];
    sprintf(buffer,"GENERATION FINISHED!!!! Sum of weights is: %e\n",weightSum);

    mvwprintwcenter(infoContainer,22,buffer);
    mvwprintwcenter(infoContainer,23,"Write this down, then press ESC to exit!");

//   refresh();
    wrefresh ( infoContainer );
    int kx=0;
    while (kx != 27)
        kx = getch();



    /*
      refresh();
      wrefresh ( titleWin );
      wrefresh ( infoContainer );*/
    clrtoeol();
    cleanUp();

    return EXIT_SUCCESS;
}


void mvwprintwcenter ( WINDOW *win, int _y, string _string ) {
    int mx, my;
    getmaxyx ( win, my, mx );

    mvwprintw ( win, _y, ( mx - _string.length() ) / 2, _string.c_str() );
//     cout << "y " << _y << " strlen " << _string.length() << " str " << _string.c_str() << endl;


}
