
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>
#include <cassert>

//reading / writing leptofiles
#include "lfread.h"
#include "lfwrite.h"

//for making histograms the hepgen-way
#include "hbooker.h"
#include "hconstants.h"
using namespace std;


int main(int argc, char** argv)
{
    printf("------- LEPTO-Analyzer \n");
    printf("------- A simple test-implementation of the lepto-read-library lfr\n");
    if (argc < 3)
    {
        printf("------- Too few arguments, use %s [INPUTFILE] [OUTPUT_HISTO_FILE] -[OPTIONS] \n",argv[0]);
        printf("------- Histograms will be created with ascii and/or root format if enabled in cmake \n");
        printf("Options:\n ");
        printf("\t -v \t Verbose mode - prints each event \n");
        return -1;
    }


    //--------------- Check if verboseMode shall be enabled
    bool verbose = false;
    if (argc > 3)
        verbose = (strcmp(argv[3],"-v")==0);

    lfr::lEvent myEvent;
    lfr::lHeader myHeader;
    fstream inputFile;
    string inputFileName = argv[1];
    int fileCounter = 0;
    int eventCounter = 0;

    //---------------- Print the library logo -----
    lfr::lPrintLogo();


    //----------------Make the file ready---------
    int returnCode = lfr::lLoadFile(inputFileName,inputFile);
    if (returnCode == -1) {
        printf("------- Error loading the file!\n");
        return -1;
    }
    else if(returnCode == 0) {
        printf("------- Found good LEPTOFILE with GFORTRAN headers \n");
    }
    else if (returnCode == 1) {
        printf("------- Found good LEPTOFILE with F77 headers \n");
    }



    bool largeHeaders = (returnCode == 0);

    //----------------Read the header of the file---------
    if (lfr::lReadHeader(inputFile,myHeader,largeHeaders) == -1)
    {
        printf("------- Error while reading header of the LEPTOFILE \n");
        return -1;
    }

    printf("------- Header has been read successfully!! -------- \n");

    int returnWriter;
    float e,x,y;

    //----------------Initialize the Histogram manager
    HBooker* myHistos = new HBooker("");
    //----------------Always add the ASCII backend
    myHistos->addASCII();
    bool useROOT=false;
    float dummyWeight = 1.0;

    //----------------Use ROOT if possible
#ifdef USE_ROOT
    useROOT = true;
    myHistos->addROOT();
#endif

    //----------------Add 2 simple histograms
    myHistos->addHBook1F(&myEvent.q2,&dummyWeight,100,0,20,"Q^{2} unweighted");
    myHistos->addHBook1F(&myEvent.q2,&myEvent.uservar[2],100,0,20,"Q^{2} weighted (sum)");



    myHistos->addHBook1F(&e,&dummyWeight,160,140.,180.,"Energy of beam particle");
    myHistos->addHBook1F(&x,&dummyWeight,200,-100.,100.,"Vertex position x");
    myHistos->addHBook1F(&y,&dummyWeight,200,-100.,100.,"Vertex position y");


    //----------------2 Dimensional also works
    myHistos->addHBook2F(&myEvent.q2,&myEvent.nu,&myEvent.uservar[9],100,100,0.5,5,80,155,"Q^{2} vs #{nu} plot");
    
    myHistos->addHBook2F(&x,&y,&dummyWeight,200,200,-100,100,-100,100,"x vs y");



    //----------------Loop over all the events in the file---------
    while (lfr::lNextEvent(inputFile,myEvent,largeHeaders)!= -1 && !inputFile.eof() )
    {

        e = myEvent.beamParts[0].p[3];
        x = myEvent.uservar[0];
        y = myEvent.uservar[1];


        if(myEvent.uservar[2] < 0.0 || myEvent.uservar[2] != myEvent.uservar[2]) {
            printf("Nan or neg weight found: %e \n",myEvent.uservar[2]);
            continue;
        }
        eventCounter++;
        if (eventCounter % 10000 == 0)
            printf("Event: %i\n",eventCounter);

        if (verbose) {
            printf("Event read!! \n");
            printf("------------------------------------------------------------------------\n");
            printf("----------Kinematics\nQ^2:\t%f\nXbj:\t%f\nNu:\t%f\n",myEvent.q2,myEvent.x_bj,myEvent.nu);
            printf("----------BeamParticles in Lujets: \n");

            for (int i = 0; i < myEvent.nBeamParticle; i++) {
                printf("Beamparticle P %i: \t %f, %f, %f, %f, %f\n",i,myEvent.beamParts[i].p[0],myEvent.beamParts[i].p[1],myEvent.beamParts[i].p[2],myEvent.beamParts[i].p[3],myEvent.beamParts[i].p[4]);
                double mass =sqrt(fabs(pow(myEvent.beamParts[i].p[3],2.0) - (pow(myEvent.beamParts[i].p[0],2.0)+pow(myEvent.beamParts[i].p[1],2.0)+(myEvent.beamParts[i].p[2],2.0))));
                printf("mass %0.2e\n",mass);
            }


            printf("----------Uservars: ");
            for (unsigned int i = 0; i < 20; i++)
                cout <<i<< ": " << myEvent.uservar[i]<< " ";
            printf("\n");


            printf("------------------------------------------------------------------------\n");
            //clean up the buffer for the event - not doing so will use new without delete - which is bad!
        }

        //---- fill the histograms
        myHistos->fill();

        lfr::lFreeBuffer(myEvent);
    }

    myHistos->dumpToFile(argv[2]);
    printf("------- DONE! All events from the file have been read! \n");

}
