
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include "lfread.h"
#include "lfwrite.h"
#include <sstream>
#include <cassert>


using namespace std;


int main(int argc, char** argv)
{
    printf("------- LEPTO-SPLITTER \n");
    printf("------- A simple test-implementation of the lepto-library lfrw\n");
    if (argc < 4)
    {
        printf("------- Too few arguments, use %s [INPUTFILE] [OUTPUT_FILE_PREFIX] [NUMBER_OF_EVENTS] \n",argv[0]);
        return -1;
    }

    lfr::lEvent myEvent;
    lfr::lHeader myHeader;
    fstream inputFile;
    fstream outputFile;
    string inputFileName = argv[1];
    string outputPrefix = argv[2];
    int eventCount = atoi(argv[3]);
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

    //----------------Loop over all the events in the file---------
    while (lfr::lNextEvent(inputFile,myEvent,largeHeaders)!= -1 && !inputFile.eof())
    {
        if (eventCounter % eventCount == 0) {
            if (outputFile.is_open()) {
                outputFile.flush();
                outputFile.close();
            }
            fileCounter++;
            stringstream fileName;
            printf("---------- Starting new File: #%i at event count: %i\n",fileCounter,eventCounter);
            fileName << outputPrefix << "_" << fileCounter;
            returnWriter = lfr::lOpenFile(fileName.str(),outputFile);
            assert(returnWriter == 0);
            lfr::lWriteHeader(outputFile,myHeader);
        }


        lfr::lWriteNextEvent(outputFile,&myEvent);
        outputFile.flush();
        eventCounter++;


        lfr::lFreeBuffer(myEvent);
    }


    printf("------- DONE! All events from the file have been read! \n");

}
