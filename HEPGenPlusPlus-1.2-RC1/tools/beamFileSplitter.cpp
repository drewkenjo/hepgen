
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <cassert>


#include "hbeamfile.h"


using namespace std;


int main(int argc, char** argv)
{
   printf("------- BEAMFILE-SPLITTER \n");
    if (argc < 5)
    {
        printf("------- Too few arguments, use %s [INPUTFILE] [OUTPUT_FILE_PREFIX] [NUMBER_OF_EVENTS] [BEAM_TYPE_WANTED]\n",argv[0]);
        return -1;
    }
    
    //open input file
    HBeamFile inputFile (argv[1]);
    HBeamType wantedType;
    
    
    //open output file
    HBeamFile outPutFile;
    
    int wantedEventsPerFile = atoi(argv[3]);
    int beamTypeWanted = atoi(argv[4]);
    switch (beamTypeWanted){
      case 1: wantedType = Beam;break;
      case 2: wantedType = Halo;break;
      case 3: wantedType = Both;break;
      default: wantedType = Both; 
    }
    
    int fileCounter = 0;
    
    char buffer[200];
    for (int i =0; i < inputFile.getNEntries();i++){
      if (i%wantedEventsPerFile == wantedEventsPerFile-1){
	sprintf(buffer,"%s_%i.bin",argv[2],fileCounter);
	outPutFile.saveFile(buffer,false);
	outPutFile.clearList();
	fileCounter++;
      }
      outPutFile.addEntry(inputFile.getNextEntry(wantedType));      
    }    
    return 0;

}