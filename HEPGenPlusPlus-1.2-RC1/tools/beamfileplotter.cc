
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <cassert>

//for making histograms the hepgen-way
#include "hbeamfile.h"


#include "TH2D.h"
#include "TH1D.h"
#include "TFile.h"

using namespace std;


int main (int argc, char** argv){
    printf("------- Beamfile-Analyzer \n");
    printf("------- A simple test-implementation of the beamfile\n");
    if (argc < 3)
    {
        printf("------- Too few arguments, use %s [INPUTFILE] [OUTPUT_HISTO_FILE] -[OPTIONS] \n",argv[0]);
        printf("------- Histograms will be created with ascii and/or root format if enabled in cmake \n");
        printf("Options:\n ");
        printf("\t -v \t Verbose mode - prints each event \n");
        return -1;
    }
    
    int beam,halo;
    beam = halo = 0;
    
    HBeamFile* myFile = new HBeamFile(argv[1]);
//     myFile->loadFile(argv[1]);
    if (!myFile->good()){
      printf("Could not open beamfile %s !\n",argv[1]);
      exit(-1);
    }
    
    TH2D* spaceHist = new TH2D("x vs y position","x vs y position",500,-50,50,500,-50,50);
    TH1D* eHist = new TH1D("energy","energy",200,0,200);
    
    
    for (int a =0; a < myFile->getNEntries();a++){
      if (a % 10000 == 0)
	printf("Event: %i\n",a);
      HBeamEntry myEntry = myFile->getNextEntry(Both);
      spaceHist->Fill(myEntry.position.X(),myEntry.position.Y());
      eHist->Fill(myEntry.energy);
      if (myEntry.type == Beam)
	beam++;
      else if (myEntry.type == Halo)
	halo++;
    }
    
    TFile* myOutFile = new TFile(argv[2],"RECREATE");
    
    if (!myOutFile->IsOpen()){
      printf("Could not open outfile %s\n",argv[2]);
      exit(-1);
    }
    spaceHist->Write();
    eHist->Write();
    myOutFile->Write();
    myOutFile->Close();
    
    printf("Beam entries: %i \n Halo entries: %i \n",beam,halo);
    
    
}