#include <iostream>
#include <fstream>
#include "hlorentzvector.h"
#include "hVGGgen.h"
#include "hMosseGen.h"
#include "hWeightingInterface.h"
#include "reweightKine.h"

int main(int argc, char** argv) {

    HLorentzVector muIn;
    muIn.setVector(HVector3(0,0,160));
    muIn.setEnergy(160);

    //this kinematic gives a VERY NARROW distribution in Herves code
    double Qsq= 4.0;
    double xbj = 0.03;
    double phi =0.0;

    double t = -0.001;
    double tmax = -2.5;
    int  nStep = 250;
    double tStep = (tmax - t)/nStep;

    hWeightInterface myWeightInterface;
    myWeightInterface.qsq = Qsq;
    myWeightInterface.xbj = xbj;
    myWeightInterface.phir = phi;
    myWeightInterface.t = t;
    myWeightInterface.beamE = 160;
    myWeightInterface.clept = +1;
    myWeightInterface.slept = -1;

    char fileName[200];
    for (int qsqBin = 0; qsqBin<2; qsqBin++) {
      myWeightInterface.phir = 0;
        for (int phiBin = 0; phiBin < 12; phiBin++) {
	 myWeightInterface.t = t;
            sprintf(fileName,"Q2_%f_phi_%f.dat",myWeightInterface.qsq,myWeightInterface.phir);
            ofstream outFile;
	    outFile.open(fileName);
	    outFile <<"Q2: " << myWeightInterface.qsq << " phi: " << myWeightInterface.phir << endl;
	    for (int i =0; i < nStep; i++) {
                double weightMoutarde[3], weightMosse[3];
                myWeightInterface.t += tStep;
                HPhysicsGenDVCSVGG::calcWeights(myWeightInterface,weightMoutarde[0],weightMoutarde[1],weightMoutarde[2]);
                HPhysicsGenDVCSMosse::calcWeights(myWeightInterface,weightMosse[0],weightMosse[1],weightMosse[2]);
		outFile << myWeightInterface.t << " " << weightMosse[0] << " " << weightMoutarde[0] << endl;
            }
            myWeightInterface.phir += 30.*(M_PI/180.);
        }
        myWeightInterface.qsq *= 2.;
    }

    cout << "FIN" << endl;

    return 0;
}

