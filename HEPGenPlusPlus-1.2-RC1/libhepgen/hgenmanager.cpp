#include "hgenmanager.h"


HGenManager* HGenManager::instance = NULL;


void HGenManager::setParamManager(HParamManager* _params)
{
    paramMan = _params;
}

void HGenManager::setupGeneratorOneShot()
{
    //we have to do this for backward compat of after cpp-2011.. one could not initialize non-consts there.
    enabledASCIIBook = false;
    enabledROOTBook = false;
    enabledBeamFile = false;


    debug = false;

    if (paramMan == NULL)
    {
        cout << "Error: Called setupGenerator before loading a datacard or providing a parammanager!" << endl;
        return;
    }
    hbookMan = new HBooker("all");
    ddBookMan = new HBooker("dd");

    randMan = new CLHEP::HepRandom;
    CLHEP::HepRandomEngine* tmp = new CLHEP::RanluxEngine;
    tmp->setSeed(1,1);
    randMan->setTheEngine(tmp);
    randEng = tmp;
    bookingEvent = new HEvent(paramMan);

    myEvent = new HEvent(paramMan);

    int physics = paramMan->getIVECM().at(0);

    printf("physics: %i \n",physics);
    switch (physics)
    {
    case 0:
        generator = new HPhysicsGenDVCS(randMan, paramMan, myEvent, hbookMan, ddBookMan);
        break; //DVCS generator chosen
    case 1:
        generator = new HPhysicsGenPI(randMan, paramMan, myEvent, hbookMan, ddBookMan);
        break; // Pi0 generator chosen
    case 2:
        generator = new HPhysicsGenRHO(randMan, paramMan, myEvent, hbookMan, ddBookMan);
        break; //Rho0 generator chosen
    case 3:
        generator = new HPhysicsGenPHI(randMan, paramMan, myEvent, hbookMan, ddBookMan);
        break; //Phi generator chosen
    case 4:
        generator = new HPhysicsGenJPsi(randMan, paramMan, myEvent, hbookMan, ddBookMan);
        break; //j/psi generator chosen
    case 6:
        generator = new HPhysicsGenOMEGA(randMan, paramMan, myEvent, hbookMan, ddBookMan);
        break; //OMEGA generator chosen

    case 7:
        generator = new HPhysicsGenRHOPlus(randMan, paramMan, myEvent, hbookMan, ddBookMan);
        break; //RHOPlus generator chosen
    case 8:
        generator = new HPhysicsGenOMEGAPiGamma(randMan, paramMan, myEvent, hbookMan, ddBookMan);
        break; //Omega to pi0+gamma
        


#ifdef USE_EXPERIMENTAL
//     case 11:
//         generator = new HPhysicsGenDVCSVGG(randMan, paramMan, myEvent, hbookMan, ddBookMan);
//         break; //VGG-DVCS generator chosen
    case 12:
        generator = new HPhysicsGenDVCSMosse(randMan, paramMan, myEvent, hbookMan, ddBookMan);
        break; //VGG-DVCS generator chosen
#endif
    case 13:
        generator = new HPhysicsGenPAMBH(randMan, paramMan, myEvent, hbookMan, ddBookMan);
        break; //VGG-DVCS generator chosen

    }

    if (generator == NULL)
    {
        std::cout << "invalid generator chosen! aborting!" << std::endl;
        return;
    }
    //set the parameters needed for correct header-writing
    generator->setParameters();

}



void HGenManager::setupGenerator()
{
    //we have to do this for backward compat of after cpp-2011.. one could not initialize non-consts there.
    enabledASCIIBook = false;
    enabledROOTBook = false;
    enabledBeamFile = false;


    debug = false;

    if (paramMan == NULL)
    {
        cout << "Error: Called setupGenerator before loading a datacard or providing a parammanager!" << endl;
        return;
    }
    hbookMan = new HBooker("all");
    ddBookMan = new HBooker("dd");

    randMan = new CLHEP::HepRandom;
    CLHEP::HepRandomEngine* tmp = new CLHEP::RanluxEngine;
    tmp->setSeed(1,1);
    randMan->setTheEngine(tmp);
    randEng = tmp;
    bookingEvent = new HEvent(paramMan);

    myEvent = new HEvent(paramMan);

    int physics = paramMan->getIVECM().at(0);

    switch (physics)
    {
    case 0:
        generator = new HPhysicsGenDVCS(randMan, paramMan, myEvent, hbookMan, ddBookMan);
        break; //DVCS generator chosen
    case 1:
        generator = new HPhysicsGenPI(randMan, paramMan, myEvent, hbookMan, ddBookMan);
        break; // Pi0 generator chosen
    case 2:
        generator = new HPhysicsGenRHO(randMan, paramMan, myEvent, hbookMan, ddBookMan);
        break; //Rho0 generator chosen
    case 3:
        generator = new HPhysicsGenPHI(randMan, paramMan, myEvent, hbookMan, ddBookMan);
        break; //Phi generator chosen
    case 4:
        generator = new HPhysicsGenJPsi(randMan, paramMan, myEvent, hbookMan, ddBookMan);
        break; //j/psi generator chosen
    case 6:
        generator = new HPhysicsGenOMEGA(randMan, paramMan, myEvent, hbookMan, ddBookMan);
        break; //OMEGA generator chosen

    case 7:
        generator = new HPhysicsGenRHOPlus(randMan, paramMan, myEvent, hbookMan, ddBookMan);
        break; //RHOPlus generator chosen
    case 8:
        generator = new HPhysicsGenOMEGAPiGamma(randMan, paramMan, myEvent, hbookMan, ddBookMan);
        break; //Omega to pi0+gamma
        
        
    case 13:
        generator = new HPhysicsGenPAMBH(randMan, paramMan, myEvent, hbookMan, ddBookMan);
        break; //RHOPlus generator chosen

#ifdef USE_EXPERIMENTAL
 /*   case 11:
        generator = new HPhysicsGenDVCSVGG(randMan, paramMan, myEvent, hbookMan, ddBookMan);
        break; //VGG-DVCS generator chosen
 */   case 12:
        generator = new HPhysicsGenDVCSMosse(randMan, paramMan, myEvent, hbookMan, ddBookMan);
        break; //VGG-DVCS generator chosen
#endif

    }

    if (generator == NULL)
    {
        std::cout << "invalid generator chosen! aborting!" << std::endl;
        return;
    }
    //set the parameters needed for correct header-writing
    generator->setParameters();

    outputMan = new HOutPut(paramMan);

    string beamFileName = paramMan->getKeyContents("BEAMFILE").at(1);

    if (paramMan->getStruct()->I_BEAMREAD != 0)
    {
        HBeamType typeWanted;
        switch (paramMan->getStruct()->I_BEAMREAD)
        {
        case 1:
            typeWanted = Beam;
            break;
        case 2:
            typeWanted = Halo;
            break;
        case 3:
            typeWanted = Both;
            break;
        }
        enableBeamFile(beamFileName, typeWanted);
        cout << "Enabled beamfile" << endl;
    }
}



void HGenManager::loadSettings(string _fileName)
{
    paramMan = new HParamManager(_fileName);
}

HGenManager::HGenManager()
{
    outputMan = NULL;
    paramMan = NULL;
    randEng = NULL;
    randMan = NULL;
}
void HGenManager::addStdHistos()
{
    hbookMan->addHBook1D(&bookingEvent->getStruct()->qsq, &bookingEvent->getStruct()->USERVAR.at(2)  ,160, 0, 80, "Q^{2}-Histogram");
    hbookMan->addHBook1D(&bookingEvent->getStruct()->nu,&bookingEvent->getStruct()->USERVAR.at(2),100, 0, 200, "Nu-Histogram");
    hbookMan->addHBook1D(&bookingEvent->getStruct()->nu,&bookingEvent->getStruct()->dummyW,100, 0, 200, "Nu-Histogram(noweight)");

    hbookMan->addHBook1D(&bookingEvent->getStruct()->xbj,&bookingEvent->getStruct()->USERVAR.at(2),50, 0, 0.5, "Xbj-Histogram");
    hbookMan->addHBook1D(&bookingEvent->getStruct()->wsq,&bookingEvent->getStruct()->USERVAR.at(2),100, 0, 300, "W^{2}-Histogram");
    hbookMan->addHBook1D(&bookingEvent->getStruct()->emiss,&bookingEvent->getStruct()->USERVAR.at(2),50, -5, 20, "Emiss-Histogram");
    hbookMan->addHBook1D(&bookingEvent->getStruct()->tprim,&bookingEvent->getStruct()->USERVAR.at(2),150, 0, 1.5, "tprim_all");
    hbookMan->addHBook1D(&bookingEvent->getStruct()->tprim,&bookingEvent->getStruct()->USERVAR.at(15),150, 0, 1.5, "tprim_BH");
    hbookMan->addHBook1D(&bookingEvent->getStruct()->tprim,&bookingEvent->getStruct()->USERVAR.at(16),150, 0, 1.5, "tprim_DVCS");
}



void HGenManager::enableBeamFile(string _fileName, HBeamType _meGusta)
{
    beamFile = new HBeamFile(_fileName);
    if (beamFile->good())
        enabledBeamFile = true;
    wantedBeamType = _meGusta;
}




void HGenManager::setSeed(int _seed)
{
    seed = _seed;
    randEng->setSeed(_seed,1);
}


void HGenManager::oneShot()
{
    myEvent->reset();
    generator->generateEvent();

    bookingEvent->copyFrom(myEvent);

    printDebug();
}


void HGenManager::saveLastEvent()
{
    outputMan->dumpEvent(bookingEvent);
    hbookMan->fill();
    if (bookingEvent->getStruct()->ddActive)
    {
        ddBookMan->fill(0);
        ddBookMan->fill(1);
    }
}


double HGenManager::startRun(const void (*_callBack)(int))
{
    int events = paramMan->getNumEvents();

    myEvent->reset();

    addStdHistos();

    if (enabledBeamFile)
    {
        generator->useBeamFile();
        //enable the histograms for the beamfile
        hbookMan->addHBook1D((double*)&(thisBeam.energy),&bookingEvent->getStruct()->dummyW,100,100.,200.,"Beamfile Energy");
        hbookMan->addHBook1D((double*)&(aux[0]),&bookingEvent->getStruct()->dummyW,100,-5.,5.,"Beamfile Position x (fine)");
        hbookMan->addHBook1D((double*)&(aux[0]),&bookingEvent->getStruct()->dummyW,100,-30.,30.,"Beamfile Position x (raw)");
        hbookMan->addHBook1D((double*)&(aux[1]),&bookingEvent->getStruct()->dummyW,100,-5.,5.,"Beamfile Position y (fine)");
        hbookMan->addHBook1D((double*)&(aux[1]),&bookingEvent->getStruct()->dummyW,100,-30.,30.,"Beamfile Position y (raw)");
        hbookMan->addHBook1D((double*)&(aux[2]),&bookingEvent->getStruct()->dummyW,100,0.,1.,"Beamfile PT");
    }
    
    addGenHistos();
    
    //--- main event loop
    for (int i =0; i < events; i++)
    {
      _callBack(i);



        myEvent->reset();
        if (enabledBeamFile)
        {
            thisBeam = beamFile->getNextEntry(wantedBeamType);
            paramMan->getStruct()->ELEPT = thisBeam.energy;
            paramMan->getStruct()->PARL.at(2) = thisBeam.energy;
            generator->setNextBeam(thisBeam);
            aux[0] = thisBeam.position.X();
            aux[1] = thisBeam.position.Y();
            aux[2] = thisBeam.momentum.getPtrans();
        }
        generator->generateEvent();
        bookingEvent->copyFrom(myEvent);
        saveLastEvent();
	//cout << getParamManager()->getStruct()->LST[24] << endl;
        printDebug();
    }

    closeOutPut();
    
   return generator->getWeightSum();


}

void HGenManager::closeOutPut()
{
    string _histFileName;
    if (paramMan->getOutName() == "UNSET")
    {
        _histFileName = "HEPGen++-histograms-";
        _histFileName.append(hephelp::currentDateTime());
        _histFileName.append("-");
        _histFileName.append(hephelp::IntToStr(seed));
    }
    else
    {
        _histFileName = paramMan->getOutName() + "_histos";
    }

     hbookMan->dumpToFile(_histFileName);
    _histFileName.append("_dd");
     ddBookMan->dumpToFile(_histFileName);
    outputMan->closeFile();
}



void HGenManager::enableLeptoOutput(bool _gfortran)
{

    outputMan->enableLEPTO(_gfortran);
}

void HGenManager::enableOutPut(string _filename)
{
    outputMan->setFileName(_filename);
    outputMan->initFile(bookingEvent);
    outputMan->writeHeader();
}



HGenManager::~HGenManager(void )
{
    if (hbookMan != NULL)
        delete hbookMan;
    if (outputMan != NULL)
        delete outputMan;
    if (randEng != NULL)
        delete randEng;
    if (randMan != NULL)
        delete randMan;
}


void HGenManager::enableASCIIBook()
{
    enabledASCIIBook = true;
    hbookMan->addASCII();
    ddBookMan->addASCII();
}

void HGenManager::enableROOTBook()
{
    enabledROOTBook = true;
    hbookMan->addROOT();
    ddBookMan->addROOT();
}


void HGenManager::printDebug()
{
    if (!debug)
        return;

    cout << "------------------------------------ NEW EVENT ------------------------------------------" << endl;
    bookingEvent->printDebug();
    cout << endl << "------- Full Particle Listing: " << endl;
    for (unsigned int i =0; i < bookingEvent->getStruct()->listOfParticles.size(); i++)
    {
        cout << "p" << i+1 << ":\t";
        bookingEvent->getStruct()->listOfParticles.at(i)->printDebug();
    }
    if (bookingEvent->getStruct()->ddActive && bookingEvent->getStruct()->listOfParticles.size()>8)
    {
        if (bookingEvent->getStruct()->amx2 > 4)
            cout << "dd: LONGITUDINAL" << endl;
        else
            cout << "dd: ISOTROPIC" << endl;
    }
}

