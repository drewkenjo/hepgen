#include "hpigen.h"
#include "hgenmanager.h"
#include "reweightKine.h"
#define DEBUG 1

void HPhysicsGenPI::generateEvent()
{
    //HPhysicsGen::generateEvent();
    bool eventOK = false;
    //HPhysicsGen::generateEvent();
    //generate Nu
    while (!eventOK)
    {
        event->getStruct()->listOfParticles.resize(8);
        event->getStruct()->recoil.setParticleAuxFlag(1);
        //generate Nu
        if (!generateNu())
            continue;

        //this will set the beam parameters in the future when beamfilereading is implemented
        if(!setBeam())
            continue;


        //generate Q^2
        if(!generateQSQ())
            continue;

        //meson mass to rho and a bit randomness
        if(!generateMesonMass())
            continue;


        if(!generatePhiGamma())
            continue;

        //maybe this will be implemented in the futute
        if(!generateElastic())
            continue;

        //generate the smearing
        if(!generateSmearing())
            continue;
        //generate the mandelstam t
        if(!generatet())
            continue;
        //generate the gamma kinematics
        if(!generateOutgoingParticle())
            continue;
        //calc the pt and transform outgoing vectors to lab system
        event->calcPT();


        generateDecay();

        if(!calculatePhir())
            continue;

        //here we need to get
        calcWeights();


        /*do some special stuff here for histogramming */

        //get the theta between outgoing and gamma*
        HLorentzVector gammaVirtLab =  event->goToLabSystem(event->getStruct()->gammaVirt.getVector(),event->getStruct()->CMS);

        //get 3-outvectors in labsystem
        HVector3 gammaVirtLabThree = gammaVirtLab.getVector();
        HVector3 outPart1LabThree = event->getStruct()->outPart1_Lab.getVector().getVector();

        //normalize them, we only want some angles!
        gammaVirtLabThree.normalize(1.0);
        outPart1LabThree.normalize(1.0);

        //get the scalar product of them

        theta_gamma_gammavirt = acos(gammaVirtLabThree.dotProduct(outPart1LabThree));

        /* done with special histogramming stuff */


        if (doDD)
            ddGotoLab();
        eventOK = true;
	if (!doDD)
	    weightCounter += weight;

    }
    event->rotXToZ();
    if (enabledBeamFile)
        event->rotateEventToBeam(beamEntry.momentum.getVector());



    if (DEBUG == 1) {

        hWeightInterface myWeights,myRoundedWeights;
        hreweightKine::setDefaultProductionValues(myWeights);
        hreweightKine::setReggeParams(myWeights);

        hreweightKine::setDefaultProductionValues(myRoundedWeights);
        hreweightKine::setReggeParams(myRoundedWeights);


        HLorentzVector incMuon = event->getBeam().getVector();
        HLorentzVector scatMuon = event->getScat().getVector();
        HLorentzVector gReal = event->getOutPart1_Lab().getVector();
        HLorentzVector pRecoil = event->getRecoil().getVector();



        myWeights.beamE = incMuon.getEnergy();
        myWeights.MuIn = incMuon;
        myWeights.qsq = event->getQsq();
        myWeights.t = event->getStruct()->t;
        myWeights.tprim = event->getStruct()->tprim;
        myWeights.w2 = event->getWsq();
        myWeights.nu = event->getStruct()->nu;
        myWeights.y = event->getStruct()->y;

// 	myWeights.qsqmin = 0.7;
// 	myWeights.qsqmax = 12;
// 	myWeights.numin = 5;
// 	myWeights.numax = 155;
// 	myWeights.tmax = 1.0;
// 	myWeights.tmin = 0.001;

        HPhysicsGenPI::calcWeights(myWeights,newWeightUnRounded,myPionic);
        incMuon.roundToFloat();
        myWeights.beamE = static_cast<float>(incMuon.getEnergy());
        myWeights.MuIn = incMuon;
        myWeights.qsq = static_cast<float>(event->getQsq());
        myWeights.t = static_cast<float>(event->getStruct()->t);
        myWeights.tprim = static_cast<float>(event->getStruct()->tprim);
        myWeights.w2 = static_cast<float>(event->getWsq());
        myWeights.nu = static_cast<float>(event->getStruct()->nu);
        myWeights.y = static_cast<float>(event->getStruct()->y);
        HPhysicsGenPI::calcWeights(myWeights,newWeightRounded,myPionic);

        newWeightRounded -= weight;
        newWeightRounded /= weight;
        newWeightUnRounded -= weight;
        newWeightUnRounded /= weight;
        newWeightRoundedAbs = fabs(newWeightRounded);





    }




    setParameters();
    setUserVars();


    generateListOfParticles();


}


void HPhysicsGenPI::addHistograms()
{
    bookMan->addHBook1D(&weight,&event->getStruct()->dummyW,100,0,0.2,"weight(Pi0)");
    bookMan->addHBook1D(&phi_out, &weight,20,0,2*M_PI,"Phi_out");
    bookMan->addHBook1D(&thetapi, &weight,100 ,0 , 3.5 , "Theta_DP1");
    bookMan->addHBook1D(&phipi, &weight,100 , 0 , 6.5 , "Phi_DP1");
    bookMan->addHBook1D(&costhetapi, &weight,50,-1,1,"Cos (Theta_DP)");
    bookMan->addHBook1D(&phir,&weight,600,0,2*M_PI,"Phi_r");
    bookMan->addHBook1D(&phi_out, &weight,100,0,6.5,"Phi_out");
    bookMan->addHBook1D(&theta_out, &weight,100,0,3.5,"Theta_out");

    bookMan->addHBook1D(&theta_gamma_gammavirt, &weight,600,-3.5,7.0,"Theta_GAMMA_GAMMAVIRT(pi-Weight)");


    bookMan->addHBook1D(&newWeightRounded,&weight,1000,-0.2,0.2,"Rounded vectors PI - relative weight diff");
    bookMan->addHBook1D(&newWeightRoundedAbs,&weight,1000,0,0.2,"Rounded vectors PI - relative weight diff (abs)");

    bookMan->addHBook1D(&newWeightUnRounded,&weight,1000,-0.2,0.2,"Original double precision vectors - PI relative weight diff");

//     bookMan->addHBook2D(&newWeightRounded,&weight,&weight,1000,1000,0,0,0.2,1000,"PI - Rounded errors weighted vs weight");

}

void HPhysicsGenPI::printAuxVars()
{
    HPhysicsGen::printAuxVars();
    cout << "'AUX1 \t /path/to/pion_table.dat'" << endl;
    cout << "'AUX2 \t -1 for autodetect LEGACY: 1' or 0 for table with transverse components" << endl;
}




HPhysicsGenPI::HPhysicsGenPI(CLHEP::HepRandom* _randGen, HParamManager* _paramMan, HEvent* _event, HBooker* _booker, HBooker* _ddBooker) : HPhysicsGen(_randGen, _paramMan, _event, _booker, _ddBooker)
{
    std::cout << "Pion0-Generator loading!! " << std::endl;
    myRandom = _randGen;
    paramMan = _paramMan;
    event = _event;
    bookMan = _booker;
    GaussEngine = new CLHEP::RanluxEngine();
    myRandomGauss = new CLHEP::RandGauss(GaussEngine);

    string piFileName = paramMan->getStruct()->aux1.at(1);
    myPionic = new HPionicData();
    paramMan->setPionicData(myPionic);

    if (paramMan->getStruct()->aux2.at(1) == "-1") {
        printf("Trying new load-method!\n");
        myPionic->loadFileNG(piFileName);
    }
    else {

	printf("Legacy table loading requested!!\n");
        
        double pi0w[11] =  {5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.};
        double pi0tpr[16] =  {0.,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75};




        myPionic->getBinningW()->insert(myPionic->getBinningW()->end(),&pi0w[0], &pi0w[11]);
        myPionic->getBinningTPR()->insert(myPionic->getBinningTPR()->end(),&pi0tpr[0], &pi0tpr[16]);


        if (paramMan->getStruct()->aux1.size() > 1 && paramMan->getStruct()->aux2.at(1) != "1")
        {
            cout << "Using standard pi0-table" << endl;
            //set the binning for the data-table
            double pi0qsq[9] =  {2.,3.,4.,5.,6.,7.,8.,9.,10.};
            myPionic->getBinningQ2()->insert(myPionic->getBinningQ2()->end(),&pi0qsq[0],&pi0qsq[9]);

        }
        else if (paramMan->getStruct()->aux1.size() > 1 && paramMan->getStruct()->aux2.at(1) != "0")
        {

            cout << "Using pi0-table with transverse E influence" << endl;
            //set the binning for the data-table
            double pi0qsq[15] =  {2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.};
            myPionic->getBinningQ2()->insert(myPionic->getBinningQ2()->end(),&pi0qsq[0],&pi0qsq[15]);

        }

        myPionic->resetData();
        myPionic->loadFile(piFileName);
	
    }
    event->getStruct()->type = hepconst::PI0;


    //build the list
    event->getStruct()->listOfParticles.clear();
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->incBeamParticle);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->targetParticle);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->gammaVirt);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->scatBeamParticle);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart1_Lab);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart2_Lab);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart3_Lab);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->recoil);


    event->getStruct()->listOfParticles.at(0)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(1)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(2)->setSharpMass(false);
    event->getStruct()->listOfParticles.at(3)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(4)->setSharpMass(false);
    event->getStruct()->listOfParticles.at(5)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(6)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(7)->setSharpMass(true);


    //now we set the origins - these mean the (index-number+1) of the particle from which the particle originates

    //target particles dont have an origin
    event->getStruct()->listOfParticles.at(0)->setParticleOrigin(0);
    event->getStruct()->listOfParticles.at(1)->setParticleOrigin(0);
    event->getStruct()->listOfParticles.at(2)->setParticleOrigin(1);
    event->getStruct()->listOfParticles.at(3)->setParticleOrigin(1);
    event->getStruct()->listOfParticles.at(4)->setParticleOrigin(3);
    event->getStruct()->listOfParticles.at(5)->setParticleOrigin(5);
    event->getStruct()->listOfParticles.at(6)->setParticleOrigin(5);
    event->getStruct()->listOfParticles.at(7)->setParticleOrigin(2);




    //set the aux-flags of the particles, 3 are dead already 1 is decayed rest is aLivE!
    event->getStruct()->listOfParticles.at(0)->setParticleAuxFlag(21);
    event->getStruct()->listOfParticles.at(1)->setParticleAuxFlag(21);
    event->getStruct()->listOfParticles.at(2)->setParticleAuxFlag(21);
    event->getStruct()->listOfParticles.at(3)->setParticleAuxFlag(1);
    event->getStruct()->listOfParticles.at(4)->setParticleAuxFlag(11);
    event->getStruct()->listOfParticles.at(5)->setParticleAuxFlag(1);
    event->getStruct()->listOfParticles.at(6)->setParticleAuxFlag(1);
    event->getStruct()->listOfParticles.at(7)->setParticleAuxFlag(1);


    //set the daughter-line-numbers
    event->getStruct()->listOfParticles.at(0)->setParticleDaughter1(3);
    event->getStruct()->listOfParticles.at(1)->setParticleDaughter1(8);
    event->getStruct()->listOfParticles.at(2)->setParticleDaughter1(5);
    event->getStruct()->listOfParticles.at(3)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(4)->setParticleDaughter1(6);
    event->getStruct()->listOfParticles.at(5)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(6)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(7)->setParticleDaughter1(0);


    //set the daughter-line-numbers #2
    event->getStruct()->listOfParticles.at(0)->setParticleDaughter2(4);
    event->getStruct()->listOfParticles.at(1)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(2)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(3)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(4)->setParticleDaughter2(7);
    event->getStruct()->listOfParticles.at(5)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(6)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(7)->setParticleDaughter2(0);



    //for (int i =1;i<10;i++)
    //  cout << "random num "<<i<<" " << myRandom->flat()<<endl;

}

bool HPhysicsGenPI::generateMesonMass()
{
    m_meson = sqrt(hepconst::w2pi);
    HPhysicsGen::generateMesonMass();
    return true;
}


void HPhysicsGenPI::setParameters()
{
    paramMan->getStruct()->PARL.at(6)=hepconst::typePi0;
    paramMan->getStruct()->PARL.at(7)=hepconst::typeGamma;
    paramMan->getStruct()->PARL.at(8)=hepconst::typeGamma;

}




void HPhysicsGenPI::generateDecay()
{

    // first we generate the angles of the first pion
    double r = myRandom->flat();
    double costh;

    //generate theta flat
    costh = 2*r - 1;

    //generate phi now
    r = myRandom->flat();
    phipi=r*2*M_PI;

    // then we generate the particle out of the angles
    thetapi = acos (costh);





    HLorentzVector decayPi1_cms;
    //now we set it in the CMS
    //momentum of the pions
    double mrhosq = pow(m_meson,2);
    double pimom = mrhosq/(2*m_meson);
    //generate the lorentz vector
    decayPi1_cms.setLVectorAngular(pimom,thetapi,phipi,pimom);
    HLorentzVector decayPi1_lab = event->goToLabSystem(decayPi1_cms,event->getOutPart1_Lab().getVector());

    //fill it in the event - positive pion
    event->getStruct()->outPart2.setVector(decayPi1_cms);
    event->getStruct()->outPart2.setParticleType(hepconst::typeGamma);

    event->getStruct()->outPart2_Lab.setVector(decayPi1_lab);
    event->getStruct()->outPart2_Lab.setParticleType(hepconst::typeGamma);


    //checked until here


    //now we generate the negative pion by just inverting the 3-vector of momentum
    HVector3 pineg_mom = HVector3(0.,0.,0.)-decayPi1_cms.getVector();
    HLorentzVector decayPi2_cms;
    decayPi2_cms.setVector(pineg_mom);
    decayPi2_cms.setEnergy(decayPi1_cms.getEnergy());

    HLorentzVector decayPi2_lab = event->goToLabSystem(decayPi2_cms,event->getOutPart1_Lab().getVector());

    //event->getOutPart1_Lab().getVector().print();
    event->getStruct()->outPart3.setVector(decayPi2_cms);
    event->getStruct()->outPart3.setParticleType(hepconst::typeGamma);

    event->getStruct()->outPart3_Lab.setVector(decayPi2_lab);
    event->getStruct()->outPart3_Lab.setParticleType(hepconst::typeGamma);


}

void HPhysicsGenPI::calcWeights(hWeightInterface _in, double& WEIGHT,HPionicData* dataTable)
{
    if (dataTable == NULL)
        return;

    int binW = dataTable->getBinW(sqrt(_in.w2));
    int binQsq = dataTable->getBinQsq(_in.qsq);
    int binTprim = dataTable->getBinTPrim(_in.tprim);

//     printf("BINS %i %i %i \n",binW,binQsq,binTprim);

    //int binW1,binW2,binQsq1,binQsq2,binTprim1,binTprim2;
    int signW=1,signQsq=1,signTprim=1;

    double delta[3];
    delta[0] = sqrt(_in.w2)-dataTable->getPi0w(binW);
    delta[1] = _in.qsq-dataTable->getPi0Qsq(binQsq);
    delta[2] = _in.tprim - dataTable->getPi0tpr(binTprim);

    if (delta[0] < 0)
        signW = -1;
    if (binW == 0)
        signW = 1;
    if (binW == dataTable->getBinningW()->size()-1)
        signW = -1;

    if (delta[1] < 0)
        signQsq = -1;
    if (binQsq == 0)
        signQsq = 1;
    if (binQsq == dataTable->getBinningQ2()->size()-1)
        signQsq = -1;

    if (delta[2] < 0)
        signTprim = -1;
    if (binTprim == 0)
        signTprim = 1;
    if (binTprim == dataTable->getBinningTPR()->size()-1)
        signTprim = -1;


    double dersl[3];
    double derst[3];

    dersl[0] = (dataTable->getPi0sigl(binW + signW ,binQsq,binTprim)-dataTable->getPi0sigl(binW,binQsq,binTprim))*signW/1.0;
    dersl[1] = (dataTable->getPi0sigl(binW,binQsq+signQsq,binTprim)-dataTable->getPi0sigl(binW,binQsq,binTprim))*signQsq/1.0;
    dersl[2] = (dataTable->getPi0sigl(binW,binQsq,binTprim+signTprim)-dataTable->getPi0sigl(binW,binQsq,binTprim))*signTprim/0.05;




    double sigl = dataTable->getPi0sigl(binW,binQsq,binTprim) + dersl[0]*delta[0] + dersl[1]*delta[1] + dersl[2] * delta[2];

//     printf("dersl %f, %f, %f \n",dersl[0],dersl[1],dersl[2]);


    derst[0] = (dataTable->getPi0sigt(binW + signW,binQsq,binTprim)-dataTable->getPi0sigt(binW,binQsq,binTprim))*signW/1.0;
    derst[1] = (dataTable->getPi0sigt(binW,binQsq + signQsq,binTprim)-dataTable->getPi0sigt(binW,binQsq,binTprim))*signQsq/1.0;
    derst[2] = (dataTable->getPi0sigt(binW,binQsq,binTprim+ signTprim)-dataTable->getPi0sigt(binW,binQsq,binTprim))*signTprim/0.05;
    double sigt = dataTable->getPi0sigt(binW,binQsq,binTprim) + derst[0]*delta[0] + derst[1]*delta[1] + derst[2] * delta[2];
//     printf("derst %f, %f, %f \n",derst[0],derst[1],derst[2]);



    double siggam = sigt + hreweightKine::getEpsilon(_in) * sigl;


    //needs to be multiplied by total Phase Factor!
    WEIGHT = siggam * hreweightKine::getFluxCompensator(_in);


    if (WEIGHT < 0)
        WEIGHT = 0;

}



void HPhysicsGenPI::calcWeights()
{
    HPionicData* pionicdata = paramMan->getPionicData();

    int binW = pionicdata->getBinW(sqrt(event->getWsq()));
    int binQsq = pionicdata->getBinQsq(event->getQsq());
    int binTprim = pionicdata->getBinTPrim(event->getStruct()->tprim);

    int signW=1,signQsq=1,signTprim=1;

    double delta[3];
    delta[0] = sqrt(event->getWsq())-pionicdata->getPi0w(binW);
    delta[1] = event->getQsq()-pionicdata->getPi0Qsq(binQsq);
    delta[2] = event->getStruct()->tprim - pionicdata->getPi0tpr(binTprim);

    if (delta[0] < 0)
        signW = -1;
    if (binW == 0)
        signW = 1;
    if (binW == pionicdata->getBinningW()->size()-1)
        signW = -1;

    if (delta[1] < 0)
        signQsq = -1;
    if (binQsq == 0)
        signQsq = 1;
    if (binQsq == pionicdata->getBinningQ2()->size()-1)
        signQsq = -1;

    if (delta[2] < 0)
        signTprim = -1;
    if (binTprim == 0)
        signTprim = 1;
    if (binTprim == pionicdata->getBinningTPR()->size()-1)
        signTprim = -1;


    double dersl[3];
    double derst[3];


    dersl[0] = (pionicdata->getPi0sigl(binW + signW ,binQsq,binTprim)-pionicdata->getPi0sigl(binW,binQsq,binTprim))*signW/1.0;
    dersl[1] = (pionicdata->getPi0sigl(binW,binQsq+signQsq,binTprim)-pionicdata->getPi0sigl(binW,binQsq,binTprim))*signQsq/1.0;
    dersl[2] = (pionicdata->getPi0sigl(binW,binQsq,binTprim+signTprim)-pionicdata->getPi0sigl(binW,binQsq,binTprim))*signTprim/0.05;




    double sigl = pionicdata->getPi0sigl(binW,binQsq,binTprim) + dersl[0]*delta[0] + dersl[1]*delta[1] + dersl[2] * delta[2];



    derst[0] = (pionicdata->getPi0sigt(binW + signW,binQsq,binTprim)-pionicdata->getPi0sigt(binW,binQsq,binTprim))*signW/1.0;
    derst[1] = (pionicdata->getPi0sigt(binW,binQsq + signQsq,binTprim)-pionicdata->getPi0sigt(binW,binQsq,binTprim))*signQsq/1.0;
    derst[2] = (pionicdata->getPi0sigt(binW,binQsq,binTprim+ signTprim)-pionicdata->getPi0sigt(binW,binQsq,binTprim))*signTprim/0.05;
    double sigt = pionicdata->getPi0sigt(binW,binQsq,binTprim) + derst[0]*delta[0] + derst[1]*delta[1] + derst[2] * delta[2];

    double siggam = sigt + event->getStruct()->epsilon * sigl;

    weight = event->getStruct()->flux * siggam * event->getTotalPhaseFactor();


    if (weight < 0)
        weight = 0;

    costhetapi = cos (theta_out);

}

void HPhysicsGenPI::generateListOfParticles()
{
    //we set the particle IDs here and make the list ready.
    event->getStruct()->incBeamParticle.setParticleType(paramMan->getBeamPart());
    event->getStruct()->scatBeamParticle.setParticleType(paramMan->getBeamPart());
    event->getStruct()->gammaVirt.setParticleType(hepconst::typeGamma);
    event->getStruct()->outPart1_Lab.setParticleType(hepconst::typePi0);
    event->getStruct()->outPart2_Lab.setParticleType(hepconst::typeGamma);
    event->getStruct()->outPart3_Lab.setParticleType(hepconst::typeGamma);



}


void HPhysicsGenPI::setUserVars()
{
    HPhysicsGen::setUserVars();
    event->getStruct()->USERVAR.at(2)=weight;
    if (phir < 0)
        phir += 2* M_PI;
    event->getStruct()->PARL.at(27) = phir;
    event->getStruct()->PARL.at(28) = thetapi;
    event->getStruct()->PARL.at(29) = phipi;
}






HPhysicsGenPI::~HPhysicsGenPI()
{
    delete GaussEngine;
    delete myRandomGauss;
}
