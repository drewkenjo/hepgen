#include "hVGGgen.h"
#define DEBUG 0
void HPhysicsGenDVCSVGG::generateEvent()
{
    bool eventOK = false;
    //HPhysicsGen::generateEvent();
    //generate Nu
    while (!eventOK)
    {

        event->getStruct()->listOfParticles.resize(6);
        event->getStruct()->recoil.setParticleAuxFlag(1);



        //event->reset();
        if (!generateNu())
            continue;


        //this will set the beam parameters in the future when beamfilereading is implemented
        if (!setBeam())
            continue;


        //generate Q^2
        if (!generateQSQ())
            continue;

        //meson mass to 0
        if (!generateMesonMass())
            continue;


        if(!generatePhiGamma())
            continue;

        //maybe this will be implemented in the futute
        if (!generateElastic())
            continue;

        //generate the smearing
        if (!generateSmearing())
            continue;

        //generate the mandelstam t
        if (!generatet())
            continue;

        //generate the gamma kinematics
        if (!generateOutgoingParticle())
            continue;

        event->getStruct()->outPart1.setParticleType(hepconst::typeGamma);
        event->getStruct()->outPart1_Lab.setParticleType(hepconst::typeGamma);


        //calc the pt and transform outgoing vectors to lab system
        event->calcPT();



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




        //calculate the phi in DVCS
        if (!calculatePhir())
            continue;


        //calculate the weights of BH, dvcs
        calcWeights();


        if (doDD)
            ddGotoLab();

        eventOK = true;


    }


    event->rotXToZ();


    setParameters();
    setUserVars();


    generateListOfParticles();

    printf("Event finished! (just for your curiosity about speed ;) \n");

    if (DEBUG == 1)
    {
        event->getStruct()->listOfParticles.at(0)->printDebugHeader();
        for (int i =0; i < event->getStruct()->listOfParticles.size(); i ++ )
        {
            cout << i+1 << " ";
            event->getStruct()->listOfParticles.at(i)->printDebug();
        }
    }

}


void HPhysicsGenDVCSVGG::addHistograms()
{
    HPhysicsGen::addHistograms();
    bookMan->addHBook1D(&weightresult0,&event->getStruct()->dummyW,100,0,1,"Weightresult0 (DVCS)");
    bookMan->addHBook1D(&weightresult1,&event->getStruct()->dummyW,100,0,1,"Weightresult1 (BH)");
    bookMan->addHBook1D(&weightresult2,&event->getStruct()->dummyW,100,-100,100,"Weightresult2 (Interference)");
    bookMan->addHBook1D(&phir,&weight,100,0,6.5,"Phi_r (weight)");
    bookMan->addHBook1D(&phir,&none,100,0,6.5,"Phi_r (no weight)");
    bookMan->addHBook1D(&phir,&weightresult0,100,0,6.5,"Phi_r (weight-DVCS)");
    bookMan->addHBook1D(&phir,&weightresult1,100,0,6.5,"Phi_r (weight-BH)");
    bookMan->addHBook1D(&phir,&weightresult2,100,0,6.5,"Phi_r (weight-Interference)");
    bookMan->addHBook1D(&phi_out, &weight,100,0,6.5,"Phi_out");
    bookMan->addHBook1D(&theta_out, &weight,100,0,3.5,"Theta_out");
    bookMan->addHBook1D(&theta_gamma_gammavirt, &weightresult1,300,0,0.3,"Theta_GAMMA_GAMMAVIRT(BH weight)");
    bookMan->addHBook1D(&event->getStruct()->tprim,&event->getStruct()->dummyW,100,0,1,"TPRIM_NOWEIGHT");
}


/* maybe this will be added later, lets first see if it works without */
void HPhysicsGenDVCSVGG::setParameters()
{
    paramMan->getStruct()->PARL.at(6)=hepconst::typeGamma;
    //no decays here so its just 0 and 0 for the decay particles
    paramMan->getStruct()->PARL.at(7)=0.;
    paramMan->getStruct()->PARL.at(8)=0.;
}

bool HPhysicsGenDVCSVGG::generateMesonMass()
{
    m_meson = 0.0;
    //this sets the PARLD correctly
    HPhysicsGen::generateMesonMass();
    return true;
}



HPhysicsGenDVCSVGG::HPhysicsGenDVCSVGG(CLHEP::HepRandom* _randGen, HParamManager* _paramMan, HEvent* _event, HBooker* _booker, HBooker* _ddBooker) : HPhysicsGen(_randGen, _paramMan, _event, _booker, _ddBooker)
{
    printf("CAUTION! Ultra-Experimental VGG-DVCS-Generator loading! This might be slow as hell!!!!! \n");
    std::cout << "VGG-DVCS-Generator loading!! " << std::endl;
    myRandom = _randGen;
    paramMan = _paramMan;
    event = _event;
    bookMan = _booker;
    GaussEngine = new CLHEP::RanluxEngine();
    myRandomGauss = new CLHEP::RandGauss(GaussEngine);
    bookMan->addHBook1D(&phir,&event->getStruct()->dummyW,100,0,2*M_PI,"Phi_r");



    none=1;


    //for (int i =1;i<10;i++)
    //  cout << "random num "<<i<<" " << myRandom->flat()<<endl;

    event->getStruct()->type = hepconst::DVCS;



    //build the list
    event->getStruct()->listOfParticles.clear();
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->incBeamParticle);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->targetParticle);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->gammaVirt);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->scatBeamParticle);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->outPart1_Lab);
    event->getStruct()->listOfParticles.push_back(&event->getStruct()->recoil);



    //now set the sharp-mass-flags
    event->getStruct()->listOfParticles.at(0)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(1)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(2)->setSharpMass(false);
    event->getStruct()->listOfParticles.at(3)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(4)->setSharpMass(true);
    event->getStruct()->listOfParticles.at(5)->setSharpMass(true);



    //now we set the origins - these mean the (index-number+1) of the particle from which the particle originates

    //target particles dont have an origin
    event->getStruct()->listOfParticles.at(0)->setParticleOrigin(0);
    event->getStruct()->listOfParticles.at(1)->setParticleOrigin(0);
    event->getStruct()->listOfParticles.at(2)->setParticleOrigin(1);
    event->getStruct()->listOfParticles.at(3)->setParticleOrigin(1);
    event->getStruct()->listOfParticles.at(4)->setParticleOrigin(3);
    event->getStruct()->listOfParticles.at(5)->setParticleOrigin(2);



    //set the aux-flags of the particles, 3 are dead alreadz 3 are living
    event->getStruct()->listOfParticles.at(0)->setParticleAuxFlag(21);
    event->getStruct()->listOfParticles.at(1)->setParticleAuxFlag(21);
    event->getStruct()->listOfParticles.at(2)->setParticleAuxFlag(21);
    event->getStruct()->listOfParticles.at(3)->setParticleAuxFlag(1);
    event->getStruct()->listOfParticles.at(4)->setParticleAuxFlag(1);
    event->getStruct()->listOfParticles.at(5)->setParticleAuxFlag(1);

    //set the daughter-line-numbers
    event->getStruct()->listOfParticles.at(0)->setParticleDaughter1(3);
    event->getStruct()->listOfParticles.at(1)->setParticleDaughter1(6);
    event->getStruct()->listOfParticles.at(2)->setParticleDaughter1(5);
    event->getStruct()->listOfParticles.at(3)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(4)->setParticleDaughter1(0);
    event->getStruct()->listOfParticles.at(5)->setParticleDaughter1(0);

    //set the daughter-line-numbers #2
    event->getStruct()->listOfParticles.at(0)->setParticleDaughter2(4);
    event->getStruct()->listOfParticles.at(1)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(2)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(3)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(4)->setParticleDaughter2(0);
    event->getStruct()->listOfParticles.at(5)->setParticleDaughter2(0);
}





void HPhysicsGenDVCSVGG::setUserVars()
{
    HPhysicsGen::setUserVars();
    event->getStruct()->USERVAR.at(2)=weight;
    event->getStruct()->USERVAR.at(15)=weightresult0;
    event->getStruct()->USERVAR.at(16)=weightresult1;

}



void HPhysicsGenDVCSVGG::generateListOfParticles()
{



    //we set the particle IDs here and make the list ready.
    event->getStruct()->incBeamParticle.setParticleType(paramMan->getBeamPart());
    event->getStruct()->scatBeamParticle.setParticleType(paramMan->getBeamPart());
    event->getStruct()->gammaVirt.setParticleType(hepconst::typeGamma);
    event->getStruct()->outPart1_Lab.setParticleType(hepconst::typeGamma);




}


HPhysicsGenDVCSVGG::~HPhysicsGenDVCSVGG()
{
    delete GaussEngine;
    delete myRandomGauss;
}



void HPhysicsGenDVCSVGG::calcWeights(hWeightInterface _in, double* BH, double* DVCS, double* INT, bool _wantLO)
{
    map<string,double> paraMap;
    map<string,double> weightMap;
    paraMap["XBJ"]=_in.xbj;
    paraMap["Q2"] =_in.qsq;
    paraMap["T"]=_in.t;
    paraMap["MU2"]=2.5; // fixed
    paraMap["PHIR"]=_in.phir;
    paraMap["ELEPT"]=_in.beamE;
    paraMap["BEAMCHARGE"]=_in.clept;
    paraMap["BEAMHELICITY"]=-_in.clept;

    if (_wantLO)
      compassInterFace::getInstance()->getCXLO(weightMap,paraMap);
    else
      compassInterFace::getInstance()->getCX(weightMap,paraMap);
    
    *DVCS = weightMap["DVCS"];
    *BH =  weightMap["BH"];
    *INT =  weightMap["INTERFERENCE"];
}




void HPhysicsGenDVCSVGG::calcWeights()
{



    if (phir < 0)
        phir += 2* M_PI;


    phir_trento = phir;
    //phir = M_PI-phir;


    map<string,double> paraMap;
    map<string,double> weightMap;
    paraMap["XBJ"]=event->getXbj();
    paraMap["Q2"] =event->getQsq();
    paraMap["T"]=event->getT();
    paraMap["MU2"]=1.0; // fixed
    paraMap["PHIR"]=phir_trento;
    paraMap["ELEPT"]=paramMan->getBeamE();
    paraMap["BEAMCHARGE"]=paramMan->getclept();
    paraMap["BEAMHELICITY"]=paramMan->getslept();
    compassInterFace::getInstance()->getCX(weightMap,paraMap);


    //set it in PARL-struct
    event->getStruct()->PARL.at(27) = phir;
    event->getStruct()->PARL.at(28) = 0.;
    event->getStruct()->PARL.at(29) = 0.;

    if (weightMap["DVCS"]!=weightMap["DVCS"] || weightMap["BH"]!=weightMap["BH"] || weightMap["INTERFERENCE"]!=weightMap["INTERFERENCE"] || weightMap["CONVERSIONFACTOR"]!=weightMap["CONVERSIONFACTOR"])
    {
        weightresult0 = 0;
        weightresult1 = 0;
        weightresult2 = 0;
        weight =0;
    }
    else
    {

        weightresult0 = weightMap["DVCS"];
        weightresult1 =  weightMap["BH"];
        weightresult2 =  weightMap["INTERFERENCE"];
        weight = weightresult0+weightresult1+weightresult2;
    }




    //TODO regenerate the event, but this is beeing kept for historical reasons
    if (weight < 0)
    {
        weightbh =0;
        weightint = 0;
        weightdvcs=0;
        weightresult0=0;
        weightresult1=0;
        weightresult2=0;
        weight = 0;
    }
    if (!doDD)
      weightCounter+=weight;

}





