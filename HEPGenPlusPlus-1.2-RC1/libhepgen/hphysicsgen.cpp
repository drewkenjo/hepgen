#include "hphysicsgen.h"
#include "hparticle.h"




HPhysicsGen::HPhysicsGen(CLHEP::HepRandom* _randGen, HParamManager* _paramMan, HEvent* _event, HBooker* _booker, HBooker* _ddBookMan)
{
    enabledBeamFile = false;
    myRandom = _randGen;
    paramMan = _paramMan;
    event = _event;
    bookMan = _booker;
    ddBookMan = _ddBookMan;
    GaussEngine = new CLHEP::RanluxEngine();
    myRandomGauss = new CLHEP::RandGauss(GaussEngine);
    weightCounter = 0.0;

}



void HPhysicsGen::addHistograms()
{

    bookMan->addHBook1D(&event->getStruct()->qsq,&pftqsq_,100,0.5,100,"qsq pf weighted");
    bookMan->addHBook1D(&event->getStruct()->tprim,&pft_,100,0.001,1,"tprime pf weighted");
    bookMan->addHBook1D(&event->getStruct()->qsq,&event->getStruct()->dummyW,100,0.5,100,"qsq non weighted");
    bookMan->addHBook1D(&event->getStruct()->tprim,&event->getStruct()->dummyW,100,0.001,1,"tprime non weighted");
    bookMan->addHBook1D(&event->getStruct()->nu,&event->getStruct()->dummyW,200,0.0,200,"nu non weighted");
    bookMan->addHBook1D(&event->getStruct()->nu,&pfnu_,200,0.0,200,"nu pf weighted");
    
    
    
    bookMan->addHBook2D((double*)&event->getStruct()->qsq,(double*)&event->getStruct()->nu,&event->getStruct()->dummyW,200,160,0.0,0.,100.,160.,"qsq vs nu unweighted");
    bookMan->addHBook2D((double*)&event->getStruct()->qsq,(double*)&event->getStruct()->nu,&event->getStruct()->USERVAR[9],200,160,0.0,0.,100.,160.,"qsq vs nu old pfweighted");
    
    bookMan->addHBook2D((double*)&event->getStruct()->qsq,(double*)&event->getStruct()->nu,&pfoldqsq,200,160,0.0,0.,100.,160.,"qsq vs nu pfweighted");
    
    


    ddBookMan->addHBook1D((double*)&ddIterations,&event->getStruct()->dummyW,100,0.,100.,"Iterations of ddDecayParticle");
    ddBookMan->addHBook1D((double*)&event->getStruct()->amx2,&event->getStruct()->USERVAR.at(2),600,0.,40.,"AMX2 (mass square of excited state)");


    ddBookMan->addHBook1D((double*)&event->getStruct()->ddNum_Charged,&event->getStruct()->dummyW,20,0.,20.,"charged ddParticles");
    ddBookMan->addHBook1D((double*)&event->getStruct()->ddNum_Neutral,&event->getStruct()->dummyW,20,0.,20.,"neutral ddParticles");
    ddBookMan->addHBook1D((double*)&event->getStruct()->ddNum_Particles,&event->getStruct()->dummyW,20,0.,20.,"total Particles");

    ddBookMan->addHBook2D((double*)&event->getStruct()->amx2,(double*)&event->getStruct()->ddPMax,&event->getStruct()->dummyW,60,30,0.,0.,30.,3.,"amx2 vs pmax");

    ddBookMan->addHBook2D((double*)&event->getStruct()->ddMX,(double*)&event->getStruct()->ddNum_Particles,&event->getStruct()->dummyW,60,15,0.,0.,6.,15.,"mx2 vs ntot");
    ddBookMan->addHBook2D((double*)&event->getStruct()->ddMX,(double*)&event->getStruct()->ddNum_Charged,&event->getStruct()->dummyW,60,15,0.,0.,6.,15.,"mx2 vs ncharged");
    ddBookMan->addHBook2D((double*)&event->getStruct()->ddMX,(double*)&event->getStruct()->ddNum_Neutral,&event->getStruct()->dummyW,60,15,0.,0.,6.,15.,"mx2 vs nneutral");
}

void HPhysicsGen::printAuxVars()
{

    cout << "Aux-Vars for this generator: " << endl;

}




bool HPhysicsGen::setBeam()
{
    if(enabledBeamFile)
    {
        //to clarify: we set the beam along the z axis with the given energy and rotate it afterwards. That is so, because some
        //of the kinematics code does not play well with vectors. One should fix this in future versions.
        event->getBeam().getVector() = HLorentzVector(sqrt(pow(beamEntry.energy,2)-hepconst::w2mu),0.,0.,beamEntry.energy);
        paramMan->getStruct()->ELEPT = beamEntry.energy;
        event->getStruct()->USERVAR.at(0) = beamEntry.position.X();
        event->getStruct()->USERVAR.at(1) = beamEntry.position.Y();
    }

    return true;
}


bool HPhysicsGen::calculatePhir()
{
    //cout << "running phidvcs" << endl;
    HLorentzVector pcmsl = event->getCMS();
    HLorentzVector ebeamlepton = event->getBeam().getVector();
    HLorentzVector pelscatteredlepton = event->getScat().getVector();
    HLorentzVector pgammavirt = event->getGammaVirt().getVector();
    HLorentzVector pgammaout = event->getOutPart1_Lab().getVector(); //we have a real gamma here
    HLorentzVector recoil = event->getRecoil().getVector();

    // transform these vectors to the CMS system again for historical reasons
    double u = sqrt(event->getWsq());



    ebeamlepton.lorenf(u,pcmsl);
    pelscatteredlepton.lorenf(u,pcmsl);
    pgammavirt.lorenf(u,pcmsl);
    pgammaout.lorenf(u,pcmsl);
    recoil.lorenf(u,pcmsl);





    HVector3 vscat = ebeamlepton.getVector().crossProduct(pelscatteredlepton.getVector());
    vscat.normalize(1.);


    HVector3 vprod = pgammavirt.getVector().crossProduct(pgammaout.getVector());
    vprod.normalize(1.);

    HVector3 vscpr = vscat.crossProduct(vprod);
    double scprod = vscpr.dotProduct(pgammaout.getVector());

    int sign = +1;
    if (scprod < 0)
        sign = -1;


    double sphir = ((double)sign) * vscpr.length();
    double cphir = vscat.dotProduct(vprod);


    phir = atan2(sphir,cphir);


    if (phir > - M_PI/2 && phir < M_PI/2)
        iphir = 1.0;
    else
        iphir = -1.0;

    return true;
}


bool HPhysicsGen::generateNu()
{
    //as this is the generation starting point, make sure the event is reset;
    event->reset();

    // get a new flat random out of our random generator and get a nu with it
    double r = myRandom->flat();
    double nu = (paramMan->getNuMax()-paramMan->getNuMin()) * r + paramMan->getNuMin();


    //calculate phase-factor
    double phaseFactor = 1./(paramMan->getNuMax()-paramMan->getNuMin());
    pfnu_ = 1./phaseFactor;


    return event->setNu(nu,phaseFactor);

}

void HPhysicsGen::generateEvent()
{
}




bool HPhysicsGen::generateQSQ()
{
    //The Q^2 upper limit is motivated by Xbj<=1
    //double qsq_max = 2 * event->getNu() * hepconst::w2proton;
    //Roll the dice!

    bool valid  = false;


    double r,qsq,phaseFactor;
    while(!valid)
    {
        r = myRandom->flat();
        //Generation in order of 1/Q^2 instead of linear
//         qsq = paramMan->getQsqMin() * pow( (paramMan->getQsqMax() / paramMan->getQsqMin()) ,r);


//         double qsqmax = paramMan->getXbjMax() * (2 * event->getStruct()->nu * hepconst::w2prma);
        double qsqmax = paramMan->getXbjMax() * (2 * event->getStruct()->nu * hepconst::w2prma);


        //this implicitly makes sure gamma virt in CMS system has positive energy
        if (paramMan->getXbjMax() > 0.5)
            qsqmax = 0.5 * (2 * event->getStruct()->nu * hepconst::w2prma);

        double ebeam = paramMan->getBeamE();
        double emu = ebeam - event->getNu();
        double pbeam = sqrt(ebeam*ebeam - hepconst::w2mu);
        double pmu = sqrt(emu*emu - hepconst::w2mu);

        double qsqmin = -2 * hepconst::w2mu + 2*(ebeam*emu - pbeam*pmu);

        if (qsqmin < paramMan->getQsqMin())
            qsqmin = paramMan->getQsqMin();
        if (qsqmax > paramMan->getQsqMax())
            qsqmax = paramMan->getQsqMax();

        event->getStruct()->USERVAR.at(7) = qsqmin;
        event->getStruct()->USERVAR.at(8) = qsqmax;

        qsq = qsqmin * pow( (qsqmax / qsqmin) ,r);


        //Calculate phaseFactor
        phaseFactor = (1./qsq) * 1. / log (qsqmax/qsqmin);

        pftqsq_ = 1./ phaseFactor;
	
	double phaseFactorOldBuggy = (1./qsq) * 1. / log (paramMan->getQsqMax()/paramMan->getQsqMin());
	pfoldqsq = (1./phaseFactorOldBuggy)*pft_ * pfnu_;


        if (qsq < (2 * event->getStruct()->nu * hepconst::w2prma))
            valid = true;
        else
            printf("QSQ reroll!\n");



    }


    return event->setQSQ(qsq,phaseFactor);

}




bool HPhysicsGen::generateElastic()
{
    amx2 = hepconst::w2proton;

    //do some checking for elastic/inelastic and diffractive dissociation here - here we just set everything to coherent standard
    //check if we
    //throw the dice for neutron or proton - wether we use it or not, we need the throw for the numerical compatibility with hepgen


    doDD = false;

    if (paramMan->getStruct()->LST.at(19) != 0)
    {
        double sigsd = hepconst::par1sd * (1+hepconst::par2sd/event->getWsq()) * log(0.6+0.1*event->getWsq());
        double rdiss = sigsd/(sigsd + hepconst::sigel);
        double r = myRandom->flat();
        if (r < rdiss)
            doDD = true;
    }



    coherent = false;

    double rand;

    if (!doDD)
    {
        rand = myRandom->flat();
        // if (NO_DIFFRACTIVE_DISSOCIATION)
        coherent = (rand <= paramMan->getStruct()->probc);
    }

    if (paramMan->getStruct()->atomas < 2.0)
        hitNeutron = false;




    //decide proton or neutron target
    if (!coherent && paramMan->getStruct()->atomas > 1.0)
    {
        rand = myRandom->flat();
        hitNeutron = (rand < 0.5);
    }

    
    isNeutronAfterHit = hitNeutron;

    if (hitNeutron)
        event->getStruct()->targetParticle.setParticleType(hepconst::typeNeutron);
    else
        event->getStruct()->targetParticle.setParticleType(hepconst::typeProton);

    event->getStruct()->recoil.setParticleType(event->getStruct()->targetParticle.getParticleType());

    event->getStruct()->ddActive = doDD;

    if (doDD)
    {
        if (sqrt(event->getWsq()) < hepconst::ddmassThreshold)
            return false;

        ddGenerateParticles();

        if (sqrt(event->getWsq()) < sqrt(amx2) + m_meson)
            return false;

	paramMan->getStruct()->LST.at(24) = 1;
	
    }
    else
        paramMan->getStruct()->LST.at(24) = 0;
	




    return true;
}

bool HPhysicsGen::generateSmearing()
{
    double emiss = (amx2 - hepconst::w2proton)/(2. * hepconst::w2prma);

    //now get a gaussian random number
    double emiss_smeared = myRandomGauss->fire()*1.00 + emiss;
    event->setEMiss(emiss_smeared);
    return true;
}



bool HPhysicsGen::generateMesonMass()
{
    //m_meson = 0;
    event->getStruct()->PARL.at(22) = pow(m_meson,2);
    if (sqrt(event->getStruct()->wsq) <= (hepconst::w2prma + m_meson))
    {
        return false;
    }
    else
        return true;
}






HPhysicsGen::~HPhysicsGen()
{
    delete GaussEngine;
    delete myRandomGauss;
}


void HPhysicsGen::setUserVars()
{
    //0 and 1 are reserved for beamfile-parameters

    event->getStruct()->USERVAR.at(4) = event->getStruct()->epp_polarized_longitudinal;
    event->getStruct()->USERVAR.at(10) = paramMan->getStruct()->atomas;
    event->getStruct()->USERVAR.at(11) = paramMan->getStruct()->probc;
    event->getStruct()->USERVAR.at(12) = paramMan->getStruct()->bcoh;
    event->getStruct()->USERVAR.at(13) = paramMan->getStruct()->bin;
    event->getStruct()->USERVAR.at(14) = paramMan->getStruct()->alf;
    //next to be set by special generator


    //continue with standard settings
    event->getStruct()->USERVAR.at(17) = paramMan->getStruct()->B0;
    event->getStruct()->USERVAR.at(18) = paramMan->getStruct()->xbj0;
    event->getStruct()->USERVAR.at(19) = paramMan->getStruct()->alphap;


    //set the phasefactor
    event->getStruct()->USERVAR.at(9) = event->getTotalPhaseFactor();

}



bool HPhysicsGen::generatet()
{
    double slope = paramMan->getSlopeIncoherent();
    double tmin = paramMan->gettMin();
    double tmax = paramMan->gettMax();
    double rand = myRandom->flat();

    //skip 1 random number for numerical compat with hepgen
    //if (paramMan->getIVECM().at(0) == 0)
    //rand = myRandom->flat();

    double pin = sqrt(pow(event->getStruct()->wsq,2)-2*event->getStruct()->wsq*(hepconst::w2proton-event->getStruct()->qsq)+pow((hepconst::w2proton+event->getStruct()->qsq),2))/(2*sqrt(event->getStruct()->wsq));
    double pout = sqrt(pow(event->getStruct()->wsq,2)-2*event->getStruct()->wsq*(amx2+pow(m_meson,2)) +pow((amx2-pow(m_meson,2)),2)) / (2*sqrt(event->getStruct()->wsq));

    double deltat = 2 * pin * pout;
    if (tmax > 2*deltat) {
     //  printf("tlim2 corrected!\n");
        tmax = 2* deltat;
    }



    double tprim = -1./slope*
                   log(exp(-slope*tmin)-(exp(-slope*tmin)-exp(-slope*tmax))*rand);
    double phaseFactor =slope*exp(-slope*tprim)/(exp(-slope*tmin)-exp(-slope*tmax));

    pft_ = 1./phaseFactor;

    return event->setT(tprim, phaseFactor, m_meson, amx2);

}

bool HPhysicsGen::generatePhiGamma()
{
    double rand = myRandom->flat();
    //real calculation of this is done in hevent class for there we have access to the other parameters
    return event->calculateMuonKinematics(rand);
}



bool HPhysicsGen::generateOutgoingParticle()
{
    //calculate the outgoing particles parameters
    double pin = sqrt(pow(event->getWsq(),2)-2*event->getWsq()*(hepconst::w2proton-event->getQsq())+pow((hepconst::w2proton+event->getQsq()),2))/(2*sqrt(event->getWsq()));
    double pout = sqrt(pow(event->getWsq(),2)-2*event->getWsq()*(amx2+pow(m_meson,2)) +pow((amx2-pow(m_meson,2)),2)) / (2*sqrt(event->getWsq()));

    double e_out = sqrt(pow(pout,2) + pow(m_meson,2));

    double egammav = sqrt(event->getWsq()) - sqrt(pow(pin,2)+hepconst::w2proton);

    double cos_out = (event->getT() + event->getQsq() - pow(m_meson,2) + 2 * egammav * e_out)/(2 * pin * pout);


    if (cos_out > 1 )
        cos_out = 1;

    double r = myRandom->flat();
    phi_out = r * 2 * M_PI;
    theta_out = acos(cos_out);


    HLorentzVector outvec = event->getOutPart1().getVector();
    outvec.setLVectorAngular(pout, theta_out, phi_out, e_out);
    event->getOutPart1().setVector(outvec);

    return true;
}






//this method combines genmx and genmul from nucdiss.h
void HPhysicsGen::ddGenerateMultiplicities()
{

    //first we generate the missing mass square of the diffractive dissociation
    double w = sqrt(event->getWsq());
    double massMax = pow(hepconst::w2prma,2) + 0.15 * event->getWsq();
    double prob;
    bool mass_ok = false;
    double mx2;




    while (!mass_ok)
    {
        double r = myRandom->flat();
        mx2 = (massMax - hepconst::ddmassMinSq)*r + hepconst::ddmassMinSq;
        if (sqrt(mx2) > w)
            continue;

        if (mx2 < hepconst::ddmassPeak)
            prob = (hepconst::ddmassPeak/mx2) * ((mx2 - hepconst::ddmassThreshold)/(hepconst::ddmassPeak - hepconst::ddmassThreshold));
        else
            prob = (hepconst::ddmassPeak/mx2);
        r = myRandom->flat();
        if (r < prob)
            mass_ok = true;
    }


    event->getStruct()->ddMX = sqrt(mx2);
    amx2 = mx2;




    mx = sqrt(mx2);
    double amean = 2. * sqrt(mx-hepconst::w2prma);

    if (hitNeutron)
        amean -= 1;



    double sigma = amean / 2.; //for some reason this works out in the data - no real physical explanation for this though
    double amean_n = (amean-0.5) /2.;
    if (hitNeutron)
        amean_n += 1;
    double sigma_n = amean_n /2.;


    bool all_ok = false;



    //set all the particle ids to zero for the reason of seeing what is beeing used and what not
    for (int i =0; i < 30; i++)
        ddParticles[i].setParticleType(0);

    while (!all_ok)
    {
        bool charged_ok = false;

        while(!charged_ok)
        {
            //throw dice for charged particles
            //assuming hydrogen target

            double r = myRandomGauss->fire();
            //   r = -0.154492974;
            double multiplicity = r*sigma + amean;
            if (hitNeutron == false)
            {
                num_charged = (int)((multiplicity/2.))*2. +1;
                charged_ok = (num_charged > 0);
            }
            else
            {
                num_charged = (int)(((multiplicity+1)/2.))*2;
                charged_ok = (num_charged >= 0);
            }

        }


        bool neutral_ok = false;





        while(!neutral_ok)
        {
            //throw dice for charged particles
            //assuming hydrogen target
            double r = myRandomGauss->fire();
            //    r=0.26774624;
            //double r=0.26774624;
            double multiplicity = r*sigma_n + amean_n;

            num_neutral = (int) (multiplicity+0.5);
            // int num_neutral = ((int)(multiplicity+1)/2)*2;
            neutral_ok = (num_neutral >= 0);
// 	    cout << endl << "neutral " << num_neutral << " " << amean_n << " " << sigma_n << " " << multiplicity << " " << hitNeutron << endl;

        }
        //check if particles exist and if energy conservation is fulfilled
        num_particles = num_charged + num_neutral;
        if (num_particles < 2)
        {
            num_particles = 2;

        }
        if( hepconst::w2neutron + (num_particles - 1)*hepconst::w2pic < mx)
            all_ok = true;

    }

    event->getStruct()->ddNum_Charged = num_charged;
    event->getStruct()->ddNum_Neutral = num_neutral;
    event->getStruct()->ddNum_Particles = num_particles;
    event->getStruct()->amx2 = amx2;







    //use proton mass here for filling the pmax into the histograms -
    //yeye i know its not correct, but i need this hotfix to make the same histograms as andrzej
    double targetMass = hepconst::w2proton;
    double pmax = sqrt(pow(amx2,2)-2*amx2*(targetMass + hepconst::w2pi)+ pow((targetMass - hepconst::w2pi),2))/(2*sqrt(amx2));
    event->getStruct()->ddPMax = pmax;


    event->getStruct()->recoil.setParticleDaughter1(event->getStruct()->listOfParticles.size()+1);
    event->getStruct()->recoil.setParticleDaughter2(event->getStruct()->listOfParticles.size()+num_charged+num_neutral);
}

void HPhysicsGen::ddGenerateParticles()
{
    ddGenerateMultiplicities();


    for (unsigned int i = 2; i < 10; i ++)
        ddBookMan->fill(i);

    if (num_particles <= 2)
    {
        ddGenerateKinematics2Body();
    }
    else
    {
        ddGenerateKinematics();
        ddBalanceMultiParticles();
    }

    ddAddParticleList();
}


void HPhysicsGen::ddGenerateKinematics()
{
    double r = myRandom->flat();


    int motherPart;

    event->getStruct()->recoil.setSharpMass(false);
    event->getStruct()->recoil.setParticleAuxFlag(21);


    if (hitNeutron)
        event->getStruct()->recoil.setParticleType(hepconst::typeNeutronExcited);
    else
        event->getStruct()->recoil.setParticleType(hepconst::typeProtonExcited);


    if (r < 0.5)
        afterDecay = hepconst::Proton;
    else
        afterDecay = hepconst::Neutron;


    if (afterDecay == hepconst::Proton && num_charged <= 0)
        afterDecay = hepconst::Neutron;
    if (afterDecay == hepconst::Neutron && num_neutral <= 0)
        afterDecay = hepconst::Proton;

    int chargedCounter = 0;
    int neutralCounter = 0;


    if (afterDecay == hepconst::Proton)
    {
        ddParticles[0].setParticleType(hepconst::typeProton);
        ddParticles[0].setParticleAuxFlag(1);
        ddParticles[0].setSharpMass(true);
        ddParticles[0].setParticleDaughter1(0);
        ddParticles[0].setParticleDaughter2(0);
        ddParticles[0].setParticleOrigin(event->getStruct()->listOfParticles.size());
        chargedCounter++;
    }
    else
    {
        neutralCounter++;
        ddParticles[num_charged].setParticleType(hepconst::typeNeutron);
        ddParticles[num_charged].setParticleAuxFlag(1);
        ddParticles[num_charged].setSharpMass(true);
        ddParticles[num_charged].setParticleDaughter1(0);
        ddParticles[num_charged].setParticleDaughter2(0);
        ddParticles[num_charged].setParticleOrigin(event->getStruct()->listOfParticles.size());
    }

    bool fillPositive = true;
    //fill the charged ones
    int index = ddGetNextFreeIndex(0);

    if (num_charged > chargedCounter)
    {
        //if we have the proton we start with pi- and
        //if we have a neutron we start with pi+
        //to not to accumulate too much charge in the end
        if (index % 2  == 0)
            fillPositive = true;
        else
            fillPositive = false;

        for (chargedCounter; chargedCounter < num_charged; chargedCounter ++)
        {
            if (fillPositive)
                ddParticles[index].setParticleType(hepconst::typePiPlus);
            else
                ddParticles[index].setParticleType(hepconst::typePiMinus);

            ddParticles[index].setSharpMass(true);
            ddParticles[index].setParticleAuxFlag(1);
            ddParticles[index].setParticleDaughter1(0);
            ddParticles[index].setParticleDaughter2(0);
            ddParticles[index].setParticleOrigin(event->getStruct()->listOfParticles.size());
            index=ddGetNextFreeIndex(index);

            //flip the sign
            fillPositive = !fillPositive;
        }

    }
    //fill neutral ones
    if (num_neutral > neutralCounter)
    {
        // the index should be free so we can just use it one time.
        // as it was checked last time in charged or if no charged exist it was checked anyways
        for (neutralCounter; neutralCounter < num_neutral; neutralCounter ++)
        {
            ddParticles[index].setParticleType(hepconst::typePi0);
            ddParticles[index].setParticleOrigin(event->getStruct()->listOfParticles.size());
            ddParticles[index].setParticleDaughter1(0);
            ddParticles[index].setParticleDaughter2(0);
            ddParticles[index].setSharpMass(true);
            index=ddGetNextFreeIndex(index);
        }
    }
    //now that we actually have made all the particles its time to decide on the kinematics phase space, longitudinal for high amx
    //and linear for low amx2

    if (amx2 >= 4)
    {
        ddGenerateLongitudinalPhaseSpace();

    }
    else
    {
        ddGeneratePhaseSpace();
    }
}





//this generates the kinematic data of many particles according to longitudinal phase space
void HPhysicsGen::ddGenerateLongitudinalPhaseSpace()
{
    double targetMass = hepconst::w2proton;
    if (hitNeutron)
        targetMass = hepconst::w2neutron;
    double pmax = sqrt(pow(amx2,2)-2*amx2*(targetMass + hepconst::w2pi)+ pow((targetMass - hepconst::w2pi),2))/(2*sqrt(amx2));
    double ptmean,ptmax,r,pt,phi,amt,ymax,yslope,y;
    HVector3 outvec;
    for (int i =0; i < num_charged+num_neutral; i ++)
    {

        if (ddParticles[i].getParticleType() == hepconst::typeProton || ddParticles[i].getParticleType() == hepconst::typeNeutron)
            ptmean = hepconst::ptpr;
        else
            ptmean = hepconst::ptpi;
        ptmax = 3*ptmean;

        bool physOK = false;
        while (!physOK)
        {
            r = myRandom->flat();
            pt = r * ptmax;
            double prob = (2./ptmean)*exp(1)*pt*exp(-2.*pt/ptmean);
            r = myRandom->flat();
            if (r < prob)
                physOK = true;
        }

        r = myRandom->flat();
        phi = 2*r*M_PI;
        //straight up put momentum in y,z direction
        outvec = HVector3(0.,pt*cos(phi),pt*sin(phi));


        amt = sqrt(pow(ddParticles[i].getMass(),2)+pow(pt,2));
        ymax = log(sqrt(amx2/hepconst::w2proton));

        if (ymax < 2)
            yslope = 1./ymax;
        else
            yslope = 0.5;

        bool yOkay = false;
        while (!yOkay)
        {
            r = myRandom->flat();
            y = (2 * r - 1) * ymax;

            //make the y-distribution flat!
            double prob;
            if((ymax-abs(y))<2)
            {
                prob=-1*(abs(y)-ymax)*yslope;
            }
            else
                prob=1.;

            yOkay = true;

            /* here we can either choose to use this probability or not to use it */
            /* it would look like this:
             *
             * r = myRandom->flat();
             * if (r > prob)
             *   yOkay = false;

             */


        }
        //generate momentum in x direction and add to vector
        outvec = outvec + HVector3((amt/2)*(exp(y)-exp(-y)),0.,0.);


        //apply to particle
        ddParticles[i].getVector().setVector(outvec);

        //update energy
        ddParticles[i].getVector().setEnergy(sqrt(pow(amt,2)+pow(outvec.X(),2)));
    }
}

//this generates the kinematic data of many particles according to linear phase space
void HPhysicsGen::ddGeneratePhaseSpace()
{
    double targetMass = hepconst::w2proton;
    if (hitNeutron)
        targetMass = hepconst::w2neutron;


    double pmax = sqrt(pow(amx2,2)-2*amx2*(targetMass + hepconst::w2pi)+ pow((targetMass - hepconst::w2pi),2))/(2*sqrt(amx2));
    //  cout << "PMAX " << pmax << endl;
    double r,pout;

    for (int i = 0; i < num_particles; i++)
    {
        bool physOK = false;
        while (!physOK)
        {
            r = myRandom->flat();
            pout = r * pmax;
            double prob = pow(r,2.);
            r = myRandom->flat();
            if (r < prob)
                physOK = true;
        }


        r = myRandom->flat();
        double theta = acos(r*2-1);
        r = myRandom->flat();
        double phi = r*2*M_PI;

        ddParticles[i].getVector().setLVectorAngular(pout,theta,phi,sqrt(pow(ddParticles[i].getMass(),2)+pow(pout,2)));

    }
}

HLorentzVector HPhysicsGen::ddSumKinematics()
{
    HLorentzVector result = HLorentzVector(HVector3(0.,0.,0.),0.);
    for (int i =0; i < num_particles; i++)
        result  = result + ddParticles[i].getVector();
    return result;
}



void HPhysicsGen::ddBalanceMultiParticles()
{
    double usePower = 3.0;


    //start with balancing the momenta
    double sumx,sumy,sumz;
    sumx=0.;
    sumy=0.;
    sumz=0.;
    double amx = sqrt(amx2);

    HLorentzVector summedKin = ddSumKinematics();

    for (int i =0; i < num_particles; i ++)
    {
        sumx += abs(pow(ddParticles[i].getTVector().X(),usePower));
        sumy += abs(pow(ddParticles[i].getTVector().Y(),usePower));
        sumz += abs(pow(ddParticles[i].getTVector().Z(),usePower));
    }

    for (int i =0; i < num_particles; i++)
    {
        HVector3 mom = HVector3(abs(pow(ddParticles[i].getTVector().X(),usePower))* summedKin.getVector().X()/sumx,
                                abs(pow(ddParticles[i].getTVector().Y(),usePower))* summedKin.getVector().Y()/sumy,
                                abs(pow(ddParticles[i].getTVector().Z(),usePower))* summedKin.getVector().Z()/sumz);


        ddParticles[i].getVector() = ddParticles[i].getVector() -
                                     HLorentzVector(mom,0.);

        ddParticles[i].getVector().setEnergy(sqrt(pow(ddParticles[i].getMass(),2)+ddParticles[i].getTVector().dotProduct(ddParticles[i].getTVector())));
    }



    summedKin = ddSumKinematics();



    ddIterations = 0;
    double amt = 0.0;
    for (int i = 0; i < num_particles; i ++)
        amt += sqrt(pow(ddParticles[i].getMass(),2)+pow(ddParticles[i].getTVector().Y(),2)+pow(ddParticles[i].getTVector().Z(),2));

    while (abs(amx - summedKin.getEnergy()) > 0.01)
    {
        ddIterations++;
        summedKin = ddSumKinematics();



        if (amx2 <= 4.0 || amt > amx)
        {
            //use normal phase space
            double correction = 0.0;
            for (int i = 0; i < num_particles; i ++)
            {
                correction += pow(ddParticles[i].getTVector().length(),2)/ sqrt(pow(ddParticles[i].getMass(),2)+pow(ddParticles[i].getTVector().length(),2));

            }
            for (int i = 0; i < num_particles; i ++)
            {
                double deltaP = (amx - summedKin.getEnergy())/(0.5*correction);
                double momentumFactor = sqrt(1+deltaP);
                ddParticles[i].scaleMomentum(momentumFactor);
            }
        }
        else
        {
            //use longitudinal phase space
            double pzpos = 0.;
            double pzneg = 0.;
            for (int i = 0; i < num_particles; i ++)
            {
                if (ddParticles[i].getTVector().X()<0)
                    pzneg += pow(abs(ddParticles[i].getTVector().X()),usePower);
                else
                    pzpos += pow(abs(ddParticles[i].getTVector().X()),usePower);
            }

            for (int i = 0; i < num_particles; i ++)
            {
                double px = ddParticles[i].getTVector().X();
                HVector3 tmp = ddParticles[i].getTVector();
                if (px<0)
                    tmp.setX(px - pow(abs(px),usePower) * (amx-summedKin.getEnergy())/(2.*pzneg));
                else
                    tmp.setX(px + pow(abs(px),usePower) * (amx-summedKin.getEnergy())/(2.*pzpos));
                ddParticles[i].setTVector(tmp);
                ddParticles[i].getVector().setEnergy(sqrt(pow(ddParticles[i].getMass(),2.)+tmp.dotProduct(tmp)));
            }
        }
    }
}




int HPhysicsGen::ddGetNextFreeIndex(int _startIndex=0)
{
    int index = _startIndex;
    while (ddParticles[index].getParticleType() != 0)
    {
        if (index >= 30)
            cout << "ERROR: ddParticles index variable in hphysicsgen overflow! crashing!" << endl;
        index++;
    }

    return index;
}


void HPhysicsGen::ddGenerateKinematics2Body()
{
    //first we choose proton or neutron after decay of proton/lambda
    //this also sets our second particle
    //p*->p+pi0
    //p*->n+pi+ because of conservation of charge
    //this needs to be a bit more complicated for switchable target

    afterDecay = hepconst::Neutron;
    if (hitNeutron)
        afterDecay = hepconst::Proton;



    double r = myRandom->flat();
    if (r < 0.33)
    {
        if(hitNeutron)
            afterDecay = hepconst::Neutron;
        else if(!hitNeutron)
            afterDecay = hepconst::Proton;
    }
    else // r>0.66
    {
        if (hitNeutron)
            afterDecay = hepconst::Proton;
        else if (!hitNeutron)
            afterDecay = hepconst::Neutron;
    }

    if (!hitNeutron && afterDecay == hepconst::Proton)
    {   // p* -> p+ + pi0
        num_charged = 1;
        num_neutral = 1;
        ddParticles[0].setParticleType(hepconst::typeNeutron);
        ddParticles[0].setSharpMass(true);
        ddParticles[1].setParticleType(hepconst::typePi0);
        ddParticles[1].setSharpMass(true);
    }
    else if (!hitNeutron && afterDecay == hepconst::Neutron)
    {   // p* -> n + pi+
        num_charged = 1;
        num_neutral = 1;
        ddParticles[0].setParticleType(hepconst::typePiPlus);
        ddParticles[0].setSharpMass(true);
        ddParticles[1].setParticleType(hepconst::typeNeutron);
        ddParticles[1].setSharpMass(true);
    }
    else if (hitNeutron && afterDecay == hepconst::Neutron)
    {   //n* -> n + pi0
        num_neutral = 2;
        num_charged = 0;
        ddParticles[0].setParticleType(hepconst::typeNeutron);
        ddParticles[0].setSharpMass(true);
        ddParticles[1].setParticleType(hepconst::typePi0);
        ddParticles[1].setSharpMass(true);
    }
    else if (hitNeutron && afterDecay == hepconst::Proton)
    {   //n*->p+ + pi-
        num_charged = 2;
        num_neutral = 0;
        ddParticles[0].setParticleType(hepconst::typePiMinus);
        ddParticles[0].setSharpMass(true);
        ddParticles[1].setParticleType(hepconst::typeProton);
        ddParticles[1].setSharpMass(true);
    }


    if (!hitNeutron)
        event->getStruct()->recoil.setParticleType(hepconst::typeProtonExcited);
    else
        event->getStruct()->recoil.setParticleType(hepconst::typeNeutronExcited);

    event->getStruct()->recoil.setSharpMass(false);
    event->getStruct()->recoil.setParticleAuxFlag(21);


    //set the correct number straight
    num_particles = 2;

    double pout = sqrt(pow(amx2,2) - 2 * amx2 * (pow(ddParticles[0].getMass(),2) + pow(ddParticles[1].getMass(),2))
                       + pow((pow(ddParticles[0].getMass(),2) - pow(ddParticles[1].getMass(),2)),2))/(2*sqrt(amx2));

    r = myRandom->flat();
    double theta = acos(2*r-1);
    r = myRandom->flat();
    double phi = 2*M_PI*r;

    ddParticles[0].getVector().setLVectorAngular(pout,theta,phi,sqrt(pow(ddParticles[0].getMass(),2)+pout*pout));
    ddParticles[1].getVector().setVector(HVector3(0,0,0) -ddParticles[0].getTVector());
    ddParticles[1].getVector().setEnergy(sqrt(pow(ddParticles[1].getMass(),2)+pout*pout));
    ddParticles[0].setParticleAuxFlag(1);
    ddParticles[1].setParticleAuxFlag(1);


    ddParticles[0].setParticleOrigin(event->getStruct()->listOfParticles.size());
    ddParticles[1].setParticleOrigin(event->getStruct()->listOfParticles.size());


    ddParticles[0].setParticleDaughter1(0);
    ddParticles[1].setParticleDaughter1(0);
    ddParticles[0].setParticleDaughter2(0);
    ddParticles[1].setParticleDaughter2(0);



}


void HPhysicsGen::ddDecayParticle(int _index, int _startIndex)
{
    // first we generate the angles
    double r = myRandom->flat();
    double costh;
    //generate theta flat
    costh = 2*r - 1;

    //generate phi now
    r = myRandom->flat();
    double phiGamma1 =r*2*M_PI;

    // then we generate the particle out of the angles
    double thetaGamma = acos (costh);

    HLorentzVector decayGamma1_cms;
    //now we set it in the CMS
    //momentum of the pions
    double gammamom = hepconst::w2pi/2.;
    //generate the lorentz vector
    decayGamma1_cms.setLVectorAngular(gammamom,thetaGamma,phiGamma1,gammamom);
    HVector3 pineg_mom = HVector3(0.,0.,0.)-decayGamma1_cms.getVector();
    HLorentzVector decayGamma2_cms;
    decayGamma2_cms.setVector(pineg_mom);
    decayGamma2_cms.setEnergy(decayGamma1_cms.getEnergy());


    decayGamma1_cms.boost(2*gammamom,ddParticles[_index].getVector());
    decayGamma2_cms.boost(2*gammamom,ddParticles[_index].getVector());

    //now we generate the negative pion by just inverting the 3-vector of momentum
    //HLorentzVector decayPi2_lab = event->goToLabSystem(decayGamma2_cms,event->getOutPart1_Lab().getVector());


    //okay all done - now set all the flags and mark the particles according to lepto standard
    //then add them to the list
    ddParticles[_index].setParticleAuxFlag(11);

    int newIndex_gamma1 = ddGetNextFreeIndex(_index);

    ddParticles[newIndex_gamma1].setParticleType(hepconst::typeGamma);
    ddParticles[newIndex_gamma1].setParticleAuxFlag(1);
    ddParticles[newIndex_gamma1].setParticleOrigin(baseNumberOfParticles+_index+1);
    ddParticles[newIndex_gamma1].setParticleType(hepconst::typeGamma);
    ddParticles[newIndex_gamma1].setVector(decayGamma1_cms);
    ddParticles[newIndex_gamma1].setParticleDaughter1(0);
    ddParticles[newIndex_gamma1].setParticleDaughter2(0);

    int newIndex_gamma2 = ddGetNextFreeIndex(_index);
    ddParticles[newIndex_gamma2].setParticleAuxFlag(1);
    ddParticles[newIndex_gamma2].setParticleOrigin(baseNumberOfParticles+_index+1);
    ddParticles[newIndex_gamma2].setParticleType(hepconst::typeGamma);
    ddParticles[newIndex_gamma2].setVector(decayGamma2_cms);
    ddParticles[newIndex_gamma2].setParticleDaughter1(0);
    ddParticles[newIndex_gamma2].setParticleDaughter2(0);

    ddParticles[_index].setParticleDaughter1(baseNumberOfParticles+newIndex_gamma1+1);
    ddParticles[_index].setParticleDaughter2(baseNumberOfParticles+newIndex_gamma2+1);

    num_neutral += 2;
    num_particles += 2;

}



void HPhysicsGen::ddAddParticleList()
{
    //get our offset for any further decay particles to do numbers
    baseNumberOfParticles = event->getStruct()->listOfParticles.size();

    int num_particles_old = num_particles;

    //now do the pi0 decay
    for (int i=0; i < num_particles_old; i ++)
        if (ddParticles[i].getParticleType() == hepconst::typePi0)
            ddDecayParticle(i,num_particles_old);


    //then add everything to the list of generated particles
    for (int i =0; i < num_particles; i++)
        event->getStruct()->listOfParticles.push_back(&ddParticles[i]);
}






void HPhysicsGen::ddGotoLab()
{
    //transform all the particles to the lab system
    for (int i = 0; i < num_particles; i ++)
        ddParticles[i].getVector().boost(sqrt(amx2),event->getRecoil().getVector());
}





