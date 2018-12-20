#include "hevent.h"

HEvent::HEvent(HParamManager* _paramMan)
{
    getStruct()->dummyW =1;
    paramMan = _paramMan;
    // set the incBeamParticle parameters
    reset();
}


void HEvent::copyFrom(HEvent* _myEvent)
{
    dataStruct = *(_myEvent->getStruct());
}





void HEvent::reset()
{
    //cout << "resetting event" << endl;
    getStruct()->nu = 0;
    getStruct()->qsq = 0;
    getStruct()->t =0;
    getStruct()->y=0;
    getStruct()->xbj=0;

    getStruct()->ddMX =0;
    getStruct()->ddNum_Charged=0;
    getStruct()->ddNum_Neutral=0;
    getStruct()->ddNum_Particles=0;
    getStruct()->ddActive = false;



    getStruct()->incBeamParticle.setParticleType(paramMan->getBeamPart());
    getStruct()->scatBeamParticle.setParticleType(paramMan->getBeamPart());

    cout.precision(8);

    //set beam directly on z-axis
    HLorentzVector tmp = HLorentzVector((HVector3(sqrt(pow(paramMan->getBeamE(),2)-hepconst::w2mu),0.,0.)),paramMan->getBeamE());
    getStruct()->incBeamParticle.setVector(tmp);

    //this is already in lab system since we dont need it right now but later for the dumps
    getStruct()->targetParticle.setVector(HLorentzVector(HVector3(0.,0.,0.),hepconst::w2prma));


    //now we need to setup the parl-vector
    getStruct()->PARL = vector<double> (paramMan->getStruct()->PARL);

    //then we add some new elements to it for the event-specific configuration
    for (int i = 0; i < 11; i++)
        getStruct()->PARL.push_back(0.0);

    static const double uservartmp[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                                        0.,0.,0.,0.,0.,0.,0.,0.,0.,0.
                                       };

    getStruct()->USERVAR = vector<double> (uservartmp, uservartmp + sizeof(uservartmp) / sizeof(uservartmp[0]));


}

bool HEvent::setT(double _t, double _phaseFactor, double _mesonMass, double _amx2)
{
    getStruct()->t = _t;
    getStruct()->PFt = _phaseFactor;

    getStruct()->m_meson = _mesonMass;
    getStruct()->amx2 = _amx2;

    double pin = sqrt(pow(getStruct()->wsq,2)-2*getStruct()->wsq*(hepconst::w2proton-getStruct()->qsq)+pow((hepconst::w2proton+getStruct()->qsq),2))/(2*sqrt(getStruct()->wsq));
    double pout = sqrt(pow(getStruct()->wsq,2)-2*getStruct()->wsq*(_amx2+pow(_mesonMass,2)) +pow((_amx2-pow(_mesonMass,2)),2)) / (2*sqrt(getStruct()->wsq));
    double egammav = sqrt(getStruct()->wsq) - sqrt(pow(pin,2)+hepconst::w2proton);
    
    
    double e_out = sqrt(pow(pout,2) + pow(_mesonMass,2));

    //check for sanity of the t by checking the limits of t
    double tlimit = - getStruct()->qsq +  pow(_mesonMass,2) - 2 * egammav * e_out;

    double deltat = 2 * pin * pout;

    double tlimit1 = tlimit - deltat;
    double tlimit2 = tlimit + deltat;


    double tzero = -1. * tlimit2;

    getStruct()->t +=  tzero;
    getStruct()->t *= -1;


    //set the PARLs here
    getStruct()->PARL.at(23) = _amx2;
    getStruct()->PARL.at(24) = _t; //'tprim'
    getStruct()->PARL.at(25) = getStruct()->t;
    getStruct()->PARL.at(26) = tzero;


    //check if everything is ok
     if (egammav < 0 || getStruct()->t < tlimit1 || getStruct()->t > tlimit2)
     {
         printf("TFAIL!\n");
         return false;
     }

    getStruct()->tprim = _t;
    return true;
}




bool HEvent::setNu(double _nu, double _phaseFactor)
{
    //set the kinematics
    getStruct()->nu = _nu;
    getStruct()->PFnu = _phaseFactor;

    getStruct()->y = getStruct()->nu / paramMan->getBeamE();

    if (getStruct()->y < paramMan->getyMin() || getStruct()->y > paramMan->getyMax())
        return false;

    //take the nu from the incoming muon and give it to the scattered one as per definition of mu
    getStruct()->scatBeamParticle.getVector().setEnergy(paramMan->getBeamE()-getStruct()->nu);





    getStruct()->gammaVirt.getVector().setEnergy(getStruct()->nu);

    double ebeampart = paramMan->getBeamE() - getStruct()->nu;


    getStruct()->gammaVirt.setParticleType(hepconst::typeGamma);
    getStruct()->scatBeamParticle.setParticleType(getStruct()->incBeamParticle.getParticleType());


    // check if nu is sane by checking if it exceeds mass of muon, if not sane just do it again!
    if (ebeampart < sqrt(hepconst::w2mu))
        return false;
    else
        return true;



}


bool HEvent::setQSQ(double _qsq, double _phaseFactor)
{
    getStruct()->qsq = _qsq;
    getStruct()->PFqsq = _phaseFactor;

    //these are prelim because we didnt to the kinematics yet
    double ebeam = paramMan->getBeamE();
    double emu = ebeam - getNu();
    double pbeam = sqrt(ebeam*ebeam - hepconst::w2mu);
    double pmu = sqrt(emu*emu - hepconst::w2mu);

    getStruct()->xbj = getStruct()->qsq / (2* hepconst::w2prma * getStruct()->nu);
    getStruct()->wsq = hepconst::w2proton - getStruct()->qsq + 2 * hepconst::w2prma * getStruct()->nu;

    //check if qsq is sane if not, we just start over
    double qsq_min = -2 * hepconst::w2mu + 2*(ebeam*emu - pbeam*pmu);

//     cout << "WSQ " << getStruct()->wsq << endl;
//     cout << "QSQ_MIN, QSQ_MAX, QSQ " << qsq_min << " " << (2 * getStruct()->nu * hepconst::w2prma) << " " << getStruct()->qsq << endl;

    if (getStruct()->qsq <= qsq_min || getStruct()->wsq <= hepconst::w2proton || getStruct()->qsq > (2 * getStruct()->nu * hepconst::w2prma)){
        printf("QSQ skipped! %e %e \n",_qsq,getStruct()->wsq);
        return false;
    }
    else
        return true;
}








bool HEvent::setOutgoingParticle(double _phi, double _theta, double _pout)
{
    //make the lorentz vector from the angles
    HLorentzVector tmp;
    tmp.setLVectorAngular(_pout, _theta, _phi, sqrt(pow(_pout,2) + pow(getStruct()->m_meson,2) ));

    getStruct()->outPart1.setVector(tmp);
    //setting particle identity is performed by each specialized generator
    return true;

}


void HEvent::calcPT()
{
    HLorentzVector outPart1Lab = goToLabSystem(getStruct()->outPart1.getVector(),getStruct()->CMS);
    //calculate z now cause we already did the transformation anyways
    getStruct()->z = outPart1Lab.getEnergy()/getStruct()->nu;
    getStruct()->outPart1_Lab.setVector(outPart1Lab);

    double pout2 = outPart1Lab.getVector().dotProduct(outPart1Lab.getVector());

    HVector3 pgammavirtn = getStruct()->gammaVirt.getTVector();
    pgammavirtn.normalize(1.);
    double scprod = outPart1Lab.getVector().dotProduct(pgammavirtn);
    getStruct()->ptprim = pout2 - pow(scprod,2);



    //do the recoil now too
    getStruct()->recoil.setVector(HLorentzVector (getStruct()->gammaVirt.getTVector() - outPart1Lab.getVector(),    //the threeVector
                                  getStruct()->nu + hepconst::w2prma - outPart1Lab.getEnergy()));  // the energy



}





HLorentzVector HEvent::goToLabSystem(HLorentzVector& _particle, HLorentzVector& _rot)
{
// cout << endl <<"debug for going to lab system" << endl;
    //generate the rotation matrix
//  cout << "_particle "; _particle.print();
//  cout << "_rot "; _rot.print();
//  cout << "_beam "; getBeam().getTVector().print();


    HVector3 px = _rot.getVector();
    px.normalize(1.);
    HVector3 pz = px.crossProduct(getBeam().getTVector());
    pz.normalize(1.);
    HVector3 py = pz.crossProduct(px);

    //this matrix is actually in the cern-notation of writing a matrix G(i,j) as G(j, i)..
    //this dates back to fortran libs beeing able to multiply faster if used this way

    std::vector< std::vector<double> > rotMatrix;
    rotMatrix.push_back(px.toStdVector());
    rotMatrix.push_back(py.toStdVector());
    rotMatrix.push_back(pz.toStdVector());

    //now we multiply this to the outgoing particle vector.
    vector<double> newVec;
    for (int j = 0; j<3; j++)
    {
        double sum = 0;
        for (int i = 0; i<3; i++)
        {
            sum += rotMatrix.at(i).at(j) * _particle.getVector().toStdVector().at(i);
        }
        newVec.push_back(sum);
    }

    HVector3 newOutPart(newVec);

    HLorentzVector labVec(newOutPart,_particle.getEnergy());


//  cout << " Before Boost: ";
// labVec.print();

    //we boost with wsq along the CMS vector now
    //labVec.boost(sqrt(getStruct()->wsq),getStruct()->cms);
    labVec.boost(sqrt(_rot.getQuare()),_rot);


    //cout << sqrt(getStruct()->wsq) << endl;

// cout << " After Boost: ";
    //labVec.print();





    return labVec;

}



void HEvent::rotXToZ()
{   /*
       getBeam().RotXToZ();
       getScat().RotXToZ();
       getGammaVirt().RotXToZ();
       getCMS().rotXToZ();
       getOutPart1().RotXToZ();
       getOutPart2().RotXToZ();
       getOutPart3().RotXToZ();

       getOutPart1_Lab().RotXToZ();
       getOutPart2_Lab().RotXToZ();
       getOutPart3_Lab().RotXToZ();

       getRecoil().RotXToZ();*/
    for (unsigned int i =0; i < getStruct()->listOfParticles.size(); i++)
        getStruct()->listOfParticles.at(i)->RotXToZ();
//     cout << "called rot" << endl;

}




void HEvent::printDebug()
{
    cout << " ---- Debug Info ---- " << endl;
    cout << "nu: " << getNu() << endl;
    cout << "qsq: " << getQsq() << endl;
    cout << "wsq: " << getWsq() << endl;
    cout << "Xbj: " << getXbj() << endl;
    cout << "Emiss: " << getEmiss() << endl;
    cout << "z: " << getStruct()->z << endl;
    cout << "y: " << getY() << endl;
    cout << "ptprim: " << getStruct()->ptprim << endl;
    cout << "t' " << getStruct()->tprim << endl;
    cout << " ---- final vectors ----" << endl;
    cout  <<"!mu:\t";
    getBeam().getVector().print();
    cout << "!p:\t";
    getStruct()->targetParticle.getVector().print();
    cout << "!gv:\t";
    getStruct()->gammaVirt.getVector().print();
    cout << "!mu':\t";
    getScat().getVector().print();
    cout << "!p':\t";
    getRecoil().getVector().print();
    cout << "!op1:\t";
    getStruct()->outPart1_Lab.getVector().print();
    cout << "!op2:\t";
    getStruct()->outPart2_Lab.getVector().print();
    cout << "!op3:\t";
    getStruct()->outPart3_Lab.getVector().print();
    cout << "!CMS:\t";
    getCMS().print();
    cout << " --------- Phase Factors -----------" << endl;
    cout << "PFnu: " << getStruct()->PFnu << endl;
    cout << "PFqsq: " << getStruct()->PFqsq << endl;
    cout << "PFtprim: " << getStruct()->PFt << endl;
    cout << "---------- Weight ----------" << endl;
    cout << "GeneratorWeight: " << getStruct()->USERVAR.at(1) << " " << getStruct()->USERVAR.at(2) << endl;
    
    cout << "------------- USERVARS -------" << endl;
    for (int i = 0; i < 4; i++){
      printf("%e\t%e\t%e\t%e\t%e\n",dataStruct.USERVAR.at(5*i+0),dataStruct.USERVAR.at(5*i+1),dataStruct.USERVAR.at(5*i+2),dataStruct.USERVAR.at(5*i+3),dataStruct.USERVAR.at(5*i+4));
    }
//


    cout << "---------------------" << endl;
}






bool HEvent::calculateMuonKinematics(double _flatrandphi)
{


    //set properties for gammaVirt
    getStruct()->gammaVirt.getVector().setEnergy(getStruct()->nu);
    //set the momentum
    double pgamma = sqrt((pow(getStruct()->nu,2)+getStruct()->qsq));
    //calculate the angles of the vectors
    double pbe = getStruct()->incBeamParticle.getVector().getVector().length();
    double pmu = sqrt(pow(getStruct()->scatBeamParticle.getVector().getEnergy(),2)-hepconst::w2mu);

    double qsq_mn = hepconst::w2mu*getY()*getY()/(1-getY());


    if ((getStruct()->qsq-qsq_mn)/(4. * pbe * pmu) > 1.0)
        return false;

    double costh = sqrt(1. - (getStruct()->qsq-qsq_mn)/(4. * pbe * pmu));



    double themu = acos(costh) * 2.;
    double singa = sin (themu) * pmu/pgamma;
    double cosga = sqrt(1.-pow(singa,2));
    double thega = acos(cosga);

    double y = getY();




    //Roll the dice for Phi-distribution
    double phiga = _flatrandphi * 2 * M_PI;
    double qsq = getStruct()->qsq;

    //epsilon, delta and flux
    getStruct()->gam2 = getStruct()->qsq / pow(getStruct()->nu,2);
    double gam2 = getStruct()->gam2;
    getStruct()->epsilon = (1-y- y*y  *gam2/4)/ (1-y+y*y/2+y*y*gam2/4-( hepconst::w2mu /qsq)*y*y*(1+gam2));
    getStruct()->delta = 2 * hepconst::w2mu * (1.-getStruct()->epsilon)/getStruct()->qsq;
    getStruct()->flux = hepconst::alfa * (getStruct()->nu - getStruct()->qsq / (2 * hepconst::w2prma))
                        / (2 * M_PI * getStruct()->qsq * pow(getStruct()->incBeamParticle.getVector().getEnergy(),2) * pow(getStruct()->y,2))
                        * ( pow(getStruct()->y,2) * (1 - 2 * hepconst::w2mu / getStruct()->qsq)
                            + (1 - getStruct()->y - 0.25 * getStruct()->gam2 * pow(getStruct()->y,2)) * 2 / (1+getStruct()->gam2));

    //set gammaVirt-vector from angles
    getStruct()->gammaVirt.getVector().setLVectorAngular(pgamma,thega,phiga,getStruct()->nu);

    //set the scattered muon-4vector
    getStruct()->scatBeamParticle.setVector(HLorentzVector( getStruct()->incBeamParticle.getVector() - getStruct()->gammaVirt.getVector()) );

    //Calculate boost-4vector to CMS
    getStruct()->CMS = HLorentzVector(getStruct()->gammaVirt.getVector().getVector(),getStruct()->gammaVirt.getVector().getEnergy()+hepconst::w2prma);
    getStruct()->s = 2 * hepconst::w2prma * paramMan->getBeamE() + hepconst::w2proton + hepconst::w2mu;
    //set the PARLS
    getStruct()->PARL.at(20)=getBeam().getEnergy();
    getStruct()->PARL.at(21)=themu;
    return true;
}


void HEvent::rotateEventToBeam(HVector3 _beam)
{
    //when we call this, everything is already at the z axis, so first rotate the beam too
    _beam.rotXToZ();
    //then we calculate the rotation angles
    double dPTrans = sqrt(_beam.X()*_beam.X()+_beam.Y()*_beam.Y());
    double dPTot = _beam.length();
    double thetaBeam =  acos(_beam.Z()/dPTot);
    double phiBeam;
    if (dPTrans > 0.0)
    {
        phiBeam = asin(_beam.Y()/dPTrans);
        if(_beam.X()<0)
            phiBeam = M_PI - phiBeam;
    }
    else
        phiBeam = 0;


    double rotmat[3][3];
    rotmat[0][0] = cos(thetaBeam)*cos(phiBeam);
    rotmat[0][1] = -sin(phiBeam);
    rotmat[0][2] = sin(thetaBeam)*cos(phiBeam);
    rotmat[1][0] = cos(thetaBeam)*sin(phiBeam);
    rotmat[1][1] = cos(phiBeam);
    rotmat[1][2] = sin(thetaBeam)*sin(phiBeam);
    rotmat[2][0] = -sin(thetaBeam);
    rotmat[2][1] = 0.;
    rotmat[2][2] = cos(thetaBeam);

    for (unsigned int i = 0; i < getStruct()->listOfParticles.size(); i++)
    {
        const HVector3& tmp = getStruct()->listOfParticles.at(i)->getTVector();
        HVector3 newVec;
        newVec.setX(rotmat[0][0]*tmp.X()+rotmat[0][1]*tmp.Y()+rotmat[0][2]*tmp.Z());
        newVec.setY(rotmat[1][0]*tmp.X()+rotmat[1][1]*tmp.Y()+rotmat[1][2]*tmp.Z());
        newVec.setZ(rotmat[2][0]*tmp.X()+rotmat[2][1]*tmp.Y()+rotmat[2][2]*tmp.Z());



        getStruct()->listOfParticles.at(i)->getVector().setVector(newVec);
    }

}


