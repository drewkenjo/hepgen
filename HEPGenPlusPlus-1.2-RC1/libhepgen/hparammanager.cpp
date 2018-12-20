/*
 * hparammanager.cpp
 *
 *  Created on: Jan 16, 2013
 *      Author: Christopher Regali
 *
 *
 */

#include "hparammanager.h"

HParamManager::HParamManager(HPionicData* _data)
{
    resetData();
    pionicData = _data;
}

HParamManager::HParamManager(string _filename)
{
    resetData();
    readFile(_filename);
//     if (paramParser.getKeyContents("PI0_FILE").size() > 1)
//     {
//       string piFileName = paramParser.getKeyContents("PI0_FILE").at(1);
//       pionicData = new HPionicData(piFileName);
//     }
//     else
    pionicData = NULL;
}


HParamManager::HParamManager(HParams* _paramStruct)
{
    resetData();
    paramStruct = (*_paramStruct);
}




void HParamManager::readFile ( string _filename )
{
    //first we read-in the whole file with HCardParser class
    redFile = _filename;
    paramParser.reset();
    paramParser.parseFile( _filename );

    //now we check for the keys and set the changed parameters accordingly

    //NGEV the number of events we await 1 int for that key
    vector<string> result = paramParser.getKeyContents("NGEV");
    if (result.size() != 2 )
        cout << "WARNING: Unexpected count of values found for key NGEV "<< result.size() << endl;
    paramStruct.NEVENT = hephelp::StrToInt( result.at(1) );

    //BEAM we await 5 values for this key - some ints, some doubles

    result = paramParser.getKeyContents("BEAM");
    if (result.size() != 6 )
        cout << "WARNING: Unexpected count of values found for key BEAM "<< result.size() << endl;
    paramStruct.ELEPT = hephelp::StrToDouble(result.at(1));
    paramStruct.EHADR = hephelp::StrToDouble(result.at(2));
    paramStruct.IPART = hephelp::StrToInt(result.at(3));
    paramStruct.PARL.at(0) = hephelp::StrToDouble(result.at(4));
    paramStruct.PARL.at(1) = hephelp::StrToDouble(result.at(5));

    //shall we read the beamfile, this is not implemented yet in hepgen_c so do not use :)
    result = paramParser.getKeyContents("BMRD");
    if (result.size() != 2 )
        cout << "WARNING: Unexpected count of values found for key BMRD "<< result.size() << endl;
    if (result.at(1) == "UNSET")
      result.at(1) = "0";
    paramStruct.I_BEAMREAD = hephelp::StrToInt(result.at(1));

    //Use this physics programme
    result = paramParser.getKeyContents("PROC");
    if (result.size() != 2 )
        cout << "WARNING: Unexpected count of values found for key PROC "<< result.size() << endl;
    paramStruct.IPROC = hephelp::StrToInt(result.at(1));

    //Lepto soft Cuts
    result = paramParser.getKeyContents("CUTL");
    if (result.size() != 15 )
        cout << "WARNING: Unexpected count of values found for key CUTL "<< result.size() << endl;
    for (unsigned int i =1; i < result.size(); i++)
        paramStruct.CUT.at(i-1)=hephelp::StrToDouble(result.at(i));

    //Primary limits for generations
    result = paramParser.getKeyContents("TLIM");
    if (result.size() != 3 )
        cout << "WARNING: Unexpected count of values found for key TLIM "<< result.size() << endl;
    double tlimtmp1 = hephelp::StrToDouble(result.at(1));
    double tlimtmp2 = hephelp::StrToDouble(result.at(2));

    //Acceptance for the spectrometer
    result = paramParser.getKeyContents("MACC");
    if (result.size() != 2 )
        cout << "WARNING: Unexpected count of values found for key MACC "<< result.size() << endl;
    paramStruct.THMAX = hephelp::StrToDouble(result.at(1));

    //Vector Meson and its decays
    result = paramParser.getKeyContents("VMES");
    if (result.size() != 3 )
        cout << "WARNING: Unexpected count of values found for key VMES "<< result.size() << endl;
    paramStruct.IVECM.at(0) = hephelp::StrToInt(result.at(1));
    paramStruct.IVECM.at(1) = hephelp::StrToInt(result.at(2));


    //diffractive nucleon  lepto switch set?
    result = paramParser.getKeyContents("DIFF");
    if (result.size() != 2 )
        cout << "WARNING: Unexpected count of values found for key DIFF "<< result.size() << endl;
    if (result.at(1) == "UNSET")
      result.at(1) = "0";
    paramStruct.LST.at(19) = hephelp::StrToInt(result.at(1));

    //read target parameters
    result = paramParser.getKeyContents("TPAR");
    if (result.size() != 5 )
        cout << "WARNING: Unexpected count of values found for key TPAR "<< result.size() << endl;
    paramStruct.atomas = hephelp::StrToDouble(result.at(1));
    paramStruct.probc = hephelp::StrToDouble(result.at(2));
    paramStruct.bcoh = hephelp::StrToDouble(result.at(3));
    paramStruct.bin = hephelp::StrToDouble(result.at(4));


    //alpha dependance
    result = paramParser.getKeyContents("ALFA");
    if (result.size() != 2 )
        cout << "WARNING: Unexpected count of values found for key ALFA "<< result.size() << endl;
    paramStruct.alf = hephelp::StrToDouble(result.at(1));


    //lepton charge and polarisation
    result = paramParser.getKeyContents("BPAR");
    if (result.size() != 3 )
        cout << "WARNING: Unexpected count of values found for key BPAR "<< result.size() << endl;
    paramStruct.clept = hephelp::StrToDouble(result.at(1));
    paramStruct.slept = hephelp::StrToDouble(result.at(2));

    //x,t correlation for dvcs
    result = paramParser.getKeyContents("REGG");
    if (result.size() != 4 )
        cout << "WARNING: Unexpected count of values found for key REGG "<< result.size() << endl;
    paramStruct.B0 = hephelp::StrToDouble(result.at(1));
    paramStruct.xbj0 = hephelp::StrToDouble(result.at(2));
    paramStruct.alphap = hephelp::StrToDouble(result.at(3));


    result = paramParser.getKeyContents("OUTFILE");
    paramStruct.outFile = result.at(1);



    //now we just have to set the vector PARL correctly 0 and 1 are set at the target parameters so this is already done
    paramStruct.PARL.push_back(paramStruct.ELEPT); //2
    paramStruct.PARL.push_back(paramStruct.EHADR);
    paramStruct.PARL.push_back(paramStruct.IPART); //4
    paramStruct.PARL.push_back(paramStruct.IPROC);

    // to be filled by kfvm, kfdp1, kfdp2 codes for vector-meson, and the 2 decay particles
    paramStruct.PARL.push_back(0.0); //6
    paramStruct.PARL.push_back(0.0);
    paramStruct.PARL.push_back(0.0); //8

    paramStruct.PARL.push_back(paramStruct.THMAX);
    paramStruct.PARL.push_back(paramStruct.atomas); //10
    paramStruct.PARL.push_back(paramStruct.probc);
    paramStruct.PARL.push_back(paramStruct.bcoh); //12
    paramStruct.PARL.push_back(paramStruct.bin);
    paramStruct.PARL.push_back(paramStruct.alf);  //14
    paramStruct.PARL.push_back(tlimtmp1);
    paramStruct.PARL.push_back(tlimtmp2); //16
    paramStruct.PARL.push_back(paramStruct.clept);
    paramStruct.PARL.push_back(paramStruct.slept); //18
    //done constructing the parameters-vector


    //read the aux tags into the aux vector
    paramStruct.aux1 = paramParser.getKeyContents("AUX1");
    paramStruct.aux2 = paramParser.getKeyContents("AUX2");
    paramStruct.aux3 = paramParser.getKeyContents("AUX3");
    paramStruct.aux4 = paramParser.getKeyContents("AUX4");
    paramStruct.aux5 = paramParser.getKeyContents("AUX5");





}



void HParamManager::resetData()
{
    HParamManager::resetParamStruct(&paramStruct);
}




void HParamManager::resetParamStruct(HParams* toReset)
{
    toReset->NEVENT = 10;
    toReset->ELEPT = 160.0;
    toReset->EHADR = 0.9382723128;
    toReset->IPART = -13; //mu+
    toReset->I_BEAMREAD = 0;
    toReset->IPROC = 5; //exclusive VM production

    toReset->PARL.clear();
    toReset->PARL.push_back(1.0); //A of TARGET (1.0 for proton)
    toReset->PARL.push_back(1.0); //Z of TARGET (1.0 for proton)

    toReset->IVECM.clear();
    toReset->IVECM.push_back(0); //these are DVCS defaults no mesons
    toReset->IVECM.push_back(0); //will be read again from ucards datafile so its ok to just set basically anything here



    //now this part actually is really ugly - a big failing of the standard :/
    //but to do is fast this is a pretty ok method. other fast initializations include libBOOST.
    static const int leptoswitchtmp[] = {	4  ,1  ,5  ,1  ,3  ,
                                            1  ,1  ,0  ,5  ,1  ,
                                            0  ,4  ,5  ,5  ,100,
                                            100,1  ,2  ,-1,0,
                                            0,0,1,0,1,
                                            0,0,0,0,0,
                                            0,0,0,2,1,
                                            0,0,0,0,0
                                        };

    toReset->LST = vector<int> (leptoswitchtmp, leptoswitchtmp + sizeof(leptoswitchtmp) / sizeof(leptoswitchtmp[0]));



    // set the soft kinematic cuts
    static const double cutstmp[] = {0.0001, 1.0	   ,		// x - range (min, max)
                                     0.0	 , 1.0     ,		// y - range
                                     0.5   , 10000.0 ,		// Q2- range
                                     0.0   , 10000.0 ,		// W2- range
                                     0.0   , 10000.0 ,		// nu- range
                                     0.5   , 10000.0 ,		// E'- range
                                     0.0   , 3.145927		// min lepton scattering angle, pi
                                    };
    toReset->CUT = vector<double> (cutstmp, cutstmp + sizeof(cutstmp) / sizeof(cutstmp[0]));


    //assorted other parameters

    toReset->THMAX = 0.050; //scattered muon acceptance
    toReset->ivecm = 2;     //vector meson flag
    toReset->kcdec = 211;   //absolute value of decay particle id
    toReset->alf = 0.75;    //value of parameter for A dependance
    toReset->atomas = 12;   //average atomic mass for target
    toReset->probc = 0.7;   //fraction of coherent events
    toReset->bcoh = 52.5;   //slope of the nuclear formfactor
    toReset->bin = 4.5;     //slope for the production on a nucleon
    toReset->rdiss = 0.2;   //ratio of incoherent events with nucleon dissociation
    //to incoherent elastic events
    toReset->clept = -1;    //lepton beam charge - default electron
    toReset->slept = 0;     //lepton beam polarisation - default no polarisation
    toReset->B0 = 4.94116;  //parameter for x,t correlation	|
    toReset->xbj0 = 0.042;  //parameter for x,t correlation	|=> for DVCS only
    toReset->alphap = 0.8;  //parameter for x,t correlation	|

    //schildknecht et al parametrization parameters for lepton polarizsation
    toReset->aksi2 = 1.;
    toReset->amt2_rho = 0.62*pow(hepconst::mRho0,2);
    toReset->amt2_phi = 0.40*pow(hepconst::mPhi,2);

    toReset->aml2_rho = 1.5*toReset->amt2_rho;
    toReset->aml2_phi = 1.5*toReset->amt2_phi;
    toReset->slopeR = 0.07;


    toReset->targetType = hepconst::Proton;
}


