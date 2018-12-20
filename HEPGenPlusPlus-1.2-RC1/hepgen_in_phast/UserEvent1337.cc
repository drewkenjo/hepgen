#include <iostream>
#include <cmath>
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "Phast.h"
#include "PaSetup.h"
#include "PaEvent.h"
#include "PaAlgo.h"
#include "G3part.h"
#include "hepgen_in_phast.h"
#include "WeightInfo.h"
void UserEvent1337(PaEvent &e)
{
#ifdef USE_EXPERIMENTAL
  
  bool debug_print = true;
    //---------------
    Phast &ph = Phast::Ref();
    static bool first(true);
    static WeightInfo myWeightInfo;

    static double weightSum, weightBH, weightDVCS;
    static TTree *tree(NULL);
    static hWeightInterface myProduction;
    static double phi;
    static double qsq;
    static double nu;
    static double xbj, t;
    if (first) { // histograms and Ntupes booking block
        Phast::Ref().HistFileDir("UserEvent1337");
        first = false;
        printf("\n-------HEPGEN - IN - PHAST ----1337: Multi-Weighter!-----------! \n");
        printf("----------------------------------------------------------------! \n");
        printf("------------------------No Histos, just pure reweight!----------! \n");
        printf("------------------------------HEPGEN_IN_PHAST PLUSPLUS----------! \n");

        myProduction.b0 = 4.94166;
        myProduction.alphap = 0.8;
        myProduction.xbj0 = 0.042;

        //hardcoded for generation
        myProduction.numax = 170.;
        myProduction.numin = 2.;
        myProduction.qsqmax = 80.;
        myProduction.qsqmin = 0.5;
        myProduction.tmin = 0.001;
        myProduction.tmax = 1.2;
        myProduction.slpin = 5;
        myProduction.clept = +1.;
        myProduction.slept = -1.;

        static WeightInfo *ptr_userinfo = &myWeightInfo;
        const char *usrclass_name = ptr_userinfo->GetName();
        TTree *tr = ph.out_event_tree;
        if (tr == NULL) {
            return;    // no output was requested (-o option is missing)
        }
        TBranch *b_user = tr->Branch(usrclass_name, usrclass_name, &ptr_userinfo, 32000, 1);
        assert(b_user != NULL);

    } // end of histogram booking
    //---------------


    long Run;
    Run = e.RunNum();




    PYPARS myPars;
    PYSUBS mySubs;
    NLUDATA ld;

    vector<LUJET> myJets;
    int nLujets;
    e.MCgen(nLujets, myJets);
    
    //filter out dd events - we shall not reprocess them here
    if (myJets.size() > 6) {
        return;
    }


    e.MCgen(ld);
    e.MCgen(myPars);
    e.MCgen(mySubs);

    //reset the weightStruct
    myWeightInfo.reset();



    double dvcsWeightnew(0.0);
    double intWeightnew(0.0);
    double  bhWeightnew(0.0);
    //reserve some space for the weights
    double myWeights[7] = {0., 0., 0., 0., 0., 0., 0.};
    
    //******************* ANDRZEJ STYLE WEIGHTS ***************************

    //first we set the struct for this event
    //the third bool is if we want to use nicoles regge params - but we dont so now we leave it false
    prepareWeightInterface(e,myProduction,false);
    
    if (debug_print)
      printf("old_weight: %f %f %f \n",ld.uservar[16],ld.uservar[15],ld.uservar[2]-ld.uservar[16]-ld.uservar[15]);

    
    //andrzej style weights
    if (doDVCSWeights(myProduction,dvcsWeightnew,bhWeightnew,intWeightnew,phi) != 0){
       printf("Problem processing event number:%lld -- Andrzej_originalparams\n", e.UniqueEvNum());
       return;
    }
    
    double allWeightnew = (dvcsWeightnew + intWeightnew + bhWeightnew) ;
    myWeights[0] = allWeightnew;
    myWeights[1] = dvcsWeightnew;
    myWeights[2] = bhWeightnew;
    myWeights[3] = intWeightnew;
    myWeightInfo.addWeight("Andrzej_originalparams", myWeights);

    if (debug_print)
      printf("andrzej_orig: %f %f %f \n",bhWeightnew,dvcsWeightnew,intWeightnew);
    
    //******************* VGG STYLE WEIGHTS (NLO) ***************************

    if (doVGGWeights(myProduction, dvcsWeightnew, bhWeightnew, intWeightnew, phi) != 0) {
        printf("Problem processing event number:%lld -- Moutarde_NLO\n", e.UniqueEvNum());
        return;
    }
    allWeightnew = (dvcsWeightnew + intWeightnew + bhWeightnew) ;
    myWeights[0] = allWeightnew;
    myWeights[1] = dvcsWeightnew;
    myWeights[2] = bhWeightnew;
    myWeights[3] = intWeightnew;
    myWeightInfo.addWeight("VGG_weights", myWeights);

    if (debug_print)
      printf("VGG: %f %f %f \n",bhWeightnew,dvcsWeightnew,intWeightnew);


    //******************* MOUTARDE STYLE WEIGHTS (NLO) ***************************

    if (doMoutardeWeightsWithPF(myProduction, dvcsWeightnew, bhWeightnew, intWeightnew, phi) != 0) {
        printf("Problem processing event number:%lld -- Moutarde_NLO\n", e.UniqueEvNum());
        return;
    }
    allWeightnew = (dvcsWeightnew + intWeightnew + bhWeightnew) ;
    myWeights[0] = allWeightnew;
    myWeights[1] = dvcsWeightnew;
    myWeights[2] = bhWeightnew;
    myWeights[3] = intWeightnew;
    myWeightInfo.addWeight("Moutarde_CX", myWeights);

    if (debug_print)
      printf("moutarde: %f %f %f \n",bhWeightnew,dvcsWeightnew,intWeightnew);

    
    //******************* MOUTARDE STYLE WEIGHTS (LO) ***************************

    if (doMoutardeWeightsWithPF(myProduction, dvcsWeightnew, bhWeightnew, intWeightnew, phi,true) != 0) {
        printf("Problem processing event number:%lld -- Moutarde_LO \n", e.UniqueEvNum());
        return;
    }
    allWeightnew = (dvcsWeightnew + intWeightnew + bhWeightnew) ;
    myWeights[0] = allWeightnew;
    myWeights[1] = dvcsWeightnew;
    myWeights[2] = bhWeightnew;
    myWeights[3] = intWeightnew;
    myWeightInfo.addWeight("Moutarde_CX_LO", myWeights);
    if (debug_print)
      printf("moutardeLO: %f %f %f \n",bhWeightnew,dvcsWeightnew,intWeightnew);

    //******************* ANDRZEJ STYLE WEIGHTS (NICOLES REGGE) ***************************

    prepareWeightInterface(e,myProduction,true);
    if (doDVCSWeights(myProduction,dvcsWeightnew,bhWeightnew,intWeightnew,phi) != 0){
       printf("Problem processing event number:%lld -- Andrzej_originalparams\n", e.UniqueEvNum());
       return;
    }
    allWeightnew = (dvcsWeightnew + intWeightnew + bhWeightnew) ;
    myWeights[0] = allWeightnew;
    myWeights[1] = dvcsWeightnew;
    myWeights[2] = bhWeightnew;
    myWeights[3] = intWeightnew;
    myWeightInfo.addWeight("Andrzej_nicoleparams", myWeights);
    
    if (debug_print)
      printf("andrzej_nicoleregge: %f %f %f \n",bhWeightnew,dvcsWeightnew,intWeightnew);

    //**************** Tag all to save and fin!
    
    e.TagToSave();
#else
    printf("THIS USEREVENT REQUIRES EXPERIMENTAL_HEPGEN_BUILD!!!!!!\n");
#endif






}



