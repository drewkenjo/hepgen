#include "physicsEngine.h"
#include <qt4/QtGui/qapplication.h>


physicsEngineHEPGen::physicsEngineHEPGen()
{
}


physicsEngineHEPGen::~physicsEngineHEPGen()
{
}

void physicsEngineHEPGen::setSettingsAndSeed(int argc, char** argv)
{
    printf("/******---- Using HEPGen library mode! Arguments:\n");
    printf("%s -hepgen [DATACARD] [SEED]\n",argv[0]);
    if (argc < 4)
      QApplication::exit(-1);
    
    QString _fileName = argv[2];
    int seed = atoi(argv[3]);

    string fileName = _fileName.toStdString();

    cout << endl<<  endl;
    for  (int i = 0; i<hepconst::logolength; i++)
        cout << hepconst::logo[i] << endl;
    cout << endl<<  endl;

    HCardParser parserTester;

    tmp = HGenManager::getInstance();
    tmp->loadSettings(fileName);
    tmp->setupGenerator();
    HParamManager ParamManager(fileName);
    tmp->setSeed(seed);


#ifdef USE_ROOT
    tmp->enableROOTBook();
#endif

    tmp->enableASCIIBook();

    bool gfortran = false;


    //TODO: Merge this into HGenManager

    if (ParamManager.getKeyContents("ENABLE_GFORTRAN").at(1) == "1")
        gfortran = true;

    if (ParamManager.getKeyContents("ENABLE_DEBUG").at(1) == "1")
    {
        cout << " Debugging enabled!" << endl;
        tmp->enableDebug();
    }


    tmp->enableLeptoOutput(gfortran);

    string outfile;
    if (ParamManager.getOutName() == "UNSET")
    {
        outfile = "HEPGen++-outfile-";
        outfile.append(currentDateTime());
    }
    else
    {
        outfile = ParamManager.getOutName();
        outfile.append("_bin");
    }

    tmp->enableOutPut(outfile);
    tmp->addStdHistos();


}




std::vector<line> physicsEngineHEPGen::getNextEvent()
{
    std::vector<line> myLines;
    // first we make a valid event

    HGenManager::getInstance()->oneShot();

    while (HGenManager::getInstance()->getEvent()->getStruct()->USERVAR.at(2) <= 0)
        HGenManager::getInstance()->oneShot();

    //then we convert the list of particles to the c-style vertex-arrays we use to push to the GPU
    //the primary vertex will be centered at 0/0/0
    //and the primary muon will be treated differently

    vector<HParticle*>* partList = &HGenManager::getInstance()->getEvent()->getStruct()->listOfParticles;

    //this will result in a buffer-overflow, so better check-crash than to overflow!
    assert(partList->size() < 30 && partList->size() > 0);


    //we take it by half because we want maximal length to be in [-1,0]
    double scaleFac =-3./partList->at(0)->getTVector().Z();
    line tmp;
    float vertex1[3],vertex2[3];
    
    float xOffSet = HGenManager::getInstance()->getEvent()->getStruct()->USERVAR[0]/10.;
    float yOffSet = HGenManager::getInstance()->getEvent()->getStruct()->USERVAR[1]/10.;
    

    vertex1[0] = -3.0;
    vertex1[1] = partList->at(0)->getTVector().X()*scaleFac+xOffSet;
    vertex1[2] = partList->at(0)->getTVector().Y()*scaleFac+yOffSet;


    memset(vertex2,0,sizeof(float)*3);

    vertex2[1] = xOffSet;
    vertex2[2] = yOffSet;
    
    memcpy(&tmp.posVertex1,&vertex1,sizeof(float)*3);
    memcpy(&tmp.posVertex2,&vertex2,sizeof(float)*3);
    tmp.width=1.0;
    tmp.size = sizeof(tmp);
    setBlue(tmp.colorVertex);

    myLines.push_back(tmp);
    for (unsigned int i = 1; i < partList->size(); i++)
    {
        if (partList->at(i)->getParticleAuxFlag() != 1)
            continue;
        //update the scalefactor
        switch(partList->at(i)->getParticleType()) {
        case 22:
            setGreen(tmp.colorVertex);
            break;
        case 13:
            setBlue(tmp.colorVertex);
            break;
        case -13:
            setBlue(tmp.colorVertex);
            break;
        case 2212:
            setRed(tmp.colorVertex);
            break;
        case 2112:
            setPurple(tmp.colorVertex);
            break;
        default:
            setYellow(tmp.colorVertex);
            break;
        }

        scaleFac = 3./partList->at(i)->getTVector().Z();
        vertex1[0] = 3.0;
        vertex1[1] = partList->at(i)->getTVector().X()*scaleFac+xOffSet;
        if (vertex1[1] > 3.0) {
            vertex1[0] /= vertex1[1];
            vertex1[1] = 3.f;
        }

        vertex1[2] = partList->at(i)->getTVector().Y()*scaleFac+yOffSet;
//         vertex1[2] = -6.f;
        memcpy(&tmp.posVertex1[0],&vertex2[0],sizeof(float)*3);
        memcpy(&tmp.posVertex2[0],&vertex1[0],sizeof(float)*3);
        myLines.push_back(tmp);
    }

    return myLines;
}
