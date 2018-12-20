#include "physicsEngineLepto.h"
#include <qt4/QtGui/QApplication>

void physicsEngineLEPTO::closeFiles()
{
    outputFile.flush();
    outputFile.close();
}

vector< line > physicsEngineLEPTO::getNextEvent()
{
    if (myEvent.beamParts != NULL)
        lfr::lFreeBuffer(myEvent);
    assert(lfr::lNextEvent(inputFile,myEvent,largeHeaders) != -1);
    assert(!inputFile.eof());

    //print the event first!
    printf("------- Next event loaded!! --------\n");
    for (int i = 0; i < myEvent.nBeamParticle; i++) {
        if (myEvent.beamParts[i].k[0] != 1)
	  continue;
        printf(" %02i / %i: \t %.2e, %.2e, %.2e, %.2e\n",i,
	       myEvent.beamParts[i].k[1],
               myEvent.beamParts[i].p[0],
               myEvent.beamParts[i].p[1],
               myEvent.beamParts[i].p[2],
               myEvent.beamParts[i].p[3]
	      );

    }
    printf("---------------------------------------\n");

    std::vector<line> myLines;
    //then we convert the list of particles to the c-style vertex-arrays we use to push to the GPU
    //the primary vertex will be centered at 0/0/0
    //and the primary muon will be treated differently
    //we take it by half because we want maximal length to be in [-1,0]
    double scaleFac =-3./myEvent.beamParts[0].p[2];
    line tmp;
    float vertex1[3],vertex2[3];
    
    float xOffSet = myEvent.uservar[0]/10.;
    float yOffSet = myEvent.uservar[1]/10.;
    

    vertex1[0] = -3.0;
    //uservars are offsets because of beamfile
    vertex1[1] = myEvent.beamParts[0].p[0]*scaleFac+xOffSet;
    vertex1[2] = myEvent.beamParts[0].p[1]*scaleFac+yOffSet;


    memset(vertex2,0,sizeof(float)*3);

    vertex2[1] = xOffSet;
    vertex2[2] = yOffSet;
    

    memcpy(&tmp.posVertex1,&vertex1,sizeof(float)*3);
    memcpy(&tmp.posVertex2,&vertex2,sizeof(float)*3);
    tmp.width=1.0;
    tmp.size = sizeof(tmp);
    setBlue(tmp.colorVertex);

    myLines.push_back(tmp);
    for (unsigned int i = 1; i < myEvent.nBeamParticle; i++)
    {
        if (myEvent.beamParts[i].k[0] != 1)
            continue;
        //update the scalefactor
        switch(abs(myEvent.beamParts[i].k[1])) {
        case 22:
            setGreen(tmp.colorVertex);
            break;
        case 13:
            setBlue(tmp.colorVertex);
            break;
        case 2212:
            setRed(tmp.colorVertex);
            break;
        case 2112:
            setPurple(tmp.colorVertex);
            break;
	case 130: case 321: case 310: setOrange(tmp.colorVertex);break;    
        default:
            setYellow(tmp.colorVertex);
            break;
        }

        scaleFac = 3./myEvent.beamParts[i].p[2];
        vertex1[0] = 3.0;
        vertex1[1] = myEvent.beamParts[i].p[0]*scaleFac+xOffSet;
    ;
        if (vertex1[1] > 3.0) {
            vertex1[0] /= vertex1[1];
            vertex1[1] = 3.f;
        }

        vertex1[2] = myEvent.beamParts[i].p[1]*scaleFac+yOffSet;

        memcpy(&tmp.posVertex1[0],&vertex2[0],sizeof(float)*3);
        memcpy(&tmp.posVertex2[0],&vertex1[0],sizeof(float)*3);
        myLines.push_back(tmp);
    }

    return myLines;
}

void physicsEngineLEPTO::setSettingsAndSeed(int argc, char** argv)
{
    printf("/******---- Using LEPTOFILE  mode! Arguments:\n");
    printf("%s -leptofile [LEPTO_FILE_IN] [LEPTO_FILE_OUT]\n",argv[0]);
    printf("/******---- (powered by slread leptofile-reader-library)):\n");
    
    QString _inFile = argv[2];
    QString _outFile = argv[3];
    lfr::lPrintLogo();
    int returnCode = lfr::lLoadFile(_inFile.toStdString(),inputFile);
    if (returnCode == -1) {
        printf("------- Error loading the file!\n");
        QApplication::exit(-1);
    }
    else if(returnCode == 0) {
        printf("------- Found good LEPTOFILE with GFORTRAN headers \n");
    }
    else if (returnCode == 1) {
        printf("------- Found good LEPTOFILE with F77 headers \n");
    }

    largeHeaders = (returnCode == 0);

    //----------------Read the header of the file---------
    if (lfr::lReadHeader(inputFile,myHeader,largeHeaders) == -1)
    {
        printf("------- Error while reading header of the LEPTOFILE \n");
        QApplication::exit(-1);
    }

    printf("------- Header has been read successfully!! -------- \n");

    int returnWriter = lfr::lOpenFile(_outFile.toStdString(),outputFile);
    assert(returnWriter == 0);
    lfr::lWriteHeader(outputFile,myHeader);
    //important because we need to know if it has been initialized
    myEvent.beamParts = NULL;

}

void physicsEngineLEPTO::triggerTheSave()
{
    lfr::lWriteNextEvent(outputFile,&myEvent);
    outputFile.flush();
    printf("Event saved!\n");
}
