#include "houtput.h"

//TODO
void HOutPut::enableASCII()
{

}

//TODO URGENT!
void HOutPut::enableLEPTO(bool _gfortran)
{
    HOutPutBackEnd* tmp;
    if (_gfortran)
        tmp = new HOutPutBackEndLEPTOGF;
    else
        tmp = new HOutPutBackEndLEPTO;

    tmp->setParams(paramMan);
    backEndList.push_back(tmp);
    enabledLEPTO = true;

}

void HOutPut::enableROOT()
{
#ifdef USE_ROOT
    HOutPutBackEndROOT* tmp = new HOutPutBackEndROOT;
    tmp->setParams(paramMan);
    backEndList.push_back(tmp);
    enabledROOT=true;
#else
    cout << "ROOT OUTPUT DISABLED" << endl;
    enabledROOT = false;
#endif
}


void HOutPut::setFileName(std::string _fileName)
{
    fileName = _fileName;
}

void HOutPut::dumpEvent(HEvent* _event)
{
    for (unsigned int i = 0; i < backEndList.size(); i++)
        backEndList.at(i)->dumpEvent(_event);
}

void HOutPut::closeFile()
{
    for (unsigned int i = 0; i < backEndList.size(); i++)
        backEndList.at(i)->closeFile();
}

void HOutPut::initFile(HEvent* _event)
{
    for (unsigned int i = 0; i < backEndList.size(); i++)
        backEndList.at(i)->initFile(fileName, _event);
}



void HOutPut::writeHeader()
{
    for (unsigned int i = 0; i < backEndList.size(); i++)
        backEndList.at(i)->writeHeader();
}





HOutPut::HOutPut(HParamManager* _params)
{
    paramMan = _params;



    enabledASCII=false;
    enabledROOT=false;
    enabledLEPTO=false;
}

HOutPut::~HOutPut()
{
    for (unsigned int i = 1; i < backEndList.size(); i++)
        delete backEndList.at(i);
    backEndList.clear();

}

