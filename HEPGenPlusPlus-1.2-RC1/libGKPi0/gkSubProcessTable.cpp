#include "gkSubProcessTable.h"


subProcessAmplitudeLine gkSubProcessTableCache::getLineForKine(double _qsq, double _xbj,int i)
{
    double W = GKPI0::getWsq(_qsq,_xbj);
    double xi = GKPI0::compassxi(_xbj);
    int qbin, wbin;

    for (unsigned int q = 0; q < qbins.size(); q++) {
        qbin = q-1;
        if (_qsq < qbins.at(q))
            break;
    }
    int signQ;
    if (qbin == 0)
      signQ = +1;
    else if (qbin == qbins.size() -1)
      signQ = -1;
    else{
      if ( _qsq-qbins.at(qbin) >   (qbins.at(qbin)-qbins.at(qbin+1))/2.){
	qbin++;
	signQ=-1;
      }
    }
      
      
    for (unsigned int w = 0; w < qbins.size(); w++) {
        wbin = w-1;
        if (W < wbins.at(w))
            break;
    }
    int signW;
    if (wbin == 0)
      signW = +1;
    else if (wbin == wbins.size() -1)
      signW = -1;
    else{
      if ( W-wbins.at(wbin) >   (wbins.at(wbin)-wbins.at(wbin+1))/2.){
	wbin++;
	signW=-1;
      }
    }


    double delta[2] = {signQ*(_qsq-qbins.at(qbin))/(qbins.at(qbin+signQ)-qbins.at(qbin)),
		       signW*(W-wbins.at(wbin))/(wbins.at(wbin+signW)-wbins.at(wbin))};
    vector<double> xLow,xHigh;
    
        double nxlow = 2.*xi*xpseudo.at(i)-xi;
        double nxhigh = (1.-xi)*xpseudo.at(i)+xi;
	subProcessAmplitudeLine deltaLineQ = cacheMap[qbin+signQ][wbin].at(i) - cacheMap[qbin][wbin].at(i);
	subProcessAmplitudeLine deltaLineW = cacheMap[qbin][wbin+signW].at(i) - cacheMap[qbin][wbin].at(i);
	subProcessAmplitudeLine resultLine = cacheMap[qbin][wbin].at(i) + deltaLineQ*delta[0] +deltaLineW*delta[1];
	resultLine.xhigh = nxhigh;
	resultLine.xlow = nxlow;
 	resultLine.weight = weights.at(i);
	//fill the vectors for calling libgkpi0 :D
	xLow.push_back(nxlow);
	xHigh.push_back(nxhigh);
	return resultLine;
}
/*! \brief we mean from the cache, this version does forward and backwards extrapolation! and then call the library to calculate the amplitudes */
GKPI0::amplitude gkSubProcessTableCache::getAmpsForKineDoubleWay(double _qsq, double _xbj, double _t)
{
    GKPI0::amplitude myAmp;
    double W = GKPI0::getWsq(_qsq,_xbj);
    double xi = GKPI0::compassxi(_xbj);

    //check boundaries
//     if (W < wbins.front() || W > wbins.back()) {
//         printf("Fatal error in W range! Must be between %f and %f, was %f \n",wbins.front(),wbins.back(),W);
//         return myAmp;
//     }
//     if (_qsq < qbins.front() || _qsq > qbins.back()) {
//         printf("Fatal error in Qsq range! Must be between %f and %f, was %f \n",qbins.front(),qbins.back(),_qsq);
//         return myAmp;
//     }
//     int qbin, wbin;

//     for (unsigned int q = 0; q < qbins.size(); q++) {
//         qbin = q-1;
//         if (_qsq < qbins.at(q))
//             break;
//     }
//     for (unsigned int w = 0; w < qbins.size(); w++) {
//         wbin = w-1;
//         if (W < wbins.at(w))
//             break;
//     }
    
//         double W = GKPI0::getWsq(_qsq,_xbj);
//     double xi = GKPI0::compassxi(_xbj);
    int qbin, wbin;

    for (unsigned int q = 0; q < qbins.size(); q++) {
        qbin = q-1;
        if (_qsq < qbins.at(q))
            break;
    }
    int signQ;
    if (qbin == 0)
      signQ = +1;
    else if (qbin == qbins.size() -1)
      signQ = -1;
    else{
      if ( _qsq-qbins.at(qbin) >   (qbins.at(qbin)-qbins.at(qbin+1))/2.){
	qbin++;
	signQ=-1;
      }
    }
      
      
    for (unsigned int w = 0; w < qbins.size(); w++) {
        wbin = w-1;
        if (W < wbins.at(w))
            break;
    }
    int signW;
    if (wbin == 0)
      signW = +1;
    else if (wbin == wbins.size() -1)
      signW = -1;
    else{
      if ( W-wbins.at(wbin) >   (wbins.at(wbin)-wbins.at(wbin+1))/2.){
	wbin++;
	signW=-1;
      }
    }


    double delta[2] = {signQ*(_qsq-qbins.at(qbin))/(qbins.at(qbin+signQ)-qbins.at(qbin)),
		       signW*(W-wbins.at(wbin))/(wbins.at(wbin+signW)-wbins.at(wbin))};


//     double delta[2];
//     delta[0] = (_qsq-qbins.at(qbin))/(qbins.at(qbin+1)-qbins.at(qbin));
//     delta[1] = (W-wbins.at(wbin))/(wbins.at(wbin+1)-wbins.at(wbin));
    
//     cout << "Delta Q " << delta[0] << " and Delta W " << delta[1] << endl;
    vector<double> xLow,xHigh;
    vector<TComplex> xLowRes2,xLowRes3,xHighRes2,xHighRes3;
    
    for (int i =0 ; i < 128; i++){
        //calculate new xlow and xhigh for the interpolated values  
        double nxlow = 2.*xi*xpseudo.at(i)-xi;
        double nxhigh = (1.-xi)*xpseudo.at(i)+xi;
	subProcessAmplitudeLine deltaLineQ = cacheMap[qbin+signQ][wbin].at(i) - cacheMap[qbin][wbin].at(i);
	subProcessAmplitudeLine deltaLineW = cacheMap[qbin][wbin+signW].at(i) - cacheMap[qbin][wbin].at(i);
	subProcessAmplitudeLine resultLine = cacheMap[qbin][wbin].at(i);
// 	cout << "LineResult i " << i << " xlowT3 start "<< cacheMap[qbin][wbin].at(i).xlowT3  << " q+1 " << cacheMap[qbin+1][wbin].at(i).xlowT3 << " resulting in " << resultLine.xlowT3 << endl;
	resultLine = resultLine + deltaLineQ*delta[0]*signQ +deltaLineW*delta[1]*signW;
// 	cout << "LineResult i " << i << " xlowT3 start "<< cacheMap[qbin][wbin].at(i).xlowT3  << " q+1 " << cacheMap[qbin+1][wbin].at(i).xlowT3 << " resulting in " << resultLine.xlowT3 << endl;
	resultLine.xhigh = nxhigh;
	resultLine.xlow = nxlow;
 	resultLine.weight = weights.at(i);
	//fill the vectors for calling libgkpi0 :D
	xLow.push_back(nxlow);
	xHigh.push_back(nxhigh);
	xLowRes2.push_back(resultLine.xlowT2);
	xLowRes3.push_back(resultLine.xlowT3);
	xHighRes2.push_back(resultLine.xhighT2);
	xHighRes3.push_back(resultLine.xhighT3);
    }
    GKPI0::loadPreparationFromRam(_xbj,xLowRes3,xHighRes3,xLowRes2,xHighRes2,xLow,xHigh,weights,xpseudo);
    myAmp = GKPI0::getAmplitude(_qsq,xi,_xbj,_t);
    return myAmp;
}

/*! \brief we mean from the cache and then call the library to calculate the amplitudes */
GKPI0::amplitude gkSubProcessTableCache::getAmpsForKine(double _qsq, double _xbj, double _t)
{
    GKPI0::amplitude myAmp;
    double W = GKPI0::getWsq(_qsq,_xbj);
    double xi = GKPI0::compassxi(_xbj);

    //check boundaries
    if (W < wbins.front() || W > wbins.back()) {
        printf("Fatal error in W range! Must be between %f and %f, was %f \n",wbins.front(),wbins.back(),W);
        return myAmp;
    }
    if (_qsq < qbins.front() || _qsq > qbins.back()) {
        printf("Fatal error in Qsq range! Must be between %f and %f, was %f \n",qbins.front(),qbins.back(),_qsq);
        return myAmp;
    }
    int qbin, wbin;

    for (unsigned int q = 0; q < qbins.size(); q++) {
        qbin = q-1;
        if (_qsq < qbins.at(q))
            break;
    }
    for (unsigned int w = 0; w < qbins.size(); w++) {
        wbin = w-1;
        if (W < wbins.at(w))
            break;
    }


    double delta[2];
    delta[0] = (_qsq-qbins.at(qbin))/(qbins.at(qbin+1)-qbins.at(qbin));
    delta[1] = (W-wbins.at(wbin))/(wbins.at(wbin+1)-wbins.at(wbin));
    
//     cout << "Delta Q " << delta[0] << " and Delta W " << delta[1] << endl;
    vector<double> xLow,xHigh;
    vector<TComplex> xLowRes2,xLowRes3,xHighRes2,xHighRes3;
    
    for (int i =0 ; i < 128; i++){
        //calculate new xlow and xhigh for the interpolated values  
        double nxlow = 2.*xi*xpseudo.at(i)-xi;
        double nxhigh = (1.-xi)*xpseudo.at(i)+xi;
	subProcessAmplitudeLine deltaLineQ = cacheMap[qbin+1][wbin].at(i) - cacheMap[qbin][wbin].at(i);
	subProcessAmplitudeLine deltaLineW = cacheMap[qbin][wbin+1].at(i) - cacheMap[qbin][wbin].at(i);
	subProcessAmplitudeLine resultLine = cacheMap[qbin][wbin].at(i);
// 	cout << "LineResult i " << i << " xlowT3 start "<< cacheMap[qbin][wbin].at(i).xlowT3  << " q+1 " << cacheMap[qbin+1][wbin].at(i).xlowT3 << " resulting in " << resultLine.xlowT3 << endl;
	resultLine = resultLine + deltaLineQ*delta[0] +deltaLineW*delta[1];
// 	cout << "LineResult i " << i << " xlowT3 start "<< cacheMap[qbin][wbin].at(i).xlowT3  << " q+1 " << cacheMap[qbin+1][wbin].at(i).xlowT3 << " resulting in " << resultLine.xlowT3 << endl;
	resultLine.xhigh = nxhigh;
	resultLine.xlow = nxlow;
 	resultLine.weight = weights.at(i);
	//fill the vectors for calling libgkpi0 :D
	xLow.push_back(nxlow);
	xHigh.push_back(nxhigh);
	xLowRes2.push_back(resultLine.xlowT2);
	xLowRes3.push_back(resultLine.xlowT3);
	xHighRes2.push_back(resultLine.xhighT2);
	xHighRes3.push_back(resultLine.xhighT3);
    }
    GKPI0::loadPreparationFromRam(_xbj,xLowRes3,xHighRes3,xLowRes2,xHighRes2,xLow,xHigh,weights,xpseudo);
    myAmp = GKPI0::getAmplitude(_qsq,xi,_xbj,_t);
    return myAmp;
}


int gkSubProcessTableCache::loadCache()
{
  string hepPath;
    ifstream gaussFile;
    if (getenv("HEPGEN") != NULL)
        hepPath = getenv("HEPGEN");
    else {
        printf("$HEPGEN env is not set, do not know where to load gauss table values from!\n");
	return -1;
    }
    string fileName = hepPath + "/share/pi0_cache/tableIndex";
    return loadCache(fileName);
}



/*! \brief loads the cache file
 *  This should be formatted as:
 *  ---
 *  QSQ 0.5 1.0 1.5 2.0 [...]
 *  W 1 2 3 4 [...]
 *  DATAFOLDER ./tables/
 *  ---
 *
 * we decided to do 1 < W < 160 and 1 < Qsq < 25
 * but we want to stay flexible
 *
 */

int gkSubProcessTableCache::loadCache(string _fileName)
{
//     printf("Preparing weight table\n");
    string hepPath;
    ifstream gaussFile;
    if (getenv("HEPGEN") != NULL)
        hepPath = getenv("HEPGEN");
    else {
        printf("$HEPGEN env is not set, do not know where to load gauss table values from!\n");
	return -1;
    }
    string fileName = hepPath + "/share/tables/pi0-grid.dat";
    char fNameBuf[200];
    sprintf(fNameBuf,"%s",_fileName.c_str());
    char* dir = dirname(fNameBuf);
//     printf("Loading gauss integration values from %s!\n",fileName.c_str());
    gaussFile.open(fileName, ifstream::in);
    string line;
    while(!gaussFile.eof()) {
        getline(gaussFile,line);
        if (line.length() < 20)
            continue;
        string xVal = line.substr(2,14);
        string weightval = line.substr(18,14);
        // cout << weightval << endl;
        stringstream mystr;
        stringstream weightmystr;

        mystr << xVal;
        double xPseudo;
        mystr >> xPseudo;

        weightmystr << weightval;
        double weighting;
        weightmystr >> weighting;
        weights.push_back(weighting);
	xpseudo.push_back(xPseudo);
    }

//     printf("Initializing cache!\n");
    ifstream myFile;
    myFile.open(_fileName.c_str(),ios::in);
    if (!myFile.good()) {
        printf("Could not open cache index file %s\n",_fileName.c_str());
        return -2;
    }
    string dataFolder = "NA";

    while (!myFile.eof()) {
        string line;
        getline(myFile,line);
        vector<string> boomList = explodeStringWhiteSpace(line);

        if (boomList.size() < 2)
            continue;
        if (boomList.at(0) == "QSQ")
            for (unsigned qCounter = 1; qCounter < boomList.size(); qCounter++)
                qbins.push_back(StrToDouble(boomList.at(qCounter)));
        else if(boomList.at(0) == "W")
            for (unsigned qCounter = 1; qCounter < boomList.size(); qCounter++)
                wbins.push_back(StrToDouble(boomList.at(qCounter)));
        else if(boomList.at(0) == "DATAFOLDER")
            dataFolder = boomList.at(1);
    }
//     printf("Found datafolder %s, qbins %u, wbins %u!\n Trying to load the corresponding files!\n",dataFolder.c_str(),qbins.size(),wbins.size());

    for (unsigned int Q = 0; Q < qbins.size(); Q++)
        for (unsigned int W = 0; W < wbins.size(); W++) {
            double actQsq = qbins.at(Q);
            double actW = wbins.at(W);
            char buffer[200];

            sprintf(buffer,"%s/%s/preparation_%.4f_%.4f.dat",dir,dataFolder.c_str(),actQsq,actW);
            std::ifstream myInFile;
// 	    printf("loading buffer-file %s \n",buffer);
            myInFile.open(buffer,std::ios::in);

            if (!myInFile.good())
                return -3;
            int lines = 0;
            while (!myInFile.eof()) {
                subProcessAmplitudeLine thisLine;
                double buf,buf2;
                myInFile >> buf;
                thisLine.xlow = buf;
                myInFile >> buf;
                thisLine.xlowT3 = TComplex(buf,0.0);
                myInFile >> buf;
                thisLine.xhigh = buf;
                myInFile >> buf;
                myInFile >> buf2;
                thisLine.xhighT3 = TComplex(buf,buf2);
                myInFile >> buf;
                thisLine.xlowT2 = TComplex(buf,0.0);
                myInFile >> buf;
                myInFile >> buf2;
                thisLine.xhighT2 = TComplex(buf,buf2);
                myInFile>> buf;
                thisLine.weight = buf;
		myInFile>> buf;
                thisLine.xpseudo = buf;
		
                lines++;
                cacheMap[Q][W].push_back(thisLine);
            }
        }
//         printf("okay, cache seems fine!\n");
	return 0;
}



