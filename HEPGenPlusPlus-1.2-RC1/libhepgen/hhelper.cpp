
#include "hhelper.h"
using namespace std;

/**
* \file hhelper.cpp
* \brief implements some common functions used in some programs
*/
namespace hephelp
{
/**
 *  \brief    Converts the standart-arguments of a main funktion into an easy to use vector of string!
 */
void StdMainArgsToVecStr(int _argc, char* _argv[], std::vector<std::string>& _output)
{
    stringstream sstr;
    for (int i = 0; i < _argc; i++)
    {
        sstr << _argv[i];
        _output.push_back(sstr.str());
        sstr.flush();
        sstr.str(string());
    }
}


/**
 *  \brief    Returns a particle mass from a particle id
 */

double getMassByID(int _id)
{
    switch (_id)
    {
    case hepconst::typeElectronMinus:
        return sqrt(hepconst::w2electron);
        break;
    case hepconst::typeElectronPlus:
        return sqrt(hepconst::w2electron);
        break;
    case hepconst::typeMuonMinus:
        return sqrt(hepconst::w2mu);
        break;
    case hepconst::typeMuonPlus:
        return sqrt(hepconst::w2mu);
        break;
    case hepconst::typeGamma:
        return 0.0;
        break;
    case hepconst::typeProton:
        return sqrt(hepconst::w2proton);
        break;

    case hepconst::typeNeutron:
        return sqrt(hepconst::w2neutron);
        break;

    case hepconst::typePiMinus:
        return sqrt(hepconst::w2pic);
        break;
    case hepconst::typePiPlus:
        return sqrt(hepconst::w2pic);
        break;
    case hepconst::typePi0:
        return sqrt(hepconst::w2pi);
        break;

    case hepconst::typeKPlus:
        return sqrt(hepconst::w2KaonCharged);
        break;
    case hepconst::typeKMinus:
        return sqrt(hepconst::w2KaonCharged);
        break;
    case hepconst::typePhi:
        return hepconst::mPhi;
        break;


    }
    std::cout << "WARNING! PARTICLE MASS LOOKUP FAILED FOR PID: " << _id << std::endl;
    return -1;
}




/**
 * \brief Converts an Int to a String
 */

string IntToStr(int _arg)
{

    stringstream sstr;
    sstr << _arg;
    return sstr.str();
}

/**
 * \brief Converts a String to an Int
 */
int StrToInt(string _arg)
{
    stringstream sstr;
    int i;
    sstr << _arg;
    sstr >> i;
    return i;
}

/**
 * \brief Converts a String to a Double
 */
double StrToDouble(string _arg)
{
    stringstream sstr;
    double i;
    sstr << _arg;
    sstr >> i;
    return i;
}



double F2nmc(double _xbj, double _qsq)
{
    double apar[] = {-0.02778, 2.926, 1.0362, -1.840, 8.123, -13.074, 6.215};
    double bpar[] = {0.285, -2.694, 0.0188, 0.0274};
    double cpar[] = {-1.413, 9.366, -37.79, 47.10};
    double Q02 = 20.;
    //corrected lam2 factor
    double LAM2 = 0.0625;


    double A = pow(_xbj,apar[0]) * pow((1-_xbj),apar[1]) *
               (apar[2] + apar[3] * (1-_xbj) + apar[4] * pow((1-_xbj),2)
                + apar[5] * pow((1-_xbj),3) + apar[6]*pow((1-_xbj),4));
    double B = bpar[0] + bpar[1] * _xbj + bpar[2] / (_xbj + bpar[3]);
    double C = cpar[0] * _xbj + cpar[1] * pow(_xbj,2) + cpar[2] * pow(_xbj,3) + cpar[3] * pow(_xbj,4);
    double SFF2 = A * pow((log(_qsq / LAM2) / log(Q02 / LAM2)), B) * (1 + C / _qsq);
    return SFF2;
}

double F2nmc_oldLam(double _xbj, double _qsq)
{
    double apar[] = {-0.02778, 2.926, 1.0362, -1.840, 8.123, -13.074, 6.215};
    double bpar[] = {0.285, -2.694, 0.0188, 0.0274};
    double cpar[] = {-1.413, 9.366, -37.79, 47.10};
    double Q02 = 20.;
    double LAM2 = 0.250;


    double A = pow(_xbj,apar[0]) * pow((1-_xbj),apar[1]) *
               (apar[2] + apar[3] * (1-_xbj) + apar[4] * pow((1-_xbj),2)
                + apar[5] * pow((1-_xbj),3) + apar[6]*pow((1-_xbj),4));
    double B = bpar[0] + bpar[1] * _xbj + bpar[2] / (_xbj + bpar[3]);
    double C = cpar[0] * _xbj + cpar[1] * pow(_xbj,2) + cpar[2] * pow(_xbj,3) + cpar[3] * pow(_xbj,4);
    double SFF2 = A * pow((log(_qsq / LAM2) / log(Q02 / LAM2)), B) * (1 + C / _qsq);
    return SFF2;
}

const string currentDateTime()
{
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
    return buf;
}

vector<string> explodeStringWhiteSpace(string _input)
{
    vector<string> array;

    int inputlength = _input.length();
    int whitespaceLength = 1;
    int i = 0;
    int k = 0;
    while (i < inputlength) {
        int j = 0;
        while (i + j < inputlength && j < whitespaceLength
                && (_input[i + j] == ' ' || _input[i + j] == '\t'))
            j++;
        if (j == whitespaceLength) { //found whitespace
            if (_input.substr(k, i - k) != "")
                array.push_back(_input.substr(k, i - k));
            i += whitespaceLength;
            k = i;
        } else {
            i++;
        }
    }
    if (_input.substr(k, i - k) != "")
        array.push_back(_input.substr(k, i - k));
    return array;
}

vector<string> explodeStringCustomDelim ( const string& _input,const string& _delim )
{
    vector<string> arr;

    int _inputleng = _input.length();
    int delleng = _delim.length();
    if ( delleng==0 )
        return arr;//no change

    int i=0;
    int k=0;
    while ( i<_inputleng )
    {
        int j=0;
        while ( i+j<_inputleng && j<delleng && _input[i+j]==_delim[j] )
            j++;
        if ( j==delleng )
        {   //found _delim
            if (_input.substr ( k, i-k ) != "")
                arr.push_back ( _input.substr ( k, i-k ) );
            i+=delleng;
            k=i;
        }
        else
        {
            i++;
        }
    }
    if (_input.substr ( k, i-k ) != "" && _input.substr ( k, i-k ) != " " && _input.substr ( k, i-k ) != "\t")
        arr.push_back ( _input.substr ( k, i-k ) );
    return arr;
}





}
