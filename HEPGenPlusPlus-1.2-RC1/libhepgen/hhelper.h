/*!
 * \file hparammanager.h
 *
 * Created in: Hell
 * \author: Christopher Regali
 * Copyright (C) 2010 - All rights reserved
 * \brief this is actually almost a pure copy-paste of some helper functions i wrote in 2010 for a course in the university. (libvdop.so)
 * therefore the commentaries a slighty different in format - but i got better things to do now than to change commentary here.
 */

#ifndef HHELPER_H_
#define HHELPER_H_


/**
* \file hhelper.h
* \brief declares some common functions used in some programs
*/

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <vector>
#include <list>
#include <sstream>
#include <cmath>
#include <time.h>
#include "hconstants.h"
namespace hephelp
{

using namespace std;

/**
 *
 * \brief a function to return the current Date and Time
 *
 */
const string currentDateTime();



/*!
*   \brief   Converts the standart-arguments of a main funktion into an easy to use vector of string!
*/
void StdMainArgsToVecStr(int _argc, char* _argv[], vector<string>& _output);

/*!
* \brief converts a String to Int
*/
int  StrToInt(string _arg);
/*!
* \brief converts an Int to a String
*/
string IntToStr(int _arg);
/*!
* \brief converts a String to a Double
*/
double StrToDouble(string _arg);
/*!
* \brief  Sfilereader2dim reads a table from an ascii-file into vector<vector<T>!
* _dim is the columncount!
*/
template <typename T> void Sfilereader2dim(string _filename, vector<vector<T> >& _data, int _dim=2)
{
    ifstream in_file(_filename.c_str());
    T temp;
    while (1)
    {
        for (int i=0; i < _dim; i++)
        {
            in_file >> temp;
            if (!in_file.good())
                break;
            _data.at(i).push_back(temp);
        }
        if (!in_file.good())
            break;
    }
}




/*! \brief extrapolates a F2 function from a given Q^2 and Xbj */
double F2nmc(double _xbj, double _qsq);


/*! \brief the old F2 extrap function from hepgen-1.20*/
double F2nmc_oldLam(double _xbj, double _qsq);

/*! \brief returns a particle mass by its id */
double getMassByID(int _id);

/*! \brief explodes a string by its whitespaces */
vector<string> explodeStringWhiteSpace(string _input);


/*! \brief explodes a string by custom delimiter */
vector<string> explodeStringCustomDelim ( const string& _input,const string& _delim );







}



#endif
