/*!
 *  \file hcardparser.h
 *
 *  \date Created on: Jan 16, 2013
 *  \author: Christopher Regali <christopher.regali@cern.ch>
 *  Copyright (c) 2013 All Right Reserved
 */

#ifndef HCPARSER_H_
#define HCPARSER_H_
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>

#include <assert.h>

#include "hhelper.h"

using namespace std;

/*! \brief This class parses a textfile and looks for keywords in there - then returns the values of the keywords
 *
 * somehow a bit like ffread :) but it only reads doubles which is ok for this situation
 *
 */
class HCardParser
{
public:
    /*! \brief standard constructor */
    HCardParser() {};
    /*! \brief standard destructor */
    ~HCardParser() {};

    /*! \brief reads the file in _filename and then parses it. Stores the values in valueList with their respective keys.*/
    void parseFile ( string _filename );


    /*! \brief saves the file in @param _filename with separator @param _separator between segments of the keys */
    void saveToFile(string _fileName, string _separator);


    /*! \brief changes the contents of  @param _key to the new content-vector @param _newContent */
    void setKeyContents(string _key, vector<string> _newContent);


    /*! \brief returns a vector<string> with all the found keys */
    vector<string> getKeyList() const {
        return keyMap;
    };
    /*! \brief returns a vector<string> containing all the substring-values of the found key and the corresponding KEY in position 0 itsself */
    vector<string> getKeyContents(string _key) const;


    /*! \brief parses a generic ffread card without start or stop flags */
    void parseGenericFile ( string _filename );




    /*! \brief resets the class */
    void reset();



private:

    /*! \brief does the actual reading file operation */
    void readFile();

    /*! \brief removes duplicates from keymap */
    void removeDuplicates();


    /*! \brief stores the complete original file */
    vector<string> originalFile;
    /*! \brief The currently used Filename */
    string currentFile;
    /*! \brief the values to the keys, as there are multiple values to each key, each key is assigned a vector of string. (conversion to double or else can be done via sstr so no need to do it here */
    map<string, vector<string> > valueList;
    /*! \brief a list of all the found keys*/
    vector<string> keyMap;


};




#endif
