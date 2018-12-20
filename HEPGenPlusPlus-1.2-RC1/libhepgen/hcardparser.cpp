/*
 * hcardparser.cpp
 *
 *  Created on: Jan 16, 2013
 *      Author: Christopher Regali
 *
 *
 */
#include "hcardparser.h"


void HCardParser::removeDuplicates()
{
    keyMap.erase(
        unique(keyMap.begin(), keyMap.end()),
        keyMap.end());
    sort(keyMap.begin(), keyMap.end());
}


void HCardParser::setKeyContents(string _key, vector< string > _newContent)
{
    valueList[_key]= _newContent;
    keyMap.push_back(_key);
    removeDuplicates();
}





void HCardParser::saveToFile(string _fileName, string _separator)
{
    ofstream myfile (_fileName.c_str());
    assert (myfile.is_open());
    map< string,vector<string> >::iterator mapIter;

    for (mapIter=valueList.begin(); mapIter != valueList.end(); mapIter++)
    {
        for (unsigned int i = 0; i < mapIter->second.size(); i++)
            myfile << mapIter->second.at(i) << _separator;

        myfile << endl;
    }
}







void HCardParser::reset()
{
    originalFile.clear();
    valueList.clear();
    keyMap.clear();

    currentFile="";
}

void HCardParser::readFile()
{
    string line;
    ifstream myfile ( currentFile.c_str() );
    if ( myfile.is_open() ) {
        while ( myfile.good() ) {
            getline ( myfile,line );
            originalFile.push_back ( line );
        }
        myfile.close();
    } else
        cout<<"HCardParser could not open File: "<< currentFile << " - That is a big Problem!" << endl;
}




vector< string > HCardParser::getKeyContents(string _key) const
{
  if (valueList.find(_key) != valueList.end())
    return valueList.find(_key)->second;
  else{
    vector<string> tmp;
    tmp.push_back(_key);
    tmp.push_back("UNSET");
    return tmp;
  }
}



void HCardParser::parseGenericFile ( string _filename )
{

    currentFile = _filename;
    readFile();
    for (unsigned int i=0; i < originalFile.size() ; i++ ) {
        //commented lines start with a * and are completely ignored
        if ( originalFile.at ( i ).substr(0,1) == "*" )
            continue;
        if ( originalFile.at ( i ).substr(0,1) == "#" )
            continue;
        if ( originalFile.at ( i ).substr(0,2) == "//" )
            continue;
        //the rest gets exploded
        vector<string> boomResult = hephelp::explodeStringWhiteSpace(originalFile.at(i));
        //then put into the map.
        if (boomResult.size() > 0)
        {
            valueList[boomResult.at(0)]=boomResult;
            keyMap.push_back(boomResult.at(0));
        }
    }
    removeDuplicates();


}


void HCardParser::parseFile ( string _filename )
{
    currentFile = _filename;
    readFile();
    bool inList = false;
    for (unsigned int i=0; i < originalFile.size() ; i++ ) {
        //commented lines start with a * and are completely ignored
        if ( originalFile.at ( i ).substr(0,1) == "*" )
            continue;

        //we only want to parse contents from BEGIN to END
        if (originalFile.at(i) == "LIST")
        {
            inList = true;
            continue;

        }
        if (originalFile.at(i) == "END")
        {
            inList = false;
        }


        if (inList)
        {
            //the rest gets exploded
            vector<string> boomResult = hephelp::explodeStringWhiteSpace(originalFile.at(i));
            //then put into the map.
            if (boomResult.size() > 0)
            {
                valueList[boomResult.at(0)]=boomResult;
                keyMap.push_back(boomResult.at(0));
            }
        }
    }
    removeDuplicates();
}

