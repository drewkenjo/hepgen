#ifndef anaPlugin_H
#define anaPlugin_H

#include <iostream>
#include <string>

#include "lfread.h"


using namespace std;

class anaPlugin
{
public:
    anaPlugin();
    virtual ~anaPlugin();
    string getName() {
        return name;
    };
    virtual void analyze() = 0;


private:
    string name;



};


#endif

