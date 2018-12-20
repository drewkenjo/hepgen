/*!
 *  \file hbookbackendASCII.h
 *
 *  \date Created on: Jan 29, 2013
 *  \author: Christopher Regali  <christopher.regali@cern.ch>
 *
 *  Copyright (c) 2013 All Right Reserved
 */


#ifndef HBOOKBACKENDASCII_H_
#define HBOOKBACKENDASCII_H_


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <cstring>

#include "hbookbackend.h"


using namespace std;

class ASCIIhist;
class ASCIIhistd;
class ASCIIhist2d;


/*!
 * \class HBookBackEndASCII
 * \brief A simple ascii-Backend for histogram saving
 *
 * This is a simple but powerful ASCII implementation of Histograms.
 * The histograms can be used in a similar way that one would use the root histos.
 * A histogram is added to the lia via an add function.
 */
class HBookBackEndASCII: public HBookBackEnd
{
public:
    HBookBackEndASCII(std::string _nameSuffix):HBookBackEnd(_nameSuffix) {};
    /*! \brief destructor, cleans up the heap */
    ~HBookBackEndASCII();

    /*! \brief adds a double-typed histogram to the list, returns the new index else
     *  \param _var a pointer to the variable, will be read out with each fill call
     *  \param _weight a pointer to the weight-variable, will be read out with each fill call
     *  \param _binCount the number of bins
     *  \param _lowerBound the lower boundary of the histogram
     *  \param _upperBound the upper boundary of the histogram
     *  \param _title the title that this histogram is going to have
     *  \return the number of this histogram in the list after adding it. This can be used for single-filling.
     */
    int addHBook1D(double* _var, double* _weight, int _binCount, double _lowerBound, double _upperBound, string _title);

    /*! \brief add a 2-dim Double-Histogram
     *  \param _var a pointer to the variable, will be read out with each fill call
     *  \param _weight a pointer to the weight-variable, will be read out with each fill call
     *  \param _binCount the number of bins
     *  \param _lowerBound the lower boundary of the histogram
     *  \param _upperBound the upper boundary of the histogram
     *  \param _title the title that this histogram is going to have
     *  \return the number of this histogram in the list after adding it. This can be used for single-filling.
     *
     */
    int addHBook2D(double* _varX, double* _varY , double* _weight, int _binCountX,int _binCountY , double _lowerBoundX,double _lowerBoundY, double _upperBoundX,double _upperBoundY, std::string _title);

    
    
    
    /*! \brief adds a float-typed histogram to the list, returns the new index else
     *  \param _var a pointer to the variable, will be read out with each fill call
     *  \param _weight a pointer to the weight-variable, will be read out with each fill call
     *  \param _binCount the number of bins
     *  \param _lowerBound the lower boundary of the histogram
     *  \param _upperBound the upper boundary of the histogram
     *  \param _title the title that this histogram is going to have
     *  \return the number of this histogram in the list after adding it. This can be used for single-filling.
     */
    int addHBook1F(float* _var, float* _weight, int _binCount, float _lowerBound, float _upperBound, string _title);

    /*! \brief add a 2-dim Float-Histogram
     *  \param _var a pointer to the variable, will be read out with each fill call
     *  \param _weight a pointer to the weight-variable, will be read out with each fill call
     *  \param _binCount the number of bins
     *  \param _lowerBound the lower boundary of the histogram
     *  \param _upperBound the upper boundary of the histogram
     *  \param _title the title that this histogram is going to have
     *  \return the number of this histogram in the list after adding it. This can be used for single-filling.
     *
     */
    int addHBook2F(float* _varX, float* _varY , float* _weight, int _binCountX,int _binCountY , float _lowerBoundX,float _lowerBoundY, float _upperBoundX,float _upperBoundY, std::string _title);

    
    

    /*! \brief NOT IMPLEMENTED YET! (dont know if it ever will be!) */
    int addHBook1I(int* _var, double* _weight, int _binCount, int _lowerBound, int _upperBound, string _title);

    /*!  \brief saves the histograms to a file with name _fileName. Extension .dat will be added hardcoded to make the backends use different files */
    void dumpToFile(string _fileName);
    /*! \brief fills all histograms */
    void fill();
    /*! \brief fills only histogram # _histNum */
    int fill(int _histNum);
    /*! \brief returns title of histogram if _num is in range, returns "_ERROR_OUT_OF_RANGE" else */
    string getTitle(int _num);

private:
    vector<ASCIIhist*>  histList;

};


/*! \class ASCIIhist
 *  \brief A simple ascii histogram class
 *
 * This class is simple implementation of a histogram featuring only filling and saving to a file.
 * Its very basic and does not feature anything sophisticated as dividing or statistical analysis.
 * Note that this is only the base-class. All specifics histograms are derived from this.
 */

class ASCIIhist
{
public:
    /*! \brief puts the data in ascii-format to the given ofstream
     *
     */
    virtual void dumpToFile(ofstream& _file) =0;

    /*! \brief returns the title */
    string getTitle() const {
        return title;
    };
    virtual int fill()=0;

protected:
    string title;
};


/*! \class ASCIIhist2d
 * \brief a 2 dimensional double-typed histogram with ascii-writing
 * This class is an implementation of a 2 dimensional double-type histogram with write-support to an ASCII file.
 * It features weighting.
 */

class ASCIIhist2d : public ASCIIhist
{
public:

    /*! \brief the constructor
     * \param _varX the pointer to the variable on the x-axis, will be read out each fill
     * \param _varY the pointer to the variable on the y-axis, will be read out each fill
     * \param _weight the pointer to the weight, will be read out each fill
     * \param _binCountX the binCount on the x-axis
     * \param _binCountY the binCount on the y-axis
     * \param _lowerBoundX the lower boundary of the x-axis
     * \param _lowerBoundY the lower boundary of the y-axis
     * \param _upperBoundX the upper boundary of the x-axis
     * \param _upperBoundY the upper boundary of the y-axis
     * \param _title the title this histogram is going to have
     */
    ASCIIhist2d(double* _varX, double* _varY , double* _weight, int _binCountX,int _binCountY , double _lowerBoundX,double _lowerBoundY, double _upperBoundX,double _upperBoundY, std::string _title);
    ~ASCIIhist2d() {};

    /*! \brief resets the data in the dataholder */
    void reset();
    /*! \brief puts the data in ascii-format to the given ofstream */
    void dumpToFile(ofstream& _file);
    /*! \brief fills the histogram and returns 0 always! */
    int fill();


private:
    double* variableX;
    double* variableY;
    double* weight;
    vector< vector< double > > histData;
    double lowerBoundX,upperBoundX,stepLengthX;
    double lowerBoundY,upperBoundY,stepLengthY;
    int binCountX;
    int binCountY;
    string title;
    int fillCount;
};





/*!
 * \class ASCIIhistd
 * \brief a simple histogram dataholder class
 *
 * Could actually be somewhat more sophisticated, but does the job here. If ever needed there could easily be more functions implemented later on.
 * Features only weighting, titling and writing to ASCII right now
 *
 */

class ASCIIhistd : public ASCIIhist
{
public:

    ASCIIhistd(double* _var, double* _weight, int _binCount, double _lowerBound, double _upperBound, string _title);
    ~ASCIIhistd() {};

    /*! \brief resets the data in the dataholder */
    void reset(double* _var, double* _weight, int _binCount, double _lowerBound, double _upperBound, string _title);
    /*! \brief puts the data in ascii-format to the given ofstream */
    void dumpToFile(ofstream& _file);
    /*! \brief fills the histogram and returns the binnumber that was stocked */
    int fill();





private:
    double* variable;
    double* weight;
    map<int,double> histData;
    double lowerBound,upperBound,stepLength;
    int binCount;

    int fillCount;

};


/*! \class ASCIIhist2f
 * \brief a 2 dimensional float-typed histogram with ascii-writing
 * This class is an implementation of a 2 dimensional float-type histogram with write-support to an ASCII file.
 * It features weighting.
 */

class ASCIIhist2f : public ASCIIhist
{
public:

    /*! \brief the constructor
     * \param _varX the pointer to the variable on the x-axis, will be read out each fill
     * \param _varY the pointer to the variable on the y-axis, will be read out each fill
     * \param _weight the pointer to the weight, will be read out each fill
     * \param _binCountX the binCount on the x-axis
     * \param _binCountY the binCount on the y-axis
     * \param _lowerBoundX the lower boundary of the x-axis
     * \param _lowerBoundY the lower boundary of the y-axis
     * \param _upperBoundX the upper boundary of the x-axis
     * \param _upperBoundY the upper boundary of the y-axis
     * \param _title the title this histogram is going to have
     */
    ASCIIhist2f(float* _varX, float* _varY , float* _weight, int _binCountX,int _binCountY , float _lowerBoundX,float _lowerBoundY, float _upperBoundX,float _upperBoundY, std::string _title);
    ~ASCIIhist2f() {};

    /*! \brief resets the data in the dataholder */
    void reset();
    /*! \brief puts the data in ascii-format to the given ofstream */
    void dumpToFile(ofstream& _file);
    /*! \brief fills the histogram and returns 0 always! */
    int fill();


private:
    float* variableX;
    float* variableY;
    float* weight;
    vector< vector< float > > histData;
    float lowerBoundX,upperBoundX,stepLengthX;
    float lowerBoundY,upperBoundY,stepLengthY;
    int binCountX;
    int binCountY;
    string title;
    int fillCount;
};





/*!
 * \class ASCIIhistf
 * \brief a simple histogram dataholder class float type
 *
 * Could actually be somewhat more sophisticated, but does the job here. If ever needed there could easily be more functions implemented later on.
 * Features only weighting, titling and writing to ASCII right now
 *
 */

class ASCIIhistf : public ASCIIhist
{
public:

    ASCIIhistf(float* _var, float* _weight, int _binCount, float _lowerBound, float _upperBound, string _title);
    ~ASCIIhistf() {};

    /*! \brief resets the data in the dataholder */
    void reset(float* _var, float* _weight, int _binCount, float _lowerBound, float _upperBound, string _title);
    /*! \brief puts the data in ascii-format to the given ofstream */
    void dumpToFile(ofstream& _file);
    /*! \brief fills the histogram and returns the binnumber that was stocked */
    int fill();





private:
    float* variable;
    float* weight;
    map<int,float> histData;
    float lowerBound,upperBound,stepLength;
    int binCount;

    int fillCount;

};









#endif
