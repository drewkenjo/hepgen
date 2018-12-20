#ifndef lfwrite_h
#define lfwrite_h
#include "lfcommon.h"



#include <iostream>
#include <fstream>
#include <string>
#include <inttypes.h>


using namespace std;
namespace lfr {

/*! \brief this function prepares the filehandle for ascii writing of the leptofile
 *  \param _fileName the full path and name of the leptofile to open
 *  \param fileHandle the fstream to use for opening and writing
 *  \param _useLongHeaders wether to use F77 short or gfortran long headers
 *  \return 1 if success with F77 short headers, 0 if success with GFortran long headers, -1 otherwise.
 */
int lOpenFile(string _fileName, fstream& fileHandle, bool _useLongHeaders = false);

/*! \brief writes the header from the headerstruct into the file
  * \param fileHandle the fstream to use
  * \param readResult the headerstruct to write into the file
  * \param _largeHeaders wether this file uses gfortran or f77 headers
  *
  * This function writes the header into the lepto file from the struct.
  * It must be called before lWriteNextEvent once!
  */
int lWriteHeader(fstream& fileHandle, lfr::lHeader& writeMe, bool _largeHeaders = false);

/*! \brief writes the eventstruct to a file
 *  \param fileHandle the fstream to use
 *  \param event the eventBuffer to write
 *  \param _largeHeaders wether the file is using gfortran or f77 headers
 *  \return 0 if success, -1 if error, -2 if EOF
 *  This function writes the next event from the bufferstruct into the file.
 */
int lWriteNextEvent(fstream& fileHandle, lfr::lEvent* event, bool _largeHeaders = false);


void putFloat(fstream& fileHandle, float _f);
void putInt(fstream& fileHandle, int32_t _i);
void putInt64(fstream& fileHandle, int64_t _i);
}
#endif
