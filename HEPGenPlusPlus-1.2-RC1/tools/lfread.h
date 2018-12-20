#ifndef lfread_h
#define lfread_h

#include "lfcommon.h"


#include <iostream>
#include <fstream>
#include <string>
#include <inttypes.h>

using namespace std;




/*!
 * \namespace lfr
 * \brief the namespace for the lepto-file-reading library lfread
 *
 * These functions should be called in a fixed order:
 * 1.) lLoadFile()
 * 2.) lReadHeader()
 * ----- then often ------------------
 * 	3.) lNextEvent()
 * 	4.) whatEverUserRoutineOrWork
 * 	5.) lFreeBuffer()
 * -----------------------------------
 *
 *
 *
 */
namespace lfr
{


/*! \brief this function loads the file with the fstream provided - also checks for headersize
 *  \param _fileName the full path and name of the leptofile to open
 *  \param fileHandle the fstream to use for opening and reading
 *  \return 1 if success with F77 short headers, 0 if success with GFortran long headers, -1 otherwise.
 */
int lLoadFile(string _fileName, fstream& fileHandle);

/*! \brief this function can be used to check what headersize the given handle has - if it is a correct leptofile.
 *  \param fileHandle the fstream to use - make sure it is open already - also: rewind after this function!
 *  \param returnValue the bool that gets set true if the file is GFortran long headers, and false otherwise.
 *  \return 0 if success, -1 if error
 *
 * Note that one should not need to call that function by hand. It is called automatically by lLoadFile() and its bool is converted into
 * the int-return-value of this function.
 */
int lIsGFortranFile(fstream& fileHandle,bool& returnValue);

/*! \brief loads the next event into the eventStruct
 *  \param fileHandle the fstream to use
 *  \param event the eventBuffer to write into
 *  \param _largeHeaders wether the file is using gfortran or f77 headers
 *  \return 0 if success, -1 if error, -2 if EOF
 *  This function loads the next event from the file and reads it into the event struct.
 *  It allocates the memory for the lBeamParticles inside the event, so make sure to call lFreeBuffer() after use everytime!!!
 */
int lNextEvent(fstream& fileHandle, lfr::lEvent& event, bool _largeHeaders);
/*! \brief deletes the allocated memory of the eventStruct
 * \param event the struct to clean up
 * This function just deletes the allocated memory of the lBeamData array inside the struct.
 * One should always call this function after the work with the struct is finished and before the next lNextEvent() call.
 */
int lFreeBuffer(lEvent& event);
/*! \brief reads the header into the headerstruct
* \param fileHandle the fstream to use
* \param readResult the headerstruct to write into
* \param _largeHeaders wether this file uses gfortran or f77 headers
* \param beQuiet wether we print errors if we find any or just quietly return the values
*
* This function reads the header of the lepto file into the struct.
* It must be called before lNextEvent once!
*/
int lReadHeader(fstream& fileHandle, lfr::lHeader& readResult, bool _largeHeaders, bool _beQuiet = false);
/*! \brief prints the library logo!
 * Do it once, for glory and honor ;)
 */
void lPrintLogo(void);
}



#endif
