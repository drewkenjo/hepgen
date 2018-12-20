// $Id:
// -*- C++ -*-
//
// -----------------------------------------------------------------------
//                             HEP Random
//                      --- StaticRandomStates ---
//                      class implementation file
// -----------------------------------------------------------------------
//
// =======================================================================
// Mark Fischler  - Created: Dec. 21, 2004
// Mark Fischler  - Modified restore() to utilize anonymous engine input
//                  to create anonymous restore of the static distributions
//
// =======================================================================

#include "CLHEP_EMBEDDED/Random/Random/StaticRandomStates.h"
#include <string>
#include <sstream>

//======================//
//                      //
// Maintenance warning: //
//			//
//======================//
//
// Currently, only two distributions (RandFlat and RandGauss) have cached
// distribution state.  All such distributions must be saved below, so if
// another such distribution is added, this implementation file must be
// modified to reflect that.

namespace CLHEP {


std::ostream & StaticRandomStates::save(std::ostream & os) {
    return os;
}

#ifdef NOTYET
std::istream & StaticRandomStates::restore(std::istream & is) {
    RandGauss::restoreFullState(is);
    RandFlat::restoreDistState(is);
    return is;
}
#endif

std::istream & StaticRandomStates::restore(std::istream & is) {
    if ( !is ) return is;
    return is;
}

}  // namespace CLHEP
