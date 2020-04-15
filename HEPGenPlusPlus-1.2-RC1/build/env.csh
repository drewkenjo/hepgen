#!/usr/bin/tcsh

#set sourced=($_)
#set HEPGENDIR=`dirname $sourced[2]`
#set HEPGENDIR=`dirname $HEPGENDIR`

source /site/12gev_phys/softenv.csh 2.2
setenv HEPGEN $PWD/../install
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${HEPGEN}/lib
setenv PYTHONPATH      ${PYTHONPATH}:${HEPGEN}/lib
