use root

setenv CVS_RSH ssh
#setenv PYTHONPATH /group/clasdev/local/lib/python
setenv MYSQLINC /usr/include/mysql
setenv MYSQLLIB /usr/lib/mysql

#setenv ROOTINC ${ROOTSYS}/include
#setenv ROOTLIB ${ROOTSYS}/lib
#setenv CLAS_CALIB_RUNINDEX calib_user.Runindexe1_6

setenv PATH ${PATH}::$COATJAVA/bin

source /site/12gev_phys/softenv.csh 2.2

#----  hepgen
setenv HEPGENDIR       $PWD/..
setenv HEPGEN          ${HEPGENDIR}/install
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${HEPGEN}/lib
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/usr/lib64

#------PYTHON3
setenv PYTHON_EXECUTABLE  /apps/python/PRO/bin/python3                                                
setenv PYTHON_INCLUDE_DIR /apps/python/PRO/include/python3.4m/                                         
setenv PYTHON_LIBRARY     /apps/python/PRO/lib/libpython3.4m.so                                                   
setenv LD_PYTHON_LIBRARY  /apps/python/PRO/lib/                                                    
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${LD_PYTHON_LIBRARY}

setenv PYTHONPATH      ${PYTHONPATH}:${HEPGEN}/lib
