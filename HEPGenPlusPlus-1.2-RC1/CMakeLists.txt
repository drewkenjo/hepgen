project (hepgen_cpp)
cmake_minimum_required(VERSION 2.6)

#WOWWOWOW! some serious optimization here - check if machine supports that first!
#-funroll-loops -march=corei7 
#set(CMAKE_CXX_FLAGS "-pg")        ## Optimize
# set(CMAKE_CXX_FLAGS "-O1")        ## Optimize
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-write-strings -Wno-unused-variable -Wno-unused-parameter -Wno-deprecated-declarations -Wno-unused-but-set-parameter -Wno-overloaded-virtual")

#set (CMAKE_CXX_FLAGS "-O3 -funroll-loops -march=corei7")  

#Test for CXX11
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
set(USE_CXX11 0)
set(USE_CXX0X 0)


if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    set(USE_CXX11 1)
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    set(USE_CXX0X 1)
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler for the future!")
endif()


set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/Modules/)

#you will need to comment the next line out for SLC5 and older OS - as these are very old and dying this line will stay in by default - its needed on newer cmake versions
cmake_policy(SET CMP0012 NEW)


#### Print header


message("###-----------------------------##")
message("###HEPGen++ Package installation##")
message("###-----------------------------##")


#### Check input variables 

if (${HEPGEN_DISABLE_ROOT} MATCHES "YES")
message("*** WANT DISABLED ROOT")
set (ROOT_FOUND 0)
else (${HEPGEN_DISABLE_ROOT} MATCHES "YES")
find_package (ROOT)
set (USE_ROOT 1)
endif (${HEPGEN_DISABLE_ROOT} MATCHES "YES")

if (${HEPGEN_ENABLE_PYTHON} MATCHES "YES")
message("*** WANT PYTHON BINDINGS")
set (PYTHON_WANTED 1)
else (${HEPGEN_ENABLE_PYTHON} MATCHES "YES")
set (PYTHON_WANTED 0)
endif (${HEPGEN_ENABLE_PYTHON} MATCHES "YES")




if (${HEPGEN_EXPERIMENTAL_BUILD} MATCHES "YES")
message("*** WANT EXPERIMENTAL BUILD")
set (USE_EXPERIMENTAL 1)
endif (${HEPGEN_EXPERIMENTAL_BUILD} MATCHES "YES")

if (${HEPGEN_EVDISGL} MATCHES "YES")
message("*** WANT OPENGL EVENT DISPLAY")
set (USE_GL 1)
endif (${HEPGEN_EVDISGL} MATCHES "YES")

#### Set config.h file input

CONFIGURE_FILE( ${CMAKE_SOURCE_DIR}/config.h.in ${CMAKE_SOURCE_DIR}/libhepgen/config.h )


#### Set include directories to add
include_directories(hepgen_c)
include_directories(libhepgen)


#### Only can build mosse code or GK Pi0 if root is enabled
if (${ROOT_FOUND} MATCHES TRUE)
include_directories(libVGG)
add_subdirectory(libVGG)

include_directories(libGKPi0)
add_subdirectory(libGKPi0)
endif(${ROOT_FOUND} MATCHES TRUE)

#### add opengl subdirectories
if (${USE_GL})
  add_subdirectory(evDisGL)
endif(${USE_GL})


#### add standard subdirectories
add_subdirectory(libhepgen)
add_subdirectory(hepgen_c)
add_subdirectory(tools)

### add python bindings if wanted
if (${PYTHON_WANTED})
  add_subdirectory(share/swig)
endif (${PYTHON_WANTED})


#### search datacards 
FILE(GLOB datacards resources/datacards/*.data)


#### install misc stuff
install (FILES "share/beamfile.dat" DESTINATION bin)
install (DIRECTORY share/ DESTINATION share)
