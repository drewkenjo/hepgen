message("---------------------------")
message("Preparing library python bindings")
message("---------------------------")
FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})


FIND_PACKAGE(PythonLibs REQUIRED)
find_package( PythonInterp REQUIRED )

INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

include_directories (${ROOT_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/../../libhepgen/")
INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/../../libGKPi0/")
INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/../../libVGG/")




SET(CMAKE_SWIG_FLAGS "")

SET_SOURCE_FILES_PROPERTIES(libhepgen_bindings.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(libhepgen_bindings.i PROPERTIES SWIG_FLAGS "")
SWIG_ADD_MODULE(lhepgen python libhepgen_bindings.i )
SWIG_LINK_LIBRARIES(lhepgen hepgen VGG GKPi0 ${PYTHON_LIBRARIES} ${ROOT_LIBRARIES})

execute_process(COMMAND python -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
install(TARGETS _lhepgen DESTINATION $ENV{HEPGENDIR}/install/lib)
install(FILES ${CMAKE_BINARY_DIR}/share/swig/lhepgen.py DESTINATION $ENV{HEPGENDIR}/install/lib)

