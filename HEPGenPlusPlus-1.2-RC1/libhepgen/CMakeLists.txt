
message("---------------------------")
message("Preparing library libhepgen.so")
message("---------------------------")


####### make the headerlist for installation
FILE(GLOB headerlist *.h)
FILE(GLOB clhepheaderlistrand CLHEP_EMBEDDED/Random/Random/*.h)
FILE(GLOB clhepheaderlistrandicc CLHEP_EMBEDDED/Random/Random/*.icc)
FILE(GLOB clhepheaderlistut CLHEP_EMBEDDED/Utility/Utility/*.h)
FILE(GLOB clhepheaderlistuticc CLHEP_EMBEDDED/Utility/Utility/*.icc)


#### add the mytheo hh stuff
set(headerlist 
        ${headerlist} 
        myTHEO.hh
        )

#### just the normal libs
set(liblist 
	hvector.cpp 
	hlorentzvector.cpp 
	hparammanager.cpp 
	hcardparser.cpp 
	hhelper.cpp 
	hpionicdata.cpp 
	hbooker.cpp 
	hbookbackendASCII.cpp 
	hphysicsgen.cpp 
	hgenmanager.cpp 
	hevent.cpp 
	hparticle.cpp 
	hdvcsgen.cpp 
	houtput.cpp  
	houtputbackend.cpp 
	houtputbackendLEPTO.cpp 
	hrhogen.cpp 
	hpigen.cpp 
	houtputbackendLEPTOGF.cpp 
	hbeamfile.cpp 
	hphigen.cpp 
	hrhoplusgen.cpp 
	homegagenpigamma.cpp
	hjpsigen.cpp
	hrotmat.cpp 
	homegagen.cpp 
	hpamgen.cpp
	reweightKine.cc 
	myTHEO.cc
	)


#### adding CLHEP_EMBEDS
set(liblist 
	${liblist} 
	CLHEP_EMBEDDED/Random/src/Random.cc 
	CLHEP_EMBEDDED/Random/src/RandomEngine.cc 
	CLHEP_EMBEDDED/Random/src/JamesRandom.cc 
	CLHEP_EMBEDDED/Random/src/RanluxEngine.cc 
	CLHEP_EMBEDDED/Random/src/EngineFactory.cc 
	CLHEP_EMBEDDED/Random/src/DoubConv.cc 
	CLHEP_EMBEDDED/Random/src/Ranlux64Engine.cc 
	CLHEP_EMBEDDED/Random/src/StaticRandomStates.cc 
	CLHEP_EMBEDDED/Random/src/engineIDulong.cc 
	CLHEP_EMBEDDED/Random/src/RandGauss.cc 
	CLHEP_EMBEDDED/Random/src/RandGaussQ.cc 
	CLHEP_EMBEDDED/Random/src/RandGaussT.cc 
	CLHEP_EMBEDDED/Random/src/flatToGaussian.cc 
	CLHEP_EMBEDDED/Random/src/erfQ.cc)



if (${USE_EXPERIMENTAL})
set(liblist 
	${liblist}
#  	hVGGgen.cpp
)
endif (${USE_EXPERIMENTAL})



if (${ROOT_FOUND} MATCHES TRUE)
  include_directories (${ROOT_INCLUDE_DIR})
  message(STATUS "Found ROOT - Enabling root-histogram-backend!")
  set(liblist ${liblist} 
	hbookbackendROOT.cpp 
	houtputbackendROOT.cpp
	hMosseGen.cpp)
  set(linklist ${linklist} ${ROOT_LIBRARIES})
  set(USE_ROOT 1)
else()
  message(STATUS "NO ROOT FOUND! Building ultra-basic version of libhepgen.so - NO VGG-MOSSE SUPPORT WILL BE AVAILABLE!")
endif(${ROOT_FOUND} MATCHES TRUE)


#### BUILD THE LIBRARY
add_library(hepgen SHARED ${liblist})

#### LINK TO ROOT
if (${ROOT_FOUND})
  target_link_libraries(hepgen VGG ${ROOT_LIBRARIES})
endif (${ROOT_FOUND})

#### LINK TO VGG AND MOUTARDE
if (${USE_EXPERIMENTAL})
   target_link_libraries(hepgen)
endif(${USE_EXPERIMENTAL})





#### INSTALLATION #####
install (FILES ${headerlist} DESTINATION include)
install (FILES ${clhepheaderlistrand} DESTINATION include/CLHEP_EMBEDDED/Random/Random)
install (FILES ${clhepheaderlistut} DESTINATION include/CLHEP_EMBEDDED/Utility/Utility)
install (FILES ${clhepheaderlistrandicc} DESTINATION include/CLHEP_EMBEDDED/Random/Random)
install (FILES ${clhepheaderlistuticc} DESTINATION include/CLHEP_EMBEDDED/Utility/Utility)
install (TARGETS hepgen LIBRARY DESTINATION lib)
