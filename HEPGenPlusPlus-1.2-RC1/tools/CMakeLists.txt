
message("---------------------------")
message("Preparing lepto-utils")
message("---------------------------")




if (${ROOT_FOUND} MATCHES TRUE)
include_directories (${ROOT_INCLUDE_DIR})
# include_directories (${GSL_INCLUDE_DIRS})
# find_package(GSL)
add_executable(makePlotsFromFile makePlotsFromFile.cc)
add_executable(lumiCalc lumiCalc.cc)
add_executable(beamFilePlotter beamfileplotter.cc)
add_executable(leptoFileAnalyzer leptoFileAnalyzer.cpp)

target_link_libraries(beamFilePlotter ${ROOT_LIBRARIES} hepgen)
target_link_libraries(makePlotsFromFile ${ROOT_LIBRARIES})
target_link_libraries(lumiCalc ${ROOT_LIBRARIES} hepgen)
target_link_libraries(leptoFileAnalyzer slread hepgen ${ROOT_LIBRARIES})

install (TARGETS beamFilePlotter DESTINATION bin)
install (TARGETS lumiCalc DESTINATION bin)
install (TARGETS leptoFileAnalyzer DESTINATION bin)
endif (${ROOT_FOUND} MATCHES TRUE)


if (${USE_EXPERIMENTAL} MATCHES TRUE)
add_executable(xSecPlotter xSecPlotter.cc)
target_link_libraries(xSecPlotter hepgen Moutarde VGG cln)
endif(${USE_EXPERIMENTAL} MATCHES TRUE)


add_executable(beamFileSplitter beamFileSplitter.cpp)
target_link_libraries(beamFileSplitter hepgen)


add_library(slread lfread.cpp lfwrite.cpp)


add_executable(leptoFileSplitter splitter.cpp)
target_link_libraries(leptoFileSplitter slread)


install (TARGETS beamFileSplitter DESTINATION bin)
install (TARGETS slread DESTINATION lib)
install (TARGETS leptoFileSplitter DESTINATION bin)
