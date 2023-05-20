
message("---------------------------")
message("Preparing library libGKPi0.so")
message("---------------------------")

find_package(ROOT)
find_package(GSL)

include_directories (${ROOT_INCLUDE_DIR})
include_directories (${GSL_INCLUDE_DIRS})

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")


FILE(GLOB gkpi0list *.hh)
# FILE(GLOB gkpi0listData data/*)


add_library(GKPi0 SHARED libGKPi0.cpp gkSubProcessTable.cpp)
target_link_libraries(GKPi0 ${ROOT_LIBRARIES} MathMore ${GSL_LIBRARIES})

add_executable(testGKPi0 gkTest.cpp)
target_link_libraries(testGKPi0 GKPi0 ${ROOT_LIBRARIES})

add_executable(gkCrossSection gkCrossSection.cpp)
target_link_libraries(gkCrossSection GKPi0 ${ROOT_LIBRARIES})

add_executable(gkXSbatch gkXSbatch.cpp)
target_link_libraries(gkXSbatch GKPi0 ${ROOT_LIBRARIES})

add_executable(fitHT_prep fitHT_prep.cpp)
target_link_libraries(fitHT_prep GKPi0 ${ROOT_LIBRARIES})

add_executable(fitETbar_prep fitETbar_prep.cpp)
target_link_libraries(fitETbar_prep GKPi0 ${ROOT_LIBRARIES})

add_executable(fitpi0_grid fitpi0_grid.cpp)
target_link_libraries(fitpi0_grid GKPi0 ${ROOT_LIBRARIES})

add_executable(printTheoreticalComparison printTheoreticalComparison.cpp)
target_link_libraries(printTheoreticalComparison GKPi0 ${ROOT_LIBRARIES})

add_executable(printGFFs printGFFs.cpp)
target_link_libraries(printGFFs GKPi0 ${ROOT_LIBRARIES})

add_executable(printETbar printETbar.cpp)
target_link_libraries(printETbar GKPi0 ${ROOT_LIBRARIES})

add_executable(tableNator tableNator.cpp)


install (TARGETS GKPi0 DESTINATION lib)
install (TARGETS gkCrossSection DESTINATION bin)
install (FILES ${gkpi0list} DESTINATION include)
# install (FILES ${gkpi0listData} DESTINATION share/gkpi0data)
