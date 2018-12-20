
message("---------------------------")
message("Preparing library libVGG.so")
message("---------------------------")

include_directories (${ROOT_INCLUDE_DIR})
find_package(ROOT)

FILE(GLOB vgglist *.hh)
FILE(GLOB vgglistData data/*)


add_library(VGG SHARED GPDQ.cc mrst99.cc ../libhepgen/myTHEO.cc)
target_link_libraries(VGG ${ROOT_LIBRARIES})

install (TARGETS VGG DESTINATION lib)
install (FILES ${vgglist} DESTINATION include)
install (FILES ${vgglistData} DESTINATION share/vggdata)
