
message("---------------------------")
message("Preparing HEPGen++ binary")
message("---------------------------")


set(linklist hepgen)
set(CURSES_NEED_NCURSES TRUE)

find_package(Curses)

if (${ROOT_FOUND} MATCHES TRUE)
  include_directories (${ROOT_INCLUDE_DIR})
  set(linklist ${linklist} ${ROOT_LIBRARIES})
  set(USE_ROOT 1)
else()
  message(STATUS "NO ROOT FOUND! Building ultra-basic version of libhepgen.so")
endif(${ROOT_FOUND} MATCHES TRUE)

if (${CURSES_FOUND} MATCHES TRUE)
  add_executable(HEPGenCurses++ main_curses.cpp)
  target_link_libraries(HEPGenCurses++ ${linklist} ${CURSES_LIBRARIES})
  install (TARGETS HEPGenCurses++ RUNTIME DESTINATION bin)
endif(${CURSES_FOUND} MATCHES TRUE)

link_directories(${CMAKE_BINARY_DIR})


add_executable(HEPGen++ main.cpp)
target_link_libraries(HEPGen++ ${linklist})

install (TARGETS HEPGen++ RUNTIME DESTINATION bin)
