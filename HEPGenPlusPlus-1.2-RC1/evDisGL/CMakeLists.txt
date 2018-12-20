
message("---------------------------")
message("Preparing OpenGL event display")
message("---------------------------")

find_package (Qt4)
find_package(GLUT)
find_package(OpenGL)
find_package(ROOT)

set(CMAKE_AUTOMOC TRUE)
set(QT_USE_OPENGL TRUE)


include_directories(${QT_INCLUDE_DIRS})
include_directories(${ROOT_INCLUDE_DIR})
include_directories("../libhepgen/")
include_directories("../tools/")


INCLUDE(${QT_USE_FILE})

ADD_DEFINITIONS(${QT_DEFINITIONS})


add_executable(evDis evDisGL.cpp dialog.cpp glWidget.cpp physicsEngine.cpp physicsEngineLepto.cpp vStructs.h)
target_link_libraries(evDis ${QT_LIBRARY_DIRS} ${QT_QTCORE_LIBRARY}  ${QT_QTGUI_LIBRARY} ${QT_QTOPENGL_LIBRARY} ${OPENGL_LIBRARIES} hepgen slread)