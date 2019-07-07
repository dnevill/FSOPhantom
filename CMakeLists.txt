# $Id: CMakeLists.txt 86065 2014-11-07 08:51:15Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(FSOPhantom)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(FSOPhantom FSOPhantom.cc ${sources} ${headers})
target_link_libraries(FSOPhantom ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build FSOPhantom. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#

set(FSOPhantom_SCRIPTS
  init_vis.mac
  vis.mac
  )
 if (WIN32)
set(FSOPhantom_WinDLLs
  libEGL.dll
  libGLESV2.dll
  opengl32sw.dll
  Qt5Core.dll
  Qt5Gui.dll
  Qt5OpenGL.dll
  Qt5Svg.dll
  Qt5Widgets.dll
  iconengines/qsvgicon.dll
  imageformats/qgif.dll
  imageformats/qicns.dll
  imageformats/qico.dll
  imageformats/qjpeg.dll
  imageformats/qsvg.dll
  imageformats/qtga.dll
  imageformats/qtiff.dll
  imageformats/qwbmp.dll
  imageformats/qwebp.dll
  platforms/qwindows.dll
  styles/qwindowsvistastyle.dll
  translations/qt_ar.qm
  translations/qt_bg.qm
  translations/qt_ca.qm
  translations/qt_cs.qm
  translations/qt_da.qm
  translations/qt_de.qm
  translations/qt_en.qm
  translations/qt_es.qm
  translations/qt_fi.qm
  translations/qt_fr.qm
  translations/qt_gd.qm
  translations/qt_he.qm
  translations/qt_hu.qm
  translations/qt_it.qm
  translations/qt_ja.qm
  translations/qt_ko.qm
  translations/qt_lv.qm
  translations/qt_pl.qm
  translations/qt_ru.qm
  translations/qt_sk.qm
  translations/qt_uk.qm

  
  )
  foreach(_script ${FSOPhantom_WinDLLs})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/dependencies/${_script}
    COPYONLY
    )
endforeach()
  endif (WIN32)
foreach(_script ${FSOPhantom_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()
if (WIN32)
add_custom_command(TARGET FSOPhantom POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_directory
		"${PROJECT_BINARY_DIR}/dependences"
		$<TARGET_FILE_DIR:FSOPhantom>)
endif (WIN32)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS FSOPhantom DESTINATION bin)


