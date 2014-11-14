# Install script for directory: /Users/petersson1/src/visit/visit2.8.1/src/databases/sw4img

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/usr/local")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "Release")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/2.8.1/darwin-x86_64/plugins/databases" TYPE SHARED_LIBRARY FILES "/Users/petersson1/src/visit/visit2.8.1/src/plugins/databases/libIsw4imgDatabase.dylib")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.8.1/darwin-x86_64/plugins/databases/libIsw4imgDatabase.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.8.1/darwin-x86_64/plugins/databases/libIsw4imgDatabase.dylib")
    EXECUTE_PROCESS(COMMAND "/usr/bin/install_name_tool"
      -id "libIsw4imgDatabase.dylib"
      -change "/Users/petersson1/src/visit/visit2.8.1/src/lib/libvisitcommon.dylib" "libvisitcommon.dylib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.8.1/darwin-x86_64/plugins/databases/libIsw4imgDatabase.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/petersson1/src/visit/visit2.8.1/src/lib/."
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.8.1/darwin-x86_64/plugins/databases/libIsw4imgDatabase.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/petersson1/src/visit/visit2.8.1/visit/vtk/6.1.0/i386-apple-darwin13_clang/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.8.1/darwin-x86_64/plugins/databases/libIsw4imgDatabase.dylib")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.8.1/darwin-x86_64/plugins/databases/libIsw4imgDatabase.dylib")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/2.8.1/darwin-x86_64/plugins/databases" TYPE SHARED_LIBRARY FILES "/Users/petersson1/src/visit/visit2.8.1/src/plugins/databases/libMsw4imgDatabase.dylib")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.8.1/darwin-x86_64/plugins/databases/libMsw4imgDatabase.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.8.1/darwin-x86_64/plugins/databases/libMsw4imgDatabase.dylib")
    EXECUTE_PROCESS(COMMAND "/usr/bin/install_name_tool"
      -id "libMsw4imgDatabase.dylib"
      -change "/Users/petersson1/src/visit/visit2.8.1/src/lib/libavtdatabase_ser.dylib" "libavtdatabase_ser.dylib"
      -change "/Users/petersson1/src/visit/visit2.8.1/src/lib/libavtdbatts.dylib" "libavtdbatts.dylib"
      -change "/Users/petersson1/src/visit/visit2.8.1/src/lib/libavtmath.dylib" "libavtmath.dylib"
      -change "/Users/petersson1/src/visit/visit2.8.1/src/lib/libavtmir_ser.dylib" "libavtmir_ser.dylib"
      -change "/Users/petersson1/src/visit/visit2.8.1/src/lib/libavtpipeline_ser.dylib" "libavtpipeline_ser.dylib"
      -change "/Users/petersson1/src/visit/visit2.8.1/src/lib/liblightweight_visit_vtk.dylib" "liblightweight_visit_vtk.dylib"
      -change "/Users/petersson1/src/visit/visit2.8.1/src/lib/libtess2.dylib" "libtess2.dylib"
      -change "/Users/petersson1/src/visit/visit2.8.1/src/lib/libvisit_verdict.dylib" "libvisit_verdict.dylib"
      -change "/Users/petersson1/src/visit/visit2.8.1/src/lib/libvisit_vtk.dylib" "libvisit_vtk.dylib"
      -change "/Users/petersson1/src/visit/visit2.8.1/src/lib/libvisitcommon.dylib" "libvisitcommon.dylib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.8.1/darwin-x86_64/plugins/databases/libMsw4imgDatabase.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/petersson1/src/visit/visit2.8.1/src/lib/."
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.8.1/darwin-x86_64/plugins/databases/libMsw4imgDatabase.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/petersson1/src/visit/visit2.8.1/visit/vtk/6.1.0/i386-apple-darwin13_clang/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.8.1/darwin-x86_64/plugins/databases/libMsw4imgDatabase.dylib")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.8.1/darwin-x86_64/plugins/databases/libMsw4imgDatabase.dylib")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

IF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  FILE(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/2.8.1/darwin-x86_64/plugins/databases" TYPE SHARED_LIBRARY FILES "/Users/petersson1/src/visit/visit2.8.1/src/plugins/databases/libEsw4imgDatabase_ser.dylib")
  IF(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.8.1/darwin-x86_64/plugins/databases/libEsw4imgDatabase_ser.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.8.1/darwin-x86_64/plugins/databases/libEsw4imgDatabase_ser.dylib")
    EXECUTE_PROCESS(COMMAND "/usr/bin/install_name_tool"
      -id "libEsw4imgDatabase_ser.dylib"
      -change "/Users/petersson1/src/visit/visit2.8.1/src/lib/libavtdatabase_ser.dylib" "libavtdatabase_ser.dylib"
      -change "/Users/petersson1/src/visit/visit2.8.1/src/lib/libavtdbatts.dylib" "libavtdbatts.dylib"
      -change "/Users/petersson1/src/visit/visit2.8.1/src/lib/libavtmath.dylib" "libavtmath.dylib"
      -change "/Users/petersson1/src/visit/visit2.8.1/src/lib/libavtmir_ser.dylib" "libavtmir_ser.dylib"
      -change "/Users/petersson1/src/visit/visit2.8.1/src/lib/libavtpipeline_ser.dylib" "libavtpipeline_ser.dylib"
      -change "/Users/petersson1/src/visit/visit2.8.1/src/lib/liblightweight_visit_vtk.dylib" "liblightweight_visit_vtk.dylib"
      -change "/Users/petersson1/src/visit/visit2.8.1/src/lib/libtess2.dylib" "libtess2.dylib"
      -change "/Users/petersson1/src/visit/visit2.8.1/src/lib/libvisit_verdict.dylib" "libvisit_verdict.dylib"
      -change "/Users/petersson1/src/visit/visit2.8.1/src/lib/libvisit_vtk.dylib" "libvisit_vtk.dylib"
      -change "/Users/petersson1/src/visit/visit2.8.1/src/lib/libvisitcommon.dylib" "libvisitcommon.dylib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.8.1/darwin-x86_64/plugins/databases/libEsw4imgDatabase_ser.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/petersson1/src/visit/visit2.8.1/src/lib/."
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.8.1/darwin-x86_64/plugins/databases/libEsw4imgDatabase_ser.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/petersson1/src/visit/visit2.8.1/visit/vtk/6.1.0/i386-apple-darwin13_clang/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.8.1/darwin-x86_64/plugins/databases/libEsw4imgDatabase_ser.dylib")
    IF(CMAKE_INSTALL_DO_STRIP)
      EXECUTE_PROCESS(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/2.8.1/darwin-x86_64/plugins/databases/libEsw4imgDatabase_ser.dylib")
    ENDIF(CMAKE_INSTALL_DO_STRIP)
  ENDIF()
ENDIF(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")

