# Install script for directory: /mnt/c/Users/aqifa/OneDrive/Desktop/personal/projek/ANASIR/core

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/ANASIR" TYPE FILE FILES
    "/mnt/c/Users/aqifa/OneDrive/Desktop/personal/projek/ANASIR/core/include/ANASIR/Complementary_Filter.h"
    "/mnt/c/Users/aqifa/OneDrive/Desktop/personal/projek/ANASIR/core/include/ANASIR/Extended_Kalman_Filter.h"
    "/mnt/c/Users/aqifa/OneDrive/Desktop/personal/projek/ANASIR/core/include/ANASIR/Madgwick_Filter.h"
    "/mnt/c/Users/aqifa/OneDrive/Desktop/personal/projek/ANASIR/core/include/ANASIR/Mahony_Filter.h"
    "/mnt/c/Users/aqifa/OneDrive/Desktop/personal/projek/ANASIR/core/include/ANASIR/helpers.h"
    "/mnt/c/Users/aqifa/OneDrive/Desktop/personal/projek/ANASIR/core/include/ANASIR/integrator.h"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/mnt/c/Users/aqifa/OneDrive/Desktop/personal/projek/ANASIR/build/libanasirlib.a")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/ANASIR/anasir-targets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/ANASIR/anasir-targets.cmake"
         "/mnt/c/Users/aqifa/OneDrive/Desktop/personal/projek/ANASIR/build/core/CMakeFiles/Export/804b6cb4178c448c8edbf0b312ecfdc6/anasir-targets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/ANASIR/anasir-targets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/ANASIR/anasir-targets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/ANASIR" TYPE FILE FILES "/mnt/c/Users/aqifa/OneDrive/Desktop/personal/projek/ANASIR/build/core/CMakeFiles/Export/804b6cb4178c448c8edbf0b312ecfdc6/anasir-targets.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^()$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/ANASIR" TYPE FILE FILES "/mnt/c/Users/aqifa/OneDrive/Desktop/personal/projek/ANASIR/build/core/CMakeFiles/Export/804b6cb4178c448c8edbf0b312ecfdc6/anasir-targets-noconfig.cmake")
  endif()
endif()

