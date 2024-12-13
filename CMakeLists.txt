cmake_minimum_required(VERSION 3.15)

project(ANASIR VERSION 1.0.0)

add_library(anasirlib)
add_library(ANASIR::anasirlib ALIAS anasirlib)

option(BUILD_TEST "Build Test Files" OFF)

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/cmake/anasir-config.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/ANASIR")

if(MSVC)
  target_compile_options(anasirlib PRIVATE /W4 /WX)
else()
  target_compile_options(anasirlib PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()

add_subdirectory(core)

if(BUILD_TEST)
  add_subdirectory(tests)
endif()