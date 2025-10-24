cmake_minimum_required(VERSION 3.15)

project(ANASIR VERSION 1.0.0)

add_library(anasirlib)
add_library(ANASIR::anasirlib ALIAS anasirlib)

option(BUILD_TESTS "Build Test Files" OFF)

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/cmake/anasir-config.cmake"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/ANASIR")


if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "Install path" FORCE)
endif()
    

if(MSVC)
  target_compile_options(anasirlib PRIVATE /W4 /WX)
else()
  target_compile_options(anasirlib PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()

find_package(Eigen3 REQUIRED)

target_link_libraries(anasirlib PRIVATE Eigen3::Eigen)

add_subdirectory(core)

if(BUILD_TESTS)
  add_subdirectory(tests)
endif()