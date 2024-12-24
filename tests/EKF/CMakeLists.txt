project(testEKF LANGUAGES CXX)

add_executable(testEKF testEKF.cpp)

find_package(ANASIR REQUIRED)

target_link_libraries(testEKF PRIVATE ANASIR::anasirlib)