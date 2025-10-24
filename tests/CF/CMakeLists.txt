project(testCF LANGUAGES CXX)

add_executable(testCF testCF.cpp)

find_package(ANASIR REQUIRED)

target_link_libraries(testCF PRIVATE ANASIR::anasirlib)