project(testCF LANGUAGES CXX)

add_executable(main testCF.cpp)

find_package(ANASIR REQUIRED)

target_link_libraries(main PRIVATE ANASIR::anasirlib)