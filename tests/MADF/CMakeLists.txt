project(testMADF LANGUAGES CXX)

add_executable(testMADF testMADF.cpp)

find_package(ANASIR REQUIRED)

target_link_libraries(testMADF PRIVATE ANASIR::anasirlib)