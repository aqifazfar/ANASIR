get_filename_component(ANASIR_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)

if(NOT TARGET ANASIR::anasirlib)
    include("${ANASIR_CMAKE_DIR}/anasir-targets.cmake")    
endif()

