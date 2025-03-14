set_target_properties(anasirlib PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(anasirlib PUBLIC  
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>  
                           $<INSTALL_INTERFACE:include>)

target_sources(anasirlib PRIVATE 
               include/ANASIR/Complementary_Filter.h 
               include/ANASIR/Extended_Kalman_Filter.h 
               include/ANASIR/Madgwick_Filter.h 
               include/ANASIR/Mahony_Filter.h 
               include/ANASIR/helpers.h
               include/ANASIR/integrator.h)

install(FILES
        include/ANASIR/Complementary_Filter.h
        include/ANASIR/Extended_Kalman_Filter.h
        include/ANASIR/Madgwick_Filter.h
        include/ANASIR/Mahony_Filter.h
        include/ANASIR/helpers.h
        include/ANASIR/integrator.h
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ANASIR")


include(GNUInstallDirs)

install(TARGETS anasirlib
	EXPORT anasir-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(EXPORT anasir-targets
        FILE anasir-targets.cmake
        NAMESPACE ANASIR::
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/ANASIR")


