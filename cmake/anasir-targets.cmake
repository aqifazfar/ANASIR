add_library(ANASIR::anasirlib IMPORTED)

set_target_properties(ANASIR::anasirlib
                      PROPERTIES
                      INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include")