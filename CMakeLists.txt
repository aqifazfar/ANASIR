add_library(anasirlib)
add_library(BUMANTARASDK::anasirlib ALIAS anasirlib)

if(MSVC)
  target_compile_options(anasirlib PRIVATE /W4 /WX)
else()
  target_compile_options(anasirlib PRIVATE -Wall -Wextra -Wpedantic -Werror)
endif()
        
add_subdirectory(core)