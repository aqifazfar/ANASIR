file(REMOVE_RECURSE
  "libanasirlib.a"
  "libanasirlib.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/anasirlib.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
