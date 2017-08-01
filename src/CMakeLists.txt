file(GLOB_RECURSE sources  *.cpp *.h)
add_library(VinaLC ${sources})
target_link_libraries (VinaLC ${ZLIB_LIBRARIES})
install(TARGETS VinaLC DESTINATION lib)

file(GLOB_RECURSE header *.h)
install(FILES ${header} DESTINATION include)
