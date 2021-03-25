file(GLOB_RECURSE sources  *.cpp *.h)
add_library(VinaLC ${sources})
target_link_libraries (VinaLC ${ZLIB_LIBRARIES})
install(TARGETS VinaLC DESTINATION lib)

file(GLOB_RECURSE header *.h)
install(FILES ${header} DESTINATION include)

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

#add_executable(test test.cpp)
target_link_libraries(VinaLC Threads::Threads)
# mpiBparser.cpp
add_executable(vinalc vinaBMPI.cpp dockBMPI.cpp mainProcedure.cpp)
target_link_libraries(vinalc VinaLC ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES})
set_target_properties(vinalc PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install(TARGETS vinalc DESTINATION bin)

