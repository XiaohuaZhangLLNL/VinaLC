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
