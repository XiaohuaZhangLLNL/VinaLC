
add_executable(vinalc vinaBMPI.cpp dockBMPI.cpp mpiBparser.cpp mainProcedure.cpp)
target_link_libraries(vinalc VinaLC ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES})
set_target_properties(vinalc PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
install(TARGETS vinalc DESTINATION bin)

