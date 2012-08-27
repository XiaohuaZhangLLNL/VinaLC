if USE_MPI
    INCLUDES = -I. -I$(top_srcdir)/src -I/@boost_dir@/include -I/@mpi_path@/include
else
    INCLUDES = -I. -I$(top_srcdir)/src -I/@boost_dir@/include
endif

if USE_MPI
    bin_PROGRAMS = vinalc
else
endif

if USE_MPI
    vinalc_SOURCES = main.cpp docking.cpp mpiparser.cpp
    vinalc_LDFLAGS = -Wl,-rpath,@mpi_path@/lib/:@boost_dir@/lib/
    vinalc_LDADD = $(top_srcdir)/src/.libs/libVina.a -L/@boost_dir@/lib/ -lboost_system -lboost_thread -lboost_serialization -lboost_filesystem -lboost_program_options -L/@mpi_path@/lib/ -lmpi -lmpi_cxx -lz

else

endif
