if USE_MPI
    INCLUDES = -I. -I$(top_srcdir)/src -I/@boost_dir@/include -I/@mpi_path@/include
else
    INCLUDES = -I. -I$(top_srcdir)/src -I/@boost_dir@/include
endif

if USE_MPI
    bin_PROGRAMS = vinalc
else
    bin_PROGRAMS = vina
endif

if USE_MPI
    vinalc_SOURCES = vinalc.cpp docking.cpp mpiparser.cpp mainProcdure.cpp
    vinalc_LDFLAGS = -Wl,-rpath,@mpi_path@/lib/:@boost_dir@/lib/
    vinalc_LDADD = $(top_srcdir)/src/.libs/libVinaLC.a -L/@boost_dir@/lib/ -lboost_system -lboost_thread -lboost_serialization -lboost_filesystem -lboost_program_options -L/@mpi_path@/lib/ -lmpi -lmpi_cxx -lz

else
    vina_SOURCES = vina.cpp mainProcdure.cpp
    vina_LDFLAGS = -Wl,-rpath,@boost_dir@/lib/
    vina_LDADD = $(top_srcdir)/src/.libs/libVina.a -L/@boost_dir@/lib/ -lboost_system -lboost_thread -lboost_serialization -lboost_filesystem -lboost_program_options -lz

endif
