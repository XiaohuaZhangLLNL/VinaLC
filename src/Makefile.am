CFLAGS = -Wall 

if USE_MPI
INCLUDES = -I$(top_srcdir)/src -I@boost_dir@/include -I/@mpi_path@/include
else
INCLUDES = -I$(top_srcdir)/src -I@boost_dir@/include 
endif

BOOST_LIB = -L/@boost_dir@/lib/ -lboost_system -lboost_thread -lboost_serialization -lboost_filesystem -lboost_program_options -lz 

if USE_MPI
MPI_LIB= -L/@mpi_path@/lib/ -lmpi -lmpi_cxx
endif

lib_LTLIBRARIES = libVina.la

libVina_la_SOURCES = \
cache.cpp \
coords.cpp \
current_weights.cpp \
everything.cpp \
grid.cpp \
gzstream.cpp \
manifold.cpp \
model.cpp \
monte_carlo.cpp \
mutate.cpp \
my_pid.cpp \
naive_non_cache.cpp \
non_cache.cpp \
parallel_mc.cpp \
parallel_progress.cpp \
parse_pdbqt.cpp \
pdb.cpp \
quasi_newton.cpp \
quaternion.cpp \
random.cpp \
ssd.cpp \
szv_grid.cpp \
terms.cpp \
weighted_terms.cpp 

libVinaincludedir = $(includedir)/

libVinainclude_HEADERS = \
array3d.h \
atom.h \
atom_base.h \
atom_constants.h \
atom_type.h \
bfgs.h \
brick.h \
cache.h \
common.h \
conf.h \
convert_substring.h \
coords.h \
curl.h \
current_weights.h \
everything.h \
file.h \
grid.h \
grid_dim.h \
gzstream.h \
igrid.h \
incrementable.h \
int_pow.h \
macros.h \
manifold.h \
matrix.h \
model.h \
monte_carlo.h \
mutate.h \
my_pid.h \
naive_non_cache.h \
non_cache.h \
parallel.h \
parallel_mc.h \
parallel_progress.h \
parse_error.h \
parse_pdbqt.h \
pdb.h \
precalculate.h \
quasi_newton.h \
quaternion.h \
random.h \
recent_history.h \
scoring_function.h \
ssd.h \
statistics.h \
szv_grid.h \
tee.h \
terms.h \
tokenize.h \
tree.h \
triangular_matrix_index.h \
weighted_terms.h 

