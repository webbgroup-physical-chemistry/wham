AM_CPPFLAGS  = $(BOOST_CPPFLAGS)
AM_CPPFLAGS += $(HDF5_CFLAGS)
AM_LDFLAGS = 

bin_PROGRAMS = bw_wham pca_wham
bw_wham_SOURCES = bw_main.cpp bw_option_parser.cpp bw_option_parser.h bw_boltzmann_weight.cpp bw_boltzmann_weight.h bw_hdf5.cpp bw_hdf5.h bw_types.h 
bw_wham_CPPFLAGS  = $(BOOST_CPPFLAGS)
bw_wham_CPPFLAGS += $(HDF5_CFLAGS)
bw_wham_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS)
bw_wham_LDFLAGS += $(HDF5_LIBS) $(HDF5_CXXLIBS)

pca_wham_SOURCES = pca_main.cpp pca_option_parser.cpp pca_option_parser.h pca_wham.cpp pca_wham.h bw_types.h bw_hdf5.cpp bw_hdf5.h cpp_lapack.cpp cpp_lapack.h cpp_blas.cpp cpp_blas.h
pca_wham_CPPFLAGS  = $(BOOST_CPPFLAGS)
pca_wham_CPPFLAGS += $(HDF5_CFLAGS)
pca_wham_LDFLAGS = $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS)
pca_wham_LDFLAGS += $(HDF5_LIBS) $(HDF5_CXXLIBS)
pca_wham_LDFLAGS += $(BLAS_LIBS)
