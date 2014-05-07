AM_LDFLAGS = 

bin_PROGRAMS = wham
wham_SOURCES = main.cpp parse_options.cpp parse_options.h wham.cpp wham.h wham_befeus.cpp wham_befeus.h wham_fileutils.cpp wham_fileutils.h wham_functions.cpp wham_function.h wham_hdf5.cpp wham_hdf5.h wham_trig.cpp wham_trig.h wham_types.h
wham_CPPFLAGS  = $(BOOST_CPPFLAGS)
wham_CPPFLAGS += $(HDF5_CFLAGS)
wham_LDFLAGS  = $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) 
wham_LDFLAGS += $(HDF5_LIBS) $(HDF5_CXXLIBS)
wham_LDFLAGS += $(BLAS_LIBS)
