noinst_LIBRARIES = libforecaster.a
libforecaster_a_SOURCES = forecaster_methods.c

bin_PROGRAMS = forecaster_maps forecaster_maps_end asynchpersis asynchpersis_end

forecaster_maps_SOURCES = forecaster_maps.c
forecaster_maps_LDADD = libforecaster.a $(top_builddir)/src/libasynch.a -lm $(POSTGRESQL_LIBS) $(HDF5_LIBS) $(SSH2_LIBS)

forecaster_maps_end_SOURCES = forecaster_maps_end.c
forecaster_maps_end_LDADD = libforecaster.a $(top_builddir)/src/libasynch.a -lm $(POSTGRESQL_LIBS) $(HDF5_LIBS) $(SSH2_LIBS)

asynchpersis_SOURCES = asynchpersis.c
asynchpersis_LDADD = libforecaster.a $(top_builddir)/src/libasynch.a  -lm $(POSTGRESQL_LIBS) $(HDF5_LIBS) $(SSH2_LIBS)

asynchpersis_end_SOURCES = asynchpersis_end.c
asynchpersis_end_LDADD = libforecaster.a $(top_builddir)/src/libasynch.a -lm $(POSTGRESQL_LIBS) $(HDF5_LIBS) $(SSH2_LIBS)

AM_LDFLAGS = $(POSTGRESQL_LDFLAGS) $(HDF5_LDFLAGS) $(SSH2_LDFLAGS)
AM_CFLAGS = -I$(srcdir)/../src $(POSTGRESQL_CPPFLAGS) $(HDF5_CFLAGS) $(SSH2_CPPFLAGS)
