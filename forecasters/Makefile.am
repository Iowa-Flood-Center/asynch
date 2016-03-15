noinst_LIBRARIES = libforecaster.a
libforecaster_a_SOURCES = forecaster_methods.c

bin_PROGRAMS = forecaster_maps
asynch_SOURCES = forecaster_maps.c
asynch_LDADD = libforecaster.a -lm $(POSTGRESQL_LIBS) $(SSH2_LIBS)

bin_PROGRAMS = forecaster_maps_end
forecaster_maps_end_SOURCES = forecaster_maps_end.c
forecaster_maps_end_LDADD = libforecaster.a -lm $(POSTGRESQL_LIBS) $(SSH2_LIBS)

bin_PROGRAMS = asynchpersis
asynchpersis_SOURCES = asynchpersis.c
asynchpersis_LDADD = libforecaster.a -lm $(POSTGRESQL_LIBS) $(SSH2_LIBS)

bin_PROGRAMS = asynchpersis_end
asynchpersis_end_SOURCES = asynchpersis_end.c
asynchpersis_end_LDADD = libforecaster.a -lm $(POSTGRESQL_LIBS) $(SSH2_LIBS)

AM_LDFLAGS = $(POSTGRESQL_LDFLAGS) $(SSH2_LDFLAGS)
AM_CFLAGS = $(POSTGRESQL_CPPFLAGS) $(SSH2_CPPFLAGS)
