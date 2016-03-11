noinst_LIBRARIES = libasynch.a
libasynch_a_SOURCES = rkmethods.c problems.c mathmethods.c riversys.c sort.c comm.c system.c processdata.c partition.c definetype.c misc.c rainfall.c solvers.c io.c 	forcings.c compression.c date_manip.c asynch_interface.c modeloutputs.c data_types.c

bin_PROGRAMS = asynch
asynch_SOURCES = asynchdist.c
asynch_LDADD = libasynch.a $(POSTGRESQL_LIBS)
asynch_LDFLAGS = $(POSTGRESQL_LDFLAGS)

AM_CFLAGS = $(POSTGRESQL_CPPFLAGS)
