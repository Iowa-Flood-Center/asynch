noinst_LIBRARIES = libasynch.a
libasynch_a_SOURCES = rkmethods.c problems.c mathmethods.c riversys.c sort.c comm.c system.c processdata.c partition.c definetype.c misc.c rainfall.c solvers.c io.c forcings.c compression.c date_manip.c asynch_interface.c modeloutputs.c data_types.c

bin_PROGRAMS = asynch
asynch_SOURCES = optparse.c asynchdist.c
asynch_LDADD = libasynch.a $(HDF5_LIBS) $(POSTGRESQL_LIBS)
asynch_LDFLAGS = $(HDF5_LDFLAGS) $(POSTGRESQL_LDFLAGS)

include_HEADERS = rkmethods.h problems.h mathmethods.h riversys.h sort.h comm.h system.h processdata.h partition.h definetype.h misc.h rainfall.h solvers.h io.h forcings.h compression.h date_manip.h asynch_interface.h modeloutputs.h data_types.h structs.h optparse.h libpq-fwd.h

AM_CFLAGS = $(HDF5_CFLAGS) $(POSTGRESQL_CPPFLAGS)
