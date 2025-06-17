AUTOMAKE_OPTIONS = subdir-objects

SUBDIRS = src tests tools
ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST = LICENSE README.md examples

dist-hook:
	rm -f $(distdir)/examples/.gitignore
#	mkdir $(distdir)/docs
#	cp -R $(srcdir)/docs/.build/html $(distdir)/docs
