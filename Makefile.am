
if DOC_ON
DOCDIR = doc
endif

AUTOMAKE_OPTIONS = foreign
#ACLOCAL_AMFLAGS = -I m4
SUBDIRS = src apps tests $(DOCDIR) 

# tools
include_HEADERS = config.h

EXTRA_DIST = Doxyfile autogen.sh

if DOC_ON
docs:
	$(MAKE) -C doc disthtml
endif

MAINTAINERCLEANFILES = configure aclocal.m4 acinclude.m4 autom4te.cache \
        Makefile.in config.h.in
