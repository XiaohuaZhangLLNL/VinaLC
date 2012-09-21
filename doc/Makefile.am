include $(top_srcdir)/doc/local/sources.mk
include $(top_srcdir)/doc/pics/sources.mk

SUBDIRS = local pics

docdir = @docdir@

all-local: html-local

html-local: html/index.html

if GOT_DOT
html/index.html: $(top_builddir)/doc/Doxyfile \
	$(doc_local_sources)
	doxygen $< ; \
	if test -d html; then \
        cp $(top_srcdir)/doc/pics/llnl_logo.gif html/; \
    fi
else
html/index.html: $(top_builddir)/doc/Doxyfile \
        $(doc_local_sources)
	(cat Doxyfile ; echo "HAVE_DOT=NO") | doxygen -;
	if test -d html; then \
        cp $(top_srcdir)/doc/pics/llnl_logo.gif html/; \
    fi
endif

disthtml: html-local

install-data-local: html-local
	$(mkinstalldirs) $(DESTDIR)$(docdir)
	if test -d html; then \
		$(mkinstalldirs) $(DESTDIR)$(docdir)/html; \
		for f in html/*; do \
			$(INSTALL_DATA) $$f $(DESTDIR)$(docdir)/html; \
		done \
	fi

uninstall-local:
	rm -rf $(DESTDIR)$(docdir)

clean-local:
	rm -rf html *.log *.tag
