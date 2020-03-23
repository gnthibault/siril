# Trick to run the build rules of current Makefile before the subdir
# ones, so that we can move the files (see below) first.
SUBDIRS = . \
		  kplot

# This is a very ugly trick because kplot tree was just copied as is
# with a raw Makefile. So let's just copy the whole source before
# compiling if the build dir is different to the source dir.
all-local:
	if [ ! -f "$(top_builddir)/deps/kplot/Makefile" ]; then \
		cp -fr "$(top_srcdir)/deps/kplot" "$(top_builddir)/deps/kplot"; \
	fi
