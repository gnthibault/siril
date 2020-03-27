## Process this file with automake to produce Makefile.in

export CC CFLAGS CPPFLAGS

SUBDIRS = deps \
	po \
	src \
	pixmaps \
	platform-specific

# TODO add a conditionnal and a configure switch
# to be able to choose the language. Also, this should
# be a recussive make with a separate file list.
scriptsdir = $(pkgdatadir)/scripts
scripts_DATA = scripts/DSLR_preprocessing.ssf \
	3rdparty/scripts/en/DSLR_preprocessing_Cosmetic.ssf \
	3rdparty/scripts/en/DSLR_preprocessing_Drizzle_Cosmetic.ssf \
	3rdparty/scripts/en/DSLR_Preprocessing_Drizzle.ssf \
	3rdparty/scripts/en/DSLR_Preprocessing_NoDark.ssf \
	3rdparty/scripts/en/DSLR_Preprocessing_NoFlat_NoDark_NoBias.ssf \
	3rdparty/scripts/en/DSLR_Preprocessing_NoFlat_NoDark.ssf \
	3rdparty/scripts/en/DSLR_Preprocessing_NoFlat.ssf

doc_DATA = AUTHORS \
	LICENSE.md

man_MANS = siril.1

if !PLATFORM_WIN32
install-exec-hook:
	(cd $(DESTDIR)$(mandir)/man1; \
	rm -f siril-cli.1; \
	$(LN_S) siril.1 siril-cli.1)

uninstall-hook:
	(cd $(DESTDIR)$(mandir)/man1; rm -f siril-cli.1)
endif

ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST = config.rpath \
             git-version.h
AUTOMAKE_OPTIONS = foreign

generated_sources = \
	git-version.h

# Build git-version.h before anything in the subdirs as this is needed
# in the about and debug dialog (app/) and in the debug tool (tools/).
BUILT_SOURCES = $(generated_sources)

git-version.h: update-git-version-header
	@if test -e "$(top_srcdir)/.git"; then \
	  git_version="`git --git-dir=$(top_srcdir)/.git describe --always`"; \
	  git_version_abbrev="`git --git-dir=$(top_srcdir)/.git rev-parse --short HEAD`"; \
	  git_last_commit_year="`git --git-dir=$(top_srcdir)/.git log -n1 --reverse --pretty=%ci | cut -b 1-4`"; \
	elif test ! -f "$@"; then \
	  git_version="Unknown, shouldn't happen"; \
	  git_version_abbrev="$$git_version"; \
	  git_last_commit_timestamp=-1; \
	  git_last_commit_year="`date -u '+%Y'`"; \
	fi; \
	if test -n "$$git_version"; then \
	  echo "#ifndef __GIT_VERSION_H__" > "$@.tmp"; \
	  echo "#define __GIT_VERSION_H__" >> "$@.tmp"; \
	  echo "#define SIRIL_GIT_VERSION_ABBREV \"$$git_version_abbrev\"" >> "$@.tmp"; \
	  echo "#define SIRIL_GIT_LAST_COMMIT_YEAR \"$$git_last_commit_year\"" >> "$@.tmp"; \
	  echo "#endif /* __GIT_VERSION_H__ */" >> "$@.tmp"; \
	fi
	@if ( test -f "$@.tmp" && test -f "$@" && cmp "$@.tmp" "$@" > /dev/null ); then \
	  rm -f "$@.tmp"; \
	elif test -f "$@.tmp"; then \
	  mv "$@.tmp" "$@"; \
	  echo "  git HEAD changed: $@ regenerated"; \
	fi

.PHONY: update-git-version-header
