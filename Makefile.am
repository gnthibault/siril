## Process this file with automake to produce Makefile.in

export CC CFLAGS CPPFLAGS

SUBDIRS = deps/kplot \
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

ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST = config.rpath
AUTOMAKE_OPTIONS = foreign
