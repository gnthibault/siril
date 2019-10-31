## Process this file with automake to produce Makefile.in

bin_PROGRAMS = siril

AM_CPPFLAGS = \
	-DPACKAGE_DATA_DIR='"$(datadir)/${PACKAGE}"' \
	-DPACKAGE_DOC_DIR='"${datarootdir}/doc/${PACKAGE}"' \
	-DLOCALEDIR='"${localedir}"' \
	${OPENCV_CFLAGS} \
	$(FFMPEG_CFLAGS) \
	$(GTK_CFLAGS) \
	${GTK_MAC_CFLAGS} \
	$(LIBTIFF_CFLAGS) \
	$(LIBPNG_CFLAGS) \
	$(FFTW_CFLAGS) \
	$(CFITSIO_CFLAGS) \
	$(GSL_CFLAGS) \
	$(FFMS2_CFLAGS) \
	$(LIBCONFIG_CFLAGS) \
	-I$(top_builddir)/deps/kplot

siril_SOURCES = \
	main.c \
	algos/background_extraction.c \
	algos/background_extraction.h \
	algos/colors.c \
	algos/colors.h \
	algos/Def_Math.h \
	algos/Def_Wavelet.h \
	algos/Def_Mem.h \
	algos/demosaicing.c \
	algos/demosaicing.h \
	algos/geometry.c \
	algos/geometry.h \
	algos/io_wave.c \
	algos/noise.c \
	algos/noise.h \
	algos/pave.c \
	algos/photometry.c \
	algos/photometry.h \
	algos/plateSolver.c \
	algos/plateSolver.h \
	algos/PSF.c \
	algos/PSF.h \
	algos/quality.c \
	algos/quality.h \
	algos/quantize.c \
	algos/reconstr.c \
	algos/sorting.h \
	algos/sorting.c \
	algos/star_finder.c \
	algos/star_finder.h \
	algos/statistics.c \
	algos/transform.c \
	compositing/align_rgb.c \
	compositing/align_rgb.h \
	compositing/compositing.c \
	compositing/compositing.h \
	compositing/filters.c \
	compositing/filters.h \
	core/command.c \
	core/command.h \
	core/command_def.h \
	core/initfile.c \
	core/initfile.h \
	core/OS_utils.c \
	core/pipe.c \
	core/pipe.h \
	core/preprocess.c \
	core/preprocess.h \
	core/processing.c \
	core/proto.h \
	core/sequence_filtering.c \
	core/sequence_filtering.h \
	core/signals.c \
	core/signals.h \
	core/siril.c \
	core/siril.h \
	core/siril_log.c \
	core/siril_update.c \
	core/siril_update.h \
	core/undo.c \
	core/undo.h \
	core/utils.c \
	filters/asinh.c \
	filters/asinh.h \
	filters/banding.c \
	filters/banding.h \
	filters/clahe.c \
	filters/clahe.h \
	filters/cosmetic_correction.c \
	filters/cosmetic_correction.h \
	filters/deconv.c \
	filters/deconv.h \
	filters/fft.c \
	filters/fft.h \
	filters/median.c \
	filters/median.h \
	filters/rgradient.c \
	filters/rgradient.h \
	filters/saturation.c \
	filters/saturation.h \
	filters/scnr.c \
	filters/scnr.h \
	filters/wavelets.c \
	filters/wavelets.h \
	gui/about_dialog.c \
	gui/callbacks.c \
	gui/callbacks.h \
	gui/dialogs.c \
	gui/dialogs.h \
	gui/gnuplot_i/gnuplot_i.c \
	gui/gnuplot_i/gnuplot_i.h \
	gui/histogram.c \
	gui/histogram.h \
	gui/message_dialog.c \
	gui/message_dialog.h \
	gui/open_dialog.c \
	gui/open_dialog.h \
	gui/photometric_cc.c \
	gui/photometric_cc.h \
	gui/plot.c \
	gui/plot.h \
	gui/progress_and_log.c \
	gui/progress_and_log.h \
	gui/PSF_list.h \
	gui/PSF_list.c \
	gui/registration_preview.c \
	gui/save_dialog.c \
	gui/save_dialog.h \
	gui/script_menu.c \
	gui/script_menu.h \
	gui/sequence_list.c \
	gui/statistics_list.c \
	io/conversion.c \
	io/conversion.h \
	io/films.c \
	io/films.h \
	io/image_format_fits.c \
	io/image_formats_internal.c \
	io/image_formats_libraries.c \
	io/mp4_output.c \
	io/mp4_output.h \
	io/seqfile.c \
	io/sequence.c \
	io/sequence.h \
	io/sequence_export.c \
	io/ser.c \
	io/ser.h \
	io/single_image.c \
	io/single_image.h \
	registration/comet.c \
	registration/global.c \
	registration/matching/match.c \
	registration/matching/atpmatch.c \
	registration/matching/misc.c \
	registration/matching/project_coords.c \
	registration/matching/project_coords.h \
	registration/matching/degtorad.h \
	registration/matching/apply_match.c \
	registration/registration.c \
	registration/registration.h \
	stacking/median_and_mean.c \
	stacking/normalization.c \
	stacking/stacking.c \
	stacking/stacking.h \
	stacking/sum.c \
	stacking/sum.h \
	stacking/upscaling.c

if OS_WIN32
libdbghelp = -ldbghelp
endif

siril_LDADD = \
	$(top_builddir)/deps/kplot/libkplot.a \
	$(FFMPEG_LIBS) \
	$(GTK_LIBS) \
	${GTK_MAC_LIBS} \
	$(LIBRAW_LIBS) \
	$(LIBTIFF_LIBS) \
	$(LIBPNG_LIBS) \
	$(FFTW_LIBS) \
	$(CFITSIO_LIBS) \
	$(GSL_LIBS) \
	$(FFMS2_LIBS) \
	$(LIBCONFIG_LIBS) \
	$(CURL_LIBS) \
	$(EXTRA_LIBS) \
	$(libdbghelp)

siril_SOURCES += \
	io/avi_pipp/pipp_avi_write.cpp \
	io/avi_pipp/pipp_avi_write.h \
	io/avi_pipp/pipp_avi_write_dib.cpp \
	io/avi_pipp/pipp_avi_write_dib.h \
	io/avi_pipp/pipp_utf8.h \
	io/avi_pipp/pipp_buffer.cpp \
	io/avi_pipp/pipp_buffer.h \
	io/avi_pipp/pipp_video_write.h \
	io/avi_pipp/avi_writer.cpp \
	io/avi_pipp/avi_writer.h

siril_SOURCES += \
	opencv/opencv.cpp \
	opencv/opencv.h \
	opencv/ecc/ecc.cpp \
	opencv/ecc/ecc.h

siril_SOURCES += \
	opencv/findHomography/fundam.cpp \
	opencv/findHomography/calibration.cpp \
	opencv/findHomography/modelest.cpp \
	opencv/findHomography/_modelest.h \
	opencv/findHomography/calib3d.hpp \
	opencv/findHomography/precomp.hpp

siril_LDADD += $(OPENCV_LIBS)

pkgdata_DATA = \
	gui/siril3.glade \
	gui/gtk_old.css \
	gui/gtk.css

if OS_WIN32
AM_LDFLAGS += -Wl,--export-all-symbols
if !SIRIL_DEBUG
AM_LDFLAGS += -mwindows
endif #!SIRIL_DEBUG
if HAVE_WINDRES
$(top_builddir)/platform-specific/windows/siril.res.o: $(top_builddir)/platform-specific/windows/siril.rc
	$(WINDRES) $< $@
siril_LDADD += $(top_builddir)/platform-specific/windows/siril.res.o
endif #HAVE_WINDRES
else
if !OS_OSX
AM_LDFLAGS += -rdynamic
endif #!OS_OSX
endif #OS_WIN32

siril_LINK = $(CXXLINK)
