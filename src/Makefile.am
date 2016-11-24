## Process this file with automake to produce Makefile.in

bin_PROGRAMS = siril

AM_CPPFLAGS = -DPACKAGE_DATA_DIR='"$(datadir)/${PACKAGE}"' -DPACKAGE_DOC_DIR='"${datarootdir}/doc/${PACKAGE}"' -DLOCALEDIR='"${localedir}"' $(FFMPEG_CFLAGS) $(GTK_CFLAGS) ${GTK_MAC_CFLAGS} $(LIBTIFF_CFLAGS) $(LIBPNG_CFLAGS) $(FFTW_CFLAGS) $(CFITSIO_CFLAGS) $(GSL_CFLAGS) $(FFMS2_CFLAGS) $(LIBCONFIGXX_CFLAGS) -I ../deps/kplot

siril_SOURCES = \
	main.c \
	core/siril.c core/siril.h core/command.c core/command.h core/proto.h core/undo.c core/undo.h core/utils.c core/processing.c \
	core/initfile.c core/initfile.h \
	io/conversion.c io/conversion.h io/ser.c io/ser.h io/films.c io/films.h \
	io/image_formats_libraries.c io/image_formats_internal.c io/image_format_fits.c \
	io/sequence.c io/seqfile.c io/single_image.c io/single_image.h \
	io/mp4_output.h io/mp4_output.c \
	gui/callbacks.c gui/callbacks.h gui/histogram.c gui/histogram.h \
	gui/PSF_list.h gui/PSF_list.c gui/sequence_list.c gui/registration_preview.c \
	gui/statistics_list.c \
	gui/plot.c gui/plot.h \
	registration/registration.c registration/registration.h \
	registration/matching/match.c registration/matching/atpmatch.c registration/matching/misc.c \
	stacking/stacking.c stacking/stacking.h \
	algos/gradient.c algos/gradient.h algos/quality.c algos/quality.h \
	algos/statistics.c \
	algos/fft.c algos/fft.h \
	algos/colors.c algos/colors.h \
	algos/demosaicing.c algos/demosaicing.h \
	algos/pave.c algos/transform.c algos/io_wave.c algos/reconstr.c \
	algos/Def_Math.h algos/Def_Wavelet.h algos/Def_Mem.h \
	algos/PSF.c algos/PSF.h algos/star_finder.c algos/star_finder.h \
	algos/cosmetic_correction.c algos/cosmetic_correction.h \
	algos/quantize.c \
	compositing/compositing.c compositing/compositing.h compositing/filters.c compositing/filters.h compositing/align_rgb.c compositing/align_rgb.h
	
siril_LDADD = $(FFMPEG_LIBS) $(GTK_LIBS) ${GTK_MAC_LIBS} $(LIBRAW_LIBS) $(LIBTIFF_LIBS) $(LIBPNG_LIBS) $(FFTW_LIBS) $(CFITSIO_LIBS) $(GSL_LIBS) $(FFMS2_LIBS) $(LIBCONFIGXX_LIBS) ../deps/kplot/libkplot.a

siril_SOURCES += io/avi_pipp/pipp_avi_write.cpp io/avi_pipp/pipp_avi_write.h io/avi_pipp/pipp_avi_write_dib.cpp io/avi_pipp/pipp_avi_write_dib.h \
	io/avi_pipp/pipp_buffer.cpp io/avi_pipp/pipp_buffer.h io/avi_pipp/pipp_video_write.h io/avi_pipp/avi_writer.cpp io/avi_pipp/avi_writer.h

if HAVE_OPENCV
siril_SOURCES += opencv/opencv.cpp opencv/opencv.h opencv/ecc/ecc.cpp opencv/ecc/ecc.h
siril_LDADD += $(OPENCV_LIBS)
siril_LINK = $(CXXLINK)
else
siril_LINK = $(LINK)
endif
