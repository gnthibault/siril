cmake_minimum_required(VERSION 3.5)

#################
#  Generic Part #
#################
project(siril)
set(OUTPUT_NAME siril)

###############
#  Src files  #
###############
set(core_sources
  core/siril.c
  core/command.c
  core/undo.c
  core/utils.c
  core/processing.c
  core/pipe.c
  core/initfile.c
  core/siril_update.c
  core/sequence_filtering.c
  core/preprocess.c
  core/OS_utils.c
  core/siril_log.c
)
include_directories(core)

set(io_sources
  io/conversion.c
  io/ser.c
  io/films.c
  io/image_formats_libraries.c
  io/image_format_fits.c
  io/sequence.c
  io/seqfile.c
  io/single_image.c
  io/mp4_output.c
  io/sequence_export.c
  io/avi_pipp/pipp_avi_write.cpp
  io/avi_pipp/pipp_avi_write_dib.cpp
  io/avi_pipp/pipp_buffer.cpp
  io/avi_pipp/avi_writer.cpp
)
include_directories(
  io
  io/avi_pipp/
)

set(gui_sources
  gui/callbacks.c
  gui/histogram.c
  gui/save_dialog.c
  gui/open_dialog.c
  gui/message_dialog.c
  gui/PSF_list.c
  gui/sequence_list.c
  gui/registration_preview.c
  gui/statistics_list.c
  gui/script_menu.c
  gui/plot.c
  gui/asinh.c
  gui/saturation.c
  gui/wavelets.c
  gui/gnuplot_i/gnuplot_i.c
  gui/progress_and_log.c
  gui/about_dialog.c
  gui/photometric_cc.c 
)
include_directories(
  gui
  gui/gnuplot_i
)

set(registration_sources
  registration/registration.c
  registration/registration.h
  registration/global.c
  registration/comet.c
  registration/matching/match.c
  registration/matching/atpmatch.c
  registration/matching/misc.c
  registration/matching/project_coords.c
  registration/matching/project_coords.h
  registration/matching/degtorad.h
  registration/matching/apply_match.c
)
include_directories(
  registration
  registration/matching
)

set(stacking_sources
  stacking/stacking.c
  stacking/sum.c
  stacking/normalization.c 
  stacking/median_and_mean.c
  stacking/upscaling.c
)
include_directories(stacking)

set(algos_sources
  algos/statistics.c
  algos/sorting.c
  algos/background_extraction.c
  algos/quality.c
  algos/rgradient.c
  algos/fft.c
  algos/colors.c
  algos/geometry.c
  algos/demosaicing.c
  algos/pave.c
  algos/transform.c
  algos/io_wave.c
  algos/reconstr.c
  algos/PSF.c
  algos/star_finder.c
  algos/cosmetic_correction.c
  algos/quantize.c
  algos/photometry.c
  algos/plateSolver.c
)
include_directories(algos)

set(compositing_sources
  compositing/compositing.c
  compositing/filters.c
  compositing/align_rgb.c
)
include_directories(compositing)

set(opencv_sources
  opencv/opencv.cpp
  opencv/ecc/ecc.cpp
  opencv/findHomography/fundam.cpp
  opencv/findHomography/calibration.cpp
  opencv/findHomography/modelest.cpp
)
include_directories(
  opencv
  opencv/findHomography
)

set(sources
  main.c
  ${core_sources}
  ${io_sources}
  ${gui_sources}
  ${registration_sources}
  ${stacking_sources}
  ${algos_sources}
  ${compositing_sources}
  ${opencv_sources}
)

#siril_LDADD =  $(top_builddir)/deps/kplot/libkplot.a

#pkgdata_DATA = gui/siril3.glade \
#	gui/gtk_old.css \
#	gui/gtk.css

################
#  Executable  #
################

add_executable(${OUTPUT_NAME} SHARED ${sources})
target_link_libraries(${OUTPUT_NAME}
  ${boost_deps}
  ${tiff_deps}
  $(jpeg_deps)
  $(png_deps)
  $(ffmpeg_deps)
  ${ffms_deps}
  ${curl_deps}
  ${gtk_deps}
  ${libraw_deps}
  $(cfitsio_deps)
  $(fftw_deps)
  $(gsl_deps)
  ${openmp_deps}
)
set_property(TARGET ${OUTPUT_NAME}
  PROPERTY OUTPUT_NAME ${output_name})
install(TARGETS ${OUTPUT_NAME}
  LIBRARY DESTINATION ${LIBRARY_OUTPUT_DIRECTORY}
)
