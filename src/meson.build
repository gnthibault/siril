################################################################################
# Meson siril file
################################################################################

################################################################################
# Generate configuration files

configure_file(output : 'git-version.h',
               configuration : git_data)

configure_file(output : 'config.h',
               configuration : conf_data)

################################################################################
# List sources
src_files = [ 'algos/background_extraction.c',
              'algos/colors.c',
              'algos/demosaicing.c',
              'algos/demosaicing_rtp.cpp',
              'algos/geometry.c',
              'algos/io_wave.c',
              'algos/median_fast.c',
              'algos/noise.c',
              'algos/pave.c',
              'algos/photometry.c',
              'algos/plateSolver.c',
              'algos/PSF.c',
              'algos/quality.c',
              'algos/quality_float.c',
              'algos/quantize.c',
              'algos/reconstr.c',
              'algos/sorting.c',
              'algos/star_finder.c',
              'algos/statistics.c',
              'algos/statistics_float.c',
              'algos/transform.c',

              'compositing/align_rgb.c',
              'compositing/compositing.c',
              'compositing/filters.c',

              'core/arithm.c', 
              'core/siril_app_dirs.c', 
              'core/preprocess.c', 
              'core/siril_cmd_help.c', 
              'core/command.c',
              'core/processing.c',         
              'core/siril_language.c',
              'core/command_line_processor.c', 
              'core/sequence_filtering.c', 
              'core/siril_log.c',
              'core/exif.cpp',             
              'core/signals.c',            
              'core/siril_update.c',
              'core/initfile.c',                
              'core/siril_actions.c',       
              'core/undo.c',
              'core/OS_utils.c',                
              'core/siril_app_dirs.c',      
              'core/utils.c',
              'core/pipe.c',                    
              'core/siril.c', 

              'filters/asinh.c',
              'filters/banding.c',
              'filters/clahe.c',
              'filters/cosmetic_correction.c',
              'filters/deconv.c',
              'filters/fft.c',
              'filters/median.c',
              'filters/rgradient.c',
              'filters/saturation.c',
              'filters/scnr.c',
              'filters/wavelets.c',

              'io/avi_pipp/pipp_avi_write.cpp',
              'io/avi_pipp/pipp_avi_write.h',
              'io/avi_pipp/pipp_avi_write_dib.cpp',
              'io/avi_pipp/pipp_avi_write_dib.h',
              'io/avi_pipp/pipp_utf8.h',
              'io/avi_pipp/pipp_buffer.cpp',
              'io/avi_pipp/pipp_buffer.h',
              'io/avi_pipp/pipp_video_write.h',
              'io/avi_pipp/avi_writer.cpp',
              'io/avi_pipp/avi_writer.h',
              'io/conversion.c',
              'io/films.c',
              'io/image_format_fits.c',
              'io/image_formats_internal.c',
              'io/image_formats_libraries.c',
              'io/mp4_output.c',
              'io/seqfile.c',
              'io/sequence.c',
              'io/sequence_export.c',
              'io/ser.c',
              'io/single_image.c',

              'opencv/opencv.cpp',
              'opencv/opencv.h',
              'opencv/ecc/ecc.cpp',
              'opencv/ecc/ecc.h',
              'opencv/findHomography/fundam.cpp',
              'opencv/findHomography/calibration.cpp',
              'opencv/findHomography/modelest.cpp',
              'opencv/findHomography/_modelest.h',
              'opencv/findHomography/calib3d.hpp',
              'opencv/findHomography/precomp.hpp',

              'gui/about_dialog.c',
              'gui/callbacks.c',
              'gui/conversion.c',
              'gui/dialog_preview.c',
              'gui/dialogs.c',
              'gui/gnuplot_i/gnuplot_i.c',
              'gui/histogram.c',
              'gui/image_display.c',
              'gui/image_interactions.c',
              'gui/linear_match.c',
              'gui/message_dialog.c',
              'gui/open_dialog.c',
              'gui/photometric_cc.c',
              'gui/plot.c',
              'gui/preferences.c',
              'gui/progress_and_log.c',
              'gui/PSF_list.c',
              'gui/registration_preview.c',
              'gui/save_dialog.c',
              'gui/script_menu.c',
              'gui/sequence_list.c',
              'gui/siril_preview.c',
              'gui/statistics_list.c',

              'registration/comet.c',
              'registration/global.c',
              'registration/matching/match.c',
              'registration/matching/atpmatch.c',
              'registration/matching/misc.c',
              'registration/matching/project_coords.c',
              'registration/matching/apply_match.c',
              'registration/registration.c',

              'stacking/median_and_mean.c',
              'stacking/rejection_float.c',
              'stacking/normalization.c',
              'stacking/siril_fit_linear.c',
              'stacking/stacking.c',
              'stacking/sum.c',
              'stacking/upscaling.c',

              'rt/boxblur.cc',
              'rt/deconvolution.cc',
              'rt/gauss.cc',
              'rt/rt_algo.cc']

################################################################################
# Build gui executable

siril_gui = executable('siril','main.c', src_files,
                      dependencies : [m_dep, glib_dep, gtk_dep, gsl_dep, opencv_dep, rtprocess_dep, exiv_dep, libconfig_dep,fftw3f_dep, cfitsio_dep, libraw_dep, libtiff_dep, libjpeg_dep, libpng_dep, libcurl_dep, libheif_dep, ffms2_dep, libcurl_dep, openmp_dep, kplot_dep],
                      link_args : '-rdynamic',
                      c_args : ['-DHAVE_CONFIG_H'],
                      cpp_args : [siril_cpp_flag]
                      , install : true
                      )

################################################################################
# Build cli executable

siril_cli = executable('siril-cli','main-cli.c', src_files,
                      dependencies : [m_dep, glib_dep, gtk_dep, gsl_dep, opencv_dep, rtprocess_dep, exiv_dep, libconfig_dep,fftw3f_dep, cfitsio_dep, libraw_dep, libtiff_dep, libjpeg_dep, libpng_dep, libcurl_dep, libheif_dep, ffms2_dep, libcurl_dep, openmp_dep, kplot_dep],
                      link_args : '-rdynamic',
                      c_args : ['-DHAVE_CONFIG_H'],
                      cpp_args : [siril_cpp_flag]
                      , install : true
                      )

################################################################################
# Install GUI files

install_data(['gui/siril.css', 'gui/siril_old.css', 'gui/siril3.glade', 'gui/siril-shortcuts.ui'],
             install_dir : siril_datadir)