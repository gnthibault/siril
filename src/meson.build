################################################################################
# Meson siril file
################################################################################

siril_inc_dir = include_directories('.')

################################################################################
# Generate configuration files

configure_file(output : 'git-version.h',
               configuration : git_data)

configure_file(output : 'config.h',
               configuration : conf_data)

################################################################################
# Generate Siril static library with all sources

src_files = [ 
  'algos/annotate.c',
  'algos/background_extraction.c',
  'algos/colors.c',
  'algos/demosaicing.c',
  'algos/demosaicing_rtp.cpp',
  'algos/extraction.c',
  'algos/fix_xtrans_af.c',
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
  'algos/search_objects.c',
  'algos/siril_wcs.c',
  
  'compositing/align_rgb.c',
  'compositing/compositing.c',
  'compositing/filters.c',
  
  'core/arithm.c',
  'core/atomic.cpp',
  'core/command.c',
  'core/command_line_processor.c',
  'core/exif.cpp',
  'core/icc_profile.c',
  'core/initfile.c',
  'core/OS_utils.c',
  'core/pipe.c',
  'core/preprocess.c',
  'core/processing.c',
  'core/sequence_filtering.c',
  'core/signals.c',
  'core/siril.c',
  'core/siril_actions.c',
  'core/siril_app_dirs.c',
  'core/siril_cmd_help.c',
  'core/siril_date.c',
  'core/siril_language.c',
  'core/siril_log.c',
  'core/siril_update.c',
  'core/siril_world_cs.c',
  'core/undo.c',
  'core/utils.c',
  
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
  'io/fits_sequence.c',
  'io/FITS_symlink.c',
  'io/image_format_fits.c',
  'io/image_formats_internal.c',
  'io/image_formats_libraries.c',
  'io/mp4_output.c',
  'io/seqfile.c',
  'io/sequence.c',
  'io/sequence_export.c',
  'io/seqwriter.c',

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
  'gui/fix_xtrans_af.c',
  'gui/gnuplot_i/gnuplot_i.c',
  'gui/histogram.c',
  'gui/image_display.c',
  'gui/image_interactions.c',
  'gui/linear_match.c',
  'gui/message_dialog.c',
  'gui/menu_gray_psf.c',
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
  'gui/single_image.c',
  'gui/siril_css.c',
  'gui/siril_intro.c',
  'gui/siril_preview.c',
  'gui/siril-window.c',
  'gui/statistics_list.c',
  'gui/utils.c',
  
  'registration/3stars.c',
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
  'stacking/stackminmax.c',
  'stacking/sum.c',
  'stacking/upscaling.c',
  
  'rt/boxblur.cc',
  'rt/deconvolution.cc',
  'rt/gauss.cc',
  'rt/rt_algo.cc'
]

siril_lib = static_library('siril',
                           src_files,
                           dependencies: externals_dep,
                           include_directories: siril_inc_dir,
                           c_args : siril_c_flag,
                           cpp_args : siril_cpp_flag)

siril_dep = declare_dependency(link_with: siril_lib,
                               dependencies: externals_dep,
                               include_directories: siril_inc_dir)


################################################################################
# Build gui executable

siril_gui = executable('siril', 'main.c', win_res,
                       dependencies: siril_dep,
                       link_whole: siril_lib,
                       link_args : siril_link_arg,
                       c_args : siril_c_flag,
                       cpp_args : siril_cpp_flag,
                       install : true,
                      )

################################################################################
# Build cli executable

siril_cli = executable('siril-cli', 'main-cli.c', win_res,
                       dependencies: siril_dep,
                       link_whole: siril_lib,
                       link_args : siril_link_arg,
                       c_args : siril_c_flag,
                       cpp_args : siril_cpp_flag,
                       install : true
                      )

################################################################################
# Install GUI files

install_data(['gui/siril.css', 'gui/siril3.glade', 'gui/siril-shortcuts.ui'],
             install_dir : siril_datadir)


################################################################################
# Tests

subdir('tests')
