################################################################################
# Meson file for Siril scripts.
################################################################################

if platform_windows
  install_subdir('.', install_dir : 'scripts', exclude_files: 'meson.build', 
                  strip_directory : true)
else
  install_subdir('.', install_dir : siril_datadir + '/scripts', 
                exclude_files: 'meson.build', strip_directory : true)
endif