################################################################################
# Meson file
################################################################################

if platform_linux
  subdir('linux')
elif platform_windows
  subdir('windows')
elif platform_osx
  subdir('os-x')
else
  error('Unsupported OS.')
endif