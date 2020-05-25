# This file is part of Siril, an astronomy image processor.
# Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
# Copyright (C) 2012-2020 team free-astro (see more in AUTHORS file)
# Reference site is https://free-astro.org/index.php/Siril
#
# Siril is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Siril is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Siril. If not, see <http://www.gnu.org/licenses/>.


################################################################################
# Meson main file
################################################################################


################################################################################
# Define project

project('siril', ['c','cpp'], version : '0.99.1', license : ['GPL3'], meson_version: '>=0.54.0', default_options : ['c_std=gnu11','cpp_std=c++11'])

## Developement version
version_array = meson.project_version().split('.')
if (version_array[1].to_int() % 2) == 1
  message('Building development version')
  siril_unstable = true
else
  message('Building release version')
  siril_unstable = false
endif

## Get version from git
r = run_command('git', 'show', '-s', '--format="%ci %h"','-n1')

if r.returncode() != 0
  message('Git not installed or git directory not found!')
  message(r.stderr().strip())
  git_version = 'Unknown'
  git_last_commit_year = 'Unknown'
else
  git_info = r.stdout().strip().split('"')[1]
  git_info_array = git_info.split(' ')
  git_version = git_info_array[3]
  git_last_commit_year = git_info_array[0].split('-')[0]
  message('Siril version:', git_version)
  errortxt = r.stderr().strip()
endif

################################################################################
# Host system detection

host_cpu_family = host_machine.cpu_family()
host_os = host_machine.system().to_lower()

platform_linux = (
  host_os.contains('linux')
)

platform_windows = (
  host_os.contains('mingw') or
  host_os.contains('cygwin') or
  host_os.contains('windows')
)

platform_osx = (
  host_os.contains('machten') or
  host_os.contains('rhapsody') or
  host_os.contains('darwin')
)

################################################################################
# Meson modules and compilers

## Modules
cmake = import('cmake')
fs = import('fs')
i18n = import('i18n')

## Compilers
cc = meson.get_compiler('c')
cpp = meson.get_compiler('cpp')

siril_cpp_flag = []
siril_c_flag = []

################################################################################
# Directories 

siril_share = get_option('prefix') + '/' + get_option('datadir') + '/'
siril_datadir = siril_share + meson.project_name() + '/'
siril_locale = get_option('prefix') + '/' + get_option('localedir') + '/'

################################################################################
# Dependencies

## Core dependencies
glib_dep = dependency('glib-2.0', required : true)
gtk_dep = dependency('gtk+-3.0', required : true)
gsl_dep = dependency('gsl', required : true)
libconfig_dep = dependency('libconfig', required : true)
m_dep = cc.find_library('m', required : true)
exiv_dep = dependency('exiv2', required : true)
fftw3f_dep = dependency('fftw3f', required : true)
cfitsio_dep = dependency('cfitsio', required : true)
cairo_dep = dependency('cairo', required : true)
kplot_dep = dependency('kplot', fallback : ['kplot', 'kplot_dep'])

# TODO: Prefered OpenCV version???
opencv4_dep = dependency('opencv4', required : false)
opencv3_dep = dependency('opencv', required : false)
if not opencv4_dep.found()
  if opencv3_dep.found()
    message('Use OpenCV3')
    opencv_dep = opencv3_dep
  else
    error('OpenCV not found!')
  endif
else
  message('Use OpenCV4')
  opencv_dep = opencv4_dep
endif

rtprocess_dep = dependency('rtprocess', required : false)
if not rtprocess_dep.found()
  # TODO: Static link but still install the .a lib???
  sub_proj = cmake.subproject('librtprocess', cmake_options: '-DBUILD_SHARED_LIBS=OFF')
  rtprocess_dep = sub_proj.dependency('rtprocess')
  rtprocess_inc = sub_proj.include_directories('rtprocess')
else
  siril_cpp_flag += '-DHAVE_LIBRTPROCESS'
endif

## Optionnal dependencies
openmp_dep = dependency('openmp', required : false)
libraw_dep = dependency('libraw', required : false)
libtiff_dep = dependency('libtiff-4', required : false)
libjpeg_dep = dependency('libjpeg', required : false)
libpng_dep = dependency('libpng', required : false)
libheif_dep = dependency('libheif', required : false)
ffms2_dep = dependency('ffms2', required : false)
libcurl_dep = dependency('libcurl', required : false)

################################################################################
# Configuration

## General and package configuration
conf_data = configuration_data()

conf_data.set_quoted('PACKAGE', meson.project_name(), description : 'Name of package.')
conf_data.set_quoted('PACKAGE_NAME', meson.project_name(), description : 'Define to the full name of this package.')
conf_data.set_quoted('PACKAGE_VERSION', meson.project_version(), description : 'Define to the home page for this package.')
conf_data.set_quoted('VERSION', meson.project_version(), description : 'Define to the home page for this package.')
conf_data.set_quoted('PACKAGE_SOURCE_DIR', 'src', description : 'Package source dir.')
conf_data.set_quoted('PACKAGE_STRING', meson.project_name() + ' ' + meson.project_version(), description : 'Define to the full name and version of this package.')
conf_data.set_quoted('PACKAGE_TARNAME', meson.project_name(), description : 'Define to the one symbol short name of this package.')
conf_data.set_quoted('PACKAGE_URL', 'https://www.siril.org', description : 'Define to the home page for this package.')
conf_data.set_quoted('PACKAGE_BUGREPORT', 'https://gitlab.com/free-astro/siril/issues', description : 'Define to the address where bug reports for this package should be sent.')
conf_data.set_quoted('PACKAGE_DATA_DIR', siril_datadir, description : 'Define to the data directory.')
conf_data.set_quoted('LOCALEDIR', siril_locale, description : 'Define to the locale directory.')

## System configuration
if target_machine.system() == 'linux'
  conf_data.set('HAVE_BACKTRACE', 1, description : 'Define to 1 if you have the backtrace function.')
  conf_data.set('HAVE_EXECINFO_H', 1, description : 'Define to 1 if you have the backtrace function.')
else 
  conf_data.set('HAVE_BACKTRACE', 0, description : 'Define to 1 if you have the backtrace function.')
  conf_data.set('HAVE_EXECINFO_H', 0, description : 'Define to 1 if you have the backtrace function.')
endif

## Dependencies configuration
if gsl_dep.found()
  if gsl_dep.version().version_compare('>=2.0')
    conf_data.set('HAVE_GSL_1', 0, description : 'Use GSL 1.')
    conf_data.set('HAVE_GSL_2', 1, description : 'Use GSL 2.')
  else
    conf_data.set('HAVE_GSL_1', 1, description : 'Use GSL 1.')
    conf_data.set('HAVE_GSL_2', 0, description : 'Use GSL 2.')
  endif
  else
    conf_data.set('HAVE_GSL_1', 0, description : 'Use GSL 1.')
    conf_data.set('HAVE_GSL_2', 0, description : 'Use GSL 2.')
endif

if libraw_dep.found()
  conf_data.set('HAVE_LIBRAW', true, description : 'Using RAW images.')
else 
  conf_data.set('HAVE_LIBRAW', false, description : 'Using RAW images.')
  warning('libraw not found. Not using RAW importer.')
endif

if libtiff_dep.found()
  conf_data.set('HAVE_LIBTIFF', true, description : 'Using TIFF images.')
else 
  conf_data.set('HAVE_LIBTIFF', false, description : 'Using TIFF images.')
  warning('libtiff-4 not found. Not using TIFF importer.')
endif

if libjpeg_dep.found()
  conf_data.set('HAVE_JPEG', true, description : 'Using JPEG images.')
else 
  conf_data.set('HAVE_JPEG', false, description : 'Using JPEG images.')
  warning('libjpeg not found. Not using JPEG importer.')
endif

if libpng_dep.found()
  conf_data.set('HAVE_LIBPNG', true, description : 'Using PNG images.')
else 
  conf_data.set('HAVE_LIBPNG', false, description : 'Using PNG images.')
  warning('libpng not found. Not using PNG importer.')
endif

if libheif_dep.found()
  conf_data.set('HAVE_LIBHEIF', true, description : 'Using HEIF.')
else 
  conf_data.set('HAVE_LIBHEIF', false, description : 'Using HEIF.')
  warning('libheif not found. Not using HEIF importer.')
endif

if ffms2_dep.found()
  conf_data.set('HAVE_FFMS2', true, description : 'Using FFMS.')
else 
  conf_data.set('HAVE_FFMS2', false, description : 'Using FFMS.')
  warning('ffms2 not found. Not using FFMS.')
endif

if libcurl_dep.found()
  conf_data.set('HAVE_LIBCURL', true, description : 'Using Curl.')
else 
  conf_data.set('HAVE_LIBCURL', false, description : 'Using Curl.')
  warning('libcurl not found. Not using Curl.')
endif

conf_data.set('HAVE_CXX11', true, description : 'Define to 1 if cpp compiler is C++11 compliant.')

## Git configuration
git_data = configuration_data()
git_data.set_quoted('SIRIL_GIT_VERSION_ABBREV', git_version)
git_data.set_quoted('SIRIL_GIT_LAST_COMMIT_YEAR', git_last_commit_year)

################################################################################
# Sources directories

subdir('po')
subdir('src')
subdir('pixmaps')
subdir('platform-specific')

################################################################################
# Install

install_man(meson.project_name() + '.1')
install_data(['AUTHORS', 'LICENSE.md'],
             rename : ['doc/siril/AUTHORS', 'doc/siril/LICENSE.md'],
             install_dir : siril_share)

################################################################################
# Post build

warning('''
        MESON BUILD IS STILL EXPERIMENTAL!

        we encourage usage of meson if you want to test and
        report bugs. Yet if you are a packager, DO NOT USE
        meson until we fix all known issues and bugs.
        ''')