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

project('siril',
  ['c','cpp'],
  version : '0.99.7',
  license : ['GPL3'],
  meson_version: '>=0.53.0',
  default_options : [
    'c_std=gnu11',
    'cpp_std=c++11',
    'buildtype=debugoptimized',
  ],
)

enable_openmp = get_option('openmp')

# General version
siril_version = meson.project_version()

version_array = siril_version.split('.')
siril_version_major = version_array[0].to_int()
siril_version_minor = version_array[1].to_int()
siril_version_micro = version_array[2].to_int()

stable = (siril_version_minor % 2 == 0)

## Get version from git
git = find_program('git', required: false)

if not git.found()
  message('Git not installed or git directory not found!')
  git_version = 'Unknown'
  git_last_commit_year = 'Unknown'
else
  r = run_command('git', 'show', '-s', '--format="%ci %h"','-n1')
  if r.returncode() == 0
    git_info = r.stdout().strip().split('"')[1]
    git_info_array = git_info.split(' ')
    git_version = git_info_array[3]
    git_last_commit_year = git_info_array[0].split('-')[0]
    message('Siril version:' + git_version)
    errortxt = r.stderr().strip()
  else
    message('Cannot retrive git version')
    message(r.stderr().strip())
    git_version = 'Unknown'
    git_last_commit_year = 'Unknown'
  endif
endif

################################################################################
# Host system detection

host_cpu_family = host_machine.cpu_family()
host_os = host_machine.system().to_lower()
message('Host os: ' + host_os)

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

platform_bsd = (
  host_os.contains('bsd')
)

################################################################################
# Meson modules and compilers

## Modules
cmake = import('cmake')
i18n = import('i18n')
if platform_windows
  winm = import('windows')
endif
gnome = import('gnome')

## Compilers
cc = meson.get_compiler('c')
cpp = meson.get_compiler('cpp')

################################################################################
# Set compiler and linker flags

siril_cpp_flag = []
siril_c_flag = []
siril_link_arg = []

# -fno-math-errno is used to speedup code, this for example speeds up all sqrt() functions.
# There is no case where we check the math errno. So it’s safe to use this flag for siril.
siril_c_flag += ['-DHAVE_CONFIG_H', '-fcommon', '-fno-math-errno', '-D_FILE_OFFSET_BITS=64']
siril_cpp_flag += ['-DHAVE_CONFIG_H', '-fcommon', '-fno-math-errno', '-D_FILE_OFFSET_BITS=64']

if platform_linux
  siril_link_arg += '-rdynamic'
elif platform_bsd
  siril_link_arg += '-Wl,--no-as-needed' # workaround for https://gitlab.com/free-astro/siril/-/issues/614 (Bug#251206)
  siril_link_arg += '-rdynamic'
  siril_link_arg += '-lexecinfo'
elif platform_windows
  #TODO: Do Meson provide more elegant way?
  siril_link_arg += '-Wl,--export-all-symbols'
  siril_link_arg += '-lpsapi'
  siril_link_arg += '-ldbghelp'
  if stable
    siril_link_arg += '-mwindows'
  endif
  siril_c_flag += ['-D_WIN32_WINNT=0x0600', '-DWINVER=0x0600']
  siril_cpp_flag += ['-D_WIN32_WINNT=0x0600', '-DWINVER=0x0600']
elif platform_osx
  # on OSX ObjC and C sources are mixed so adding objc to the linkflags
  siril_link_arg += ['-Wl,-framework,Foundation', '-Wl,-framework,AppKit', '-ObjC']
  add_project_arguments('-ObjC', language: [ 'c'])
  add_project_link_arguments(siril_link_arg, language: [ 'c', 'cpp' ])
else
  error('Unsupported OS.')
endif

# sse2 assembly
have_sse2 = false
if cc.has_argument('-msse2') and cpp.has_argument('-msse2')
  message('sse2 assembly available')
  if cc.compiles('asm ("punpckhwd %xmm0,%xmm1");')
    have_sse2 = true
  endif
endif

# sse4.1 assembly
have_sse4_1 = false
if cc.has_argument('-msse4.1') and cpp.has_argument('-msse4.1')
  message('sse4.1 assembly available')
  if cc.compiles('asm ("pmovzxbd %xmm0,%xmm1");')
    have_sse4_1 = true
  endif
endif

################################################################################
# Directories 

siril_share = get_option('prefix') + '/' + get_option('datadir') + '/'
siril_datadir = siril_share + meson.project_name() + '/'
siril_locale = get_option('prefix') + '/' + get_option('localedir') + '/'

################################################################################
# Dependencies

if get_option('relocatable-bundle') == 'yes'
  relocatable_bundle = true
elif get_option('relocatable-bundle') == 'no'
  relocatable_bundle = false
else # == 'platform-default'
  # By default, assume building for Windows or macOS everything to be on
  # the same prefix and can be relocated.
  # On other platforms, build-time paths are meaningful.
  if platform_windows or platform_osx
      relocatable_bundle = true
  else
      relocatable_bundle = false
  endif
endif

if get_option('enable-libcurl') == 'yes'
  enable_libcurl= true
elif get_option('enable-libcurl') == 'no'
  enable_libcurl = false
else # == 'platform-default'
  if platform_windows or platform_osx
      enable_libcurl = true
  else
      enable_libcurl = false
  endif
endif

# Dependencies
no_dep = dependency('', required: false)

## Core dependencies
dbghelp = platform_windows ? cc.find_library('dbghelp') : no_dep
glib_dep = dependency('glib-2.0', required : true, version: '>= 2.56.0')
gio_dep = dependency('gio-2.0', required : true, version: '>= 2.56.0')
if platform_windows
  gio_dep = dependency('gio-windows-2.0', required : true, version: '>= 2.56.0')
else
  gio_dep = dependency('gio-unix-2.0', required : true, version: '>= 2.56.0')
endif
gtk_dep = dependency('gtk+-3.0', required : true, version: '>= 3.20.0')
gsl_dep = dependency('gsl', required : true)
libconfig_dep = dependency('libconfig', required : true)
m_dep = cc.find_library('m', required : true)
exiv_dep = dependency('exiv2', required : true, version: '>= 0.25')
fftw3f_dep = dependency('fftw3f', required : true)
cfitsio_dep = dependency('cfitsio', required : true)
cairo_dep = dependency('cairo', required : true)
kplot_dep = dependency('kplot', fallback : ['kplot', 'kplot_dep'], required : true)


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
  sub_proj = cmake.subproject('librtprocess', cmake_options: '-DBUILD_SHARED_LIBS=OFF')
  rtprocess_dep = sub_proj.dependency('rtprocess')
  rtprocess_inc = sub_proj.include_directories('rtprocess')
else
  siril_cpp_flag += '-DHAVE_LIBRTPROCESS'
endif

## Optionnal dependencies

gio_unix_dep = dependency('gio-unix-2.0', required : false, version: '>= 2.56.0')

openmp_dep = []
if enable_openmp
  if platform_osx
    message('Set OpenMP flags for macOS')
    siril_link_arg += ['-lomp']
    siril_cpp_flag += ['-Xpreprocessor', '-fopenmp']
    siril_c_flag += ['-Xpreprocessor', '-fopenmp']
  else
  openmp_dep = dependency('openmp', required : false)
  endif
else
  message('OpenMP disabled')
endif

libraw_dep = dependency('libraw', required : false)
libtiff_dep = dependency('libtiff-4', required : false)
libjpeg_dep = dependency('libjpeg', required : false)
libpng_dep = dependency('libpng', required : false, version: '>= 1.6.0')
libheif_dep = dependency('libheif', required : false)
ffms2_dep = dependency('ffms2', required : false)
ffmpeg = [
  dependency('libavformat', required : false),
  dependency('libavcodec', required : false),
  dependency('libavutil', version : '>= 55.20', required : false),
  dependency('libswscale', required : false),
  dependency('libswresample', required : false)
]
libcurl_dep = no_dep
if enable_libcurl
  libcurl_dep = dependency('libcurl', required : false)
endif
criterion_dep = dependency('criterion', required : false)

wcslib_dep = dependency('wcslib', fallback : ['wcslib', 'wcslib_dep'], required : false)

################################################################################
# Configuration

## General and package configuration
conf_data = configuration_data()

conf_data.set('SIRIL_UNSTABLE', stable ? false : 1)
conf_data.set('SIRIL_OUTPUT_DEBUG', get_option('buildtype').startswith('debug') ? 1 : false)
conf_data.set('ENABLE_RELOCATABLE_RESOURCES', relocatable_bundle ? 1 : false)
conf_data.set('OS_OSX', platform_osx ? 1 : false)

conf_data.set_quoted('PACKAGE', meson.project_name(), description : 'Name of package.')
conf_data.set_quoted('PACKAGE_NAME', meson.project_name(), description : 'Define to the full name of this package.')
conf_data.set_quoted('PACKAGE_VERSION', siril_version, description : 'Define to the home page for this package.')
conf_data.set_quoted('VERSION', siril_version, description : 'Define to the home page for this package.')
conf_data.set_quoted('PACKAGE_SOURCE_DIR', 'src', description : 'Package source dir.')
conf_data.set_quoted('PACKAGE_STRING', meson.project_name() + ' ' + siril_version, description : 'Define to the full name and version of this package.')
conf_data.set_quoted('PACKAGE_TARNAME', meson.project_name(), description : 'Define to the one symbol short name of this package.')
conf_data.set_quoted('PACKAGE_URL', 'https://www.siril.org', description : 'Define to the home page for this package.')
conf_data.set_quoted('PACKAGE_BUGREPORT', 'https://gitlab.com/free-astro/siril/issues', description : 'Define to the address where bug reports for this package should be sent.')
conf_data.set_quoted('PACKAGE_DATA_DIR', siril_datadir, description : 'Define to the data directory.')
conf_data.set_quoted('LOCALEDIR', siril_locale, description : 'Define to the locale directory.')

# Check for available headers
foreach header : [
    { 'm': 'HAVE_EXECINFO_H',     'v': 'execinfo.h' },
    { 'm': 'HAVE_SYS_STATVFS_H',  'v': 'sys/statvfs.h' },
    { 'm': 'HAVE_SYS_VFS_H',      'v': 'sys/vfs.h' },
    { 'm': 'HAVE_SYS_PARAM_H',    'v': 'sys/param.h' },
  ]
  conf_data.set(header['m'], cc.has_header(header['v']) ? 1 : false)
endforeach

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
  conf_data.set('HAVE_LIBJPEG', true, description : 'Using JPEG images.')
else 
  conf_data.set('HAVE_LIBJPEG', false, description : 'Using JPEG images.')
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

ffmpeg_found = false
foreach dep : ffmpeg
  if not dep.found()
    ffmpeg_found = false
    warning('An FFmpeg library (' + dep.name() + ') not found or too old. Not using film export formats.')
    break
  else
    ffmpeg_found = true
  endif
endforeach
conf_data.set('HAVE_FFMPEG', ffmpeg_found, description : 'Using FFMS.')

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
endif

if wcslib_dep.found()
  conf_data.set('HAVE_WCSLIB', true, description : 'Using wcslib.')
else 
  conf_data.set('HAVE_WCSLIB', false, description : 'Using wcslib.')
  warning('wcslib not found. Not using WCS API.')
endif

conf_data.set('HAVE_CXX11', true, description : 'Define to 1 if cpp compiler is C++11 compliant.')

## Git configuration
git_data = configuration_data()
git_data.set_quoted('SIRIL_GIT_VERSION_ABBREV', git_version)
git_data.set_quoted('SIRIL_GIT_LAST_COMMIT_YEAR', git_last_commit_year)

################################################################################
# Externals dependencies

externals_dep = declare_dependency(
  dependencies: [
    dbghelp,
    m_dep,
    glib_dep,
    gtk_dep,
    libcurl_dep,
    gsl_dep,
    opencv_dep,
    rtprocess_dep,
    exiv_dep,
    wcslib_dep,
    libconfig_dep,
    fftw3f_dep,
    cfitsio_dep,
    libraw_dep,
    libtiff_dep,
    libjpeg_dep,
    libpng_dep,
    libheif_dep,
    ffms2_dep,
    openmp_dep,
    kplot_dep,
    ffmpeg
  ]
)

################################################################################
# Sources directories

subdir('platform-specific')
subdir('build')
subdir('src')
subdir('pixmaps')
subdir('data')
subdir('po')
subdir('scripts')
subdir('3rdparty')

################################################################################
# Install misc

man_files = [
  configure_file(
    input : meson.project_name() + '.1',
    output: meson.project_name() + '.1',
    configuration: conf_data,
  ),
  configure_file(
    input : meson.project_name() + '.1',
    output: meson.project_name() + '-cli.1',
    configuration: conf_data,
  ),
]

install_man(man_files)
install_data(['AUTHORS', 'LICENSE.md'],
             rename : ['doc/siril/AUTHORS', 'doc/siril/LICENSE.md'],
             install_dir : siril_share)


################################################################################
# Build summary

summary(
  {
    'Building Siril with prefix': get_option('prefix'),
  }, section: 'Directories'
)

summary(
  {
    'sse2'      : have_sse2,
    'sse4.1'    : have_sse4_1,
  }, section: 'Processor extensions'
)

summary(
  {
    'openmp'    : get_option('openmp'),
  }, section: 'Optional features'
)

summary(
  {
    'libraw'    : libraw_dep.found(),
    'TIFF'      : libtiff_dep.found(),
    'JPEG'      : libjpeg_dep.found(),
    'PNG'       : libpng_dep.found(),
    'HEIF'      : libheif_dep.found(),
    'FFMS2'     : ffms2_dep.found(),
    'FFMPEG'    : ffmpeg_found,
    'libcurl'   : libcurl_dep.found(),
    'wcslib'    : wcslib_dep.found(),
    'criterion' : criterion_dep.found(),
  }, section: 'Optional dependencies'
)
