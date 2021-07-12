################################################################################
# Meson main file
################################################################################

################################################################################
# Define project

project('wcslib', 
  ['c','cpp'],
  version : '7.7.0',
  meson_version: '>=0.53.0',
  default_options : [
    'c_std=gnu11',
    'cpp_std=c++11'
    ],
)

################################################################################
# Meson modules and compilers

## Compilers
cc = meson.get_compiler('c')

################################################################################
# Dependencies

## Core dependencies
cfitsio_dep = dependency('cfitsio', required : true)


################################################################################
# Configuration

## General and package configuration
conf_data = configuration_data()

################################################################################
# Generate configuration files

################################################################################
# List sources

inc = include_directories('.')

src_files = ['cel.c',
            'dis.c',
            'fitshdr.c',
            'getwcstab.c',
            'lin.c',
            'log.c',
            'prj.c',
            'spc.c',
            'sph.c',
            'spx.c',
            'tab.c',
            'wcsbth.c',
            'wcs.c',
            'wcserr.c',
            'wcsfix.c',
            'wcshdr.c',
            'wcspih.c',
            'wcsprintf.c',
            'wcstrig.c',
            'wcsulex.c',
            'wcsunits.c',
            'wcsutil.c',
            'wcsutrn.c'
]

################################################################################
# Build kplot library

wcslib_lib = static_library('wcslib', src_files, dependencies : [cfitsio_dep], include_directories: '.')

################################################################################
# Dependency kplot

wcslib_dep = declare_dependency(include_directories: '.', link_with : wcslib_lib)
