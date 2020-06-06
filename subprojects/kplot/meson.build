################################################################################
# Meson main file
################################################################################

################################################################################
# Define project

project('kplot', ['c','cpp'], version : '0.1.14', license : ['ISC'], meson_version: '>=0.51.0', default_options : ['c_std=gnu11','cpp_std=c++11'])

################################################################################
# Meson modules and compilers

## Compilers
cc = meson.get_compiler('c')

################################################################################
# Dependencies

## Core dependencies
m_dep = cc.find_library('m', required : true)
cairo_dep = dependency('cairo', required : true)

################################################################################
# Configuration

## General and package configuration
conf_data = configuration_data()

## Test reallocarray
code = '''#include <stdlib.h>
int main(void) { return( ! reallocarray(NULL, 2, 2));}
'''

# TODO: Better crosscompilation support.
# TODO: Also has_exe_wrapper will be deperacted in meson version 0.55.
rcode = 1
if meson.has_exe_wrapper()
    result = cc.run(code, name : 'reallocarray check')
    rcode = result.returncode()
endif

if rcode == 0
  conf_data.set('HAVE_', true)
else
  conf_data.set('HAVE_', false)
endif

################################################################################
# Generate configuration files

configure_file(input : 'compat.h.in',
               output : 'compat.h',
               configuration : conf_data)

################################################################################
# List sources

inc = include_directories('.')

src_files = ['colours.c',
              'array.c',
              'border.c',
              'bucket.c',
              'buffer.c',
              'draw.c',
              'grid.c',
              'hist.c',
              'label.c',
              'kdata.c',
              'kplot.c',
              'margin.c',
              'mean.c',
              'plotctx.c',
              'reallocarray.c',
              'stddev.c',
              'tic.c',
              'vector.c']

################################################################################
# Build kplot library

kplot_lib = static_library('kplot', src_files, dependencies : [m_dep, cairo_dep], include_directories : inc)

################################################################################
# Dependency kplot

kplot_dep = declare_dependency(include_directories : inc, link_with : kplot_lib)