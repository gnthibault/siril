#!/bin/sh
CC=cc
CXX=c++
LD=c++
CFLAGS="-Wall -g -I.. -fopenmp `pkg-config --cflags gtk+-3.0` `pkg-config --cflags cfitsio` `pkg-config --cflags gsl`"
LDFLAGS="-Wl,--unresolved-symbols=ignore-all -fopenmp `pkg-config --libs gtk+-3.0` `pkg-config --libs cfitsio` `pkg-config --libs gsl` -lm"

set -x
# compile the compare_fits tool
# $CC $CFLAGS -c -o compare_fits.o compare_fits.c &&
# $CC $CFLAGS -c -o dummy.o dummy.c &&
# $LD $LDFLAGS -o compare_fits compare_fits.o dummy.o ../io/image_format_fits.o ../core/utils.o ../core/siril_log.o ../core/siril_date.o

# compile the sorting algorithm test (broken, use meson tests)
# $CC $CFLAGS -c -o sorting.o sorting.c &&
# $CC $CFLAGS -DUSE_ALL_SORTING_ALGOS -c -o ../algos/sorting.o ../algos/sorting.c &&
# $CXX $CFLAGS -c -o ../rt/rt_algo.o ../rt/rt_algo.cc
# $LD $LDFLAGS -o sorting sorting.o ../algos/sorting.o ../rt/rt_algo.o

# compile the stacking tests
$CC $CFLAGS -DWITH_MAIN -c -o stacking_blocks_test.o stacking_blocks_test.c &&
$CC $CFLAGS -DDUMMY_LOG -c -o dummy.o dummy.c &&
$LD $LDFLAGS -o stacking_blocks_test stacking_blocks_test.o ../stacking/median_and_mean.o ../core/utils.o dummy.o

# in mpp, a libsiril.a was created to not require the endless dummy list of functions:
# https://gitlab.com/free-astro/siril/-/blob/mpp/src/Makefile.am
