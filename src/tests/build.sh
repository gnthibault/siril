#!/bin/sh
CC=cc
CXX=c++
LD=c++
CFLAGS="-Wall -I.. -O2 -fopenmp `pkg-config --cflags gtk+-3.0` `pkg-config --cflags cfitsio` `pkg-config --cflags gsl`"
LDFLAGS="-fopenmp `pkg-config --libs gtk+-3.0` `pkg-config --libs cfitsio` `pkg-config --libs gsl` -lm"

set -x
# compile the compare_fits tool
#$CC $CFLAGS -c -o compare_fits.o compare_fits.c &&
#$CC $CFLAGS -c -o dummy.o dummy.c &&
#$LD $LDFLAGS -o compare_fits compare_fits.o dummy.o ../io/image_format_fits.o ../core/utils.o ../gui/progress_and_log.o

# compile the sorting algorithm test
$CC $CFLAGS -c -o sorting.o sorting.c &&
$CC $CFLAGS -DUSE_ALL_SORTING_ALGOS -c -o ../algos/sorting.o ../algos/sorting.c &&
$CXX $CFLAGS -c -o ../rt/rt_algo.o ../rt/rt_algo.cc
$LD $LDFLAGS -o sorting sorting.o ../algos/sorting.o ../rt/rt_algo.o
