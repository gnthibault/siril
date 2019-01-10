#!/bin/sh
CC=clang
LD=clang
CFLAGS="-I.. -O2 -fopenmp `pkg-config --cflags gtk+-3.0` `pkg-config --cflags cfitsio` `pkg-config --cflags gsl`"
LDFLAGS="-fopenmp `pkg-config --libs gtk+-3.0` `pkg-config --libs cfitsio` `pkg-config --libs gsl` -lm"

set -x
$CC -Wall $CFLAGS -c -o compare_fits.o compare_fits.c &&
$CC -Wall $CFLAGS -UHAVE_STATS -c -o dummy.o dummy.c &&
$LD $LDFLAGS -o compare_fits compare_fits.o dummy.o ../io/image_format_fits.o ../core/utils.o ../gui/progress_and_log.o


$CC -Wall $CFLAGS -c -o sorting.o sorting.c &&
$CC -Wall -DHAVE_STATS $CFLAGS -c -o dummy.o dummy.c &&
$LD $LDFLAGS -o sorting sorting.o dummy.o ../core/utils.o ../gui/progress_and_log.o ../algos/statistics.o
