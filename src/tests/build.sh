#!/bin/sh
CC=clang
LD=clang
CFLAGS="-I.. -U_OPENMP `pkg-config --cflags gtk+-3.0` `pkg-config --cflags cfitsio`"
LDFLAGS="`pkg-config --libs gtk+-3.0` `pkg-config --libs cfitsio`"

set -x
$CC -Wall $CFLAGS -c -o compare_fits.o compare_fits.c &&
$CC -Wall $CFLAGS -c -o dummy.o dummy.c &&
$LD $LDFLAGS -o compare_fits compare_fits.o dummy.o ../io/image_format_fits.o ../core/utils.o ../gui/progress_and_log.o


$CC -Wall $CFLAGS -c -o sorting.o sorting.c &&
$LD $LDFLAGS -o sorting sorting.o dummy.o ../core/utils.o ../gui/progress_and_log.o
