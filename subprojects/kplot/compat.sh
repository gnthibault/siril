#! /bin/sh

$CC $CFLAGS -o test-$1 test-$1.c 1>/dev/null 2>&1 || exit

./test-$1 1>/dev/null 2>&1 || exit

echo "#define " `echo HAVE_$2 | tr '[:lower:]' '[:upper:]'`
