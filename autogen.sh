#!/bin/sh
# Run this to generate all the initial makefiles, etc.

srcdir=`dirname $0`
test -z "$srcdir" && srcdir=.
ORIGDIR=`pwd`
cd $srcdir

PKG_NAME=`autoconf --trace 'AC_INIT:$1' "$srcdir/configure.ac"`

if [ "$#" = 0 -a "x$NOCONFIGURE" = "x" ]; then
        echo "**Warning**: I am going to run \`configure' with no arguments." >&2
        echo "If you wish to pass any to it, please specify them on the" >&2
        echo \`$0\'" command line." >&2
        echo "" >&2
fi

# Manage librtprocess automatically, for the first autogen.
# To update librtprocess, manual intervention, like 'git submodule update' then cmake and make, is required.
cd deps
cd librtprocess || ( echo 'Failed to get librtprocess, please download it and compile it yourself. Aborting.' && exit 1 )
if [ ! -d build ] ; then
	CMAKE_FLAGS="$(uname -s | grep '^MSYS' && echo ' -G "MSYS Makefiles ')"
	mkdir build && cd build &&
	cmake $CMAKE_FLAGS -DCMAKE_BUILD_TYPE="Release" -DWITH_STATIC_LIB=ON -DOPTION_OMP=OFF .. || exit 1
	cd ..
fi
cd ../..

set -x
aclocal --install || exit 1
intltoolize --force --copy --automake || exit 1
autoreconf --verbose --force --install -Wno-portability || exit 1
{ set +x; } 2>/dev/null

cd $ORIGDIR

if [ "$NOCONFIGURE" = "" ]; then
        set -x
        $srcdir/configure "$@" || exit 1
        { set +x; } 2>/dev/null

        if [ "$1" = "--help" ]; then exit 0 else
                echo "Now type \`make\' to compile $PKG_NAME" || exit 1
        fi
else
        echo "Skipping configure process."
fi
