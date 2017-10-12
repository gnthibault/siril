dnl Process this file with autoconf to produce a configure script.
m4_define([siril_major_version], [0])
m4_define([siril_minor_version], [9])
m4_define([siril_micro_version], [8])
m4_define([siril_prerelease_version], [alpha])
m4_define([siril_version],
          [siril_major_version.siril_minor_version.siril_micro_version])
m4_define([siril_version_full],
          [siril_major_version().siril_minor_version().siril_micro_version()m4_bpatsubst(siril_prerelease_version(), [^], [-])])

dnl change siril_version_full by siril_version for stable release
AC_INIT([siril], [siril_version_full],
dnl AC_INIT([siril], [siril_version],
        [https://free-astro.org/bugs/view_all_bug_page.php],
        [siril],
        [https://free-astro.org/index.php/Siril/fr]) 
        
AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_HEADERS(config.h)
AC_CONFIG_SRCDIR([src/main.c])
AC_CONFIG_MACRO_DIR([m4])

dnl ================================================================
dnl Gettext stuff.
dnl ================================================================
IT_PROG_INTLTOOL([0.50])
AC_SUBST([GETTEXT_PACKAGE], [siril])

remember_set_CFLAGS="$CFLAGS"

AC_PROG_CC
AC_HEADER_STDC
AC_OPENMP

AC_LANG([C])
AC_LANG([C++])

# remove -g flag from default CFLAGS
if test "x$remember_set_CFLAGS" = "x"; then
	if test "$CFLAGS" = "-g -O2"; then
		CFLAGS="-O2"
	elif test "$CFLAGS" = "-g"; then
		CFLAGS=""
	fi
fi

CFLAGS="$CFLAGS $OPENMP_CFLAGS"
CXXFLAGS="$CFLAGS"

dnl Set PACKAGE_SOURCE_DIR in config.h.
packagesrcdir=`cd $srcdir && pwd`
AC_DEFINE_UNQUOTED(PACKAGE_SOURCE_DIR, "${packagesrcdir}", [Package source dir])

dnl Set PACKAGE_DOC_DIR in config.h.
dnl AC_DEFINE_UNQUOTED(PACKAGE_DOC_DIR, "${ac_default_prefix}/share/doc/${PACKAGE}", "package doc dir")

dnl ================ MANDATORY DEPENDECIES ================ 
dnl check math lib
AC_CHECK_LIB(m, log, LIBS="$LIBS -lm")

dnl check GTK
PKG_CHECK_MODULES(GTK, gtk+-3.0 >= 3.6.0)

dnl check fftw3
PKG_CHECK_MODULES(FFTW, [fftw3])

dnl check GNU Scientific Library
PKG_CHECK_MODULES(GSL, [gsl < 2],
		  [AC_DEFINE([HAVE_GSL_1], [1], [Use GSL 1])],
		  [PKG_CHECK_MODULES(GSL, [gsl >= 2],
				     [AC_DEFINE([HAVE_GSL_2], [1], [Use GSL 2])
				     ])
		  ])

dnl check cfitsio
PKG_CHECK_MODULES(CFITSIO, [cfitsio])

dnl check libconfig
PKG_CHECK_MODULES([LIBCONFIGXX], [libconfig++ >= 1.4],,
	AC_MSG_ERROR([libconfig++ 1.4 or newer not found.]))

dnl ================ OPTIONAL DEPENDECIES ================ 
dnl check raw lib
PKG_CHECK_MODULES(LIBRAW, [libraw],
 		  [AC_DEFINE([HAVE_LIBRAW], [1], [Using RAW images])],
 		  AC_MSG_WARN([libraw not found. Not using RAW importer.]))

dnl check tiff lib
PKG_CHECK_MODULES(LIBTIFF, [libtiff-4],
		  [AC_DEFINE([HAVE_LIBTIFF], [1], [Using TIFF images])],
		  AC_MSG_WARN([libtiff not found. Not using TIFF importer and exporter.]))

dnl check jpeg lib
dnl PKG_CHECK_MODULES(JPEG, [libjpeg])
AC_CHECK_LIB(jpeg, jpeg_mem_src, [],
	     AC_MSG_WARN([libjpeg not found. Not using JPEG importer and exporter.]))

dnl check png lib
PKG_CHECK_MODULES(LIBPNG, [libpng],
		  [AC_DEFINE([HAVE_LIBPNG], [1], [Using PNG images])],
		  AC_MSG_WARN([libpng not found. Not using PNG importer and exporter.]))

dnl check ffms2
PKG_CHECK_MODULES(FFMS2, [ffms2 < 2.21],
		  [AC_DEFINE([HAVE_FFMS2_1], [1], [Use FFMS2 < 2.21])],
		  [PKG_CHECK_MODULES(FFMS2, [ffms2 >= 2.21],
		  [AC_DEFINE([HAVE_FFMS2_2], [1], [Use FFMS2 > 2.21])],
		  AC_MSG_WARN([libffms2 not found. Not using AVI importer.]))
		  ])

dnl checking several ffmpeg libraries
PKG_CHECK_MODULES(FFMPEG, [libavformat libavutil >= 55.20 libavcodec libswscale libswresample],
		  [AC_DEFINE([HAVE_FFMPEG], [1], [film export supported])],
		  AC_MSG_WARN([ffmpeg libav libraries not found or too old. Not using film export formats.]))

dnl check opencv
AC_PROG_CXX
if test -n "$CXX" && ( test "X$CXX" != "Xno" &&
		      ( (test "X$CXX" = "Xg++" && `g++ -v >/dev/null 2>&1` ) ||
		      (test "X$CXX" != "Xg++"))) ; then
	AC_PROG_CXXCPP
	
	dnl check for C++11 support
	m4_include([m4/ax_cxx_compile_stdcxx.m4])
	AX_CXX_COMPILE_STDCXX(11, noext, mandatory)
	
	dnl AC_MSG_WARN(CXX is "$CXX")
	dnl http://nerdland.net/2009/07/detecting-c-libraries-with-autotools/ solution 2
	PKG_CHECK_MODULES([OPENCV], [opencv >= 3.0],
			  [AC_DEFINE([HAVE_OPENCV], [1], [Using OpenCV])]
			  [OPENCV_LIBS="-lopencv_core -lopencv_imgproc -lopencv_calib3d"],
			  [PKG_CHECK_MODULES([OPENCV], [opencv >= 2.0],
			  [AC_DEFINE([HAVE_OPENCV], [1], [Using OpenCV])]
			  [OPENCV_LIBS="-lopencv_core -lopencv_imgproc"],
			  AC_MSG_WARN([opencv not found. Not using some image processing.]))]
			  )

else
	AC_MSG_ERROR([no C++ compiler found])
fi

AM_CONDITIONAL(HAVE_OPENCV, [test "x$OPENCV_LIBS" != "x"])

dnl Compute the canonical host-system type variable (need config.guess and config.sub)
AC_CANONICAL_HOST

#################
# Check for Win32
#################

AC_MSG_CHECKING([if compiling for Win32])
case "$host_os" in
  mingw* | cygwin*)
    platform_win32=yes
    ;;
  *)
    platform_win32=no
    ;;
esac
AC_MSG_RESULT([$platform_win32])
AM_CONDITIONAL(PLATFORM_WIN32, test "$platform_win32" = "yes")

PATHSEP=':'
os_win32=no
if test "x$platform_win32" != xno; then
  AC_MSG_CHECKING([if this is native Win32])
  case "$host_os" in
    mingw*)
      os_win32=yes
      case "$host_cpu" in
        x86_64)
        ;;
        *)
        WIN32_LARGE_ADDRESS_AWARE='-Wl,--large-address-aware'
        ;;
      esac
      PATHSEP=';'
      ;;
    *)
      ;;
  esac
  AC_MSG_RESULT([$os_win32])
fi
AC_SUBST(WIN32_LARGE_ADDRESS_AWARE)
AC_SUBST(PATHSEP)
AM_CONDITIONAL(OS_WIN32, test "$os_win32" = "yes")
AM_CONDITIONAL(OS_UNIX, test "$os_win32" != "yes")

if test "$os_win32" = "yes"; then
  AC_CHECK_PROG(ms_librarian, lib.exe, yes, no)
  AC_CHECK_TOOL(WINDRES, windres, :)
else
  WINDRES=":"
fi

AM_CONDITIONAL(MS_LIB_AVAILABLE, test "x$ms_librarian" = xyes)
AM_CONDITIONAL(HAVE_WINDRES, test "x$WINDRES" != "x:")

AC_SUBST(WINDRES)

dnl check for osx platform
AC_MSG_CHECKING([for Mac OS X platform])
case "$host" in
  *-*-darwin*)
    platform_osx=yes
    ;;
  *)
    platform_osx=no
    ;;
esac
AC_MSG_RESULT([$platform_osx])
AM_CONDITIONAL(PLATFORM_OSX, test "$platform_osx" = "yes")

dnl check for native osx
gdk_targets=`$PKG_CONFIG --variable=targets gdk-3.0`

gdk_windowing="x11"

for target in $gdk_targets;
do
	case "$target" in
		quartz)
			gdk_windowing=quartz
			;;
		win32)
			gdk_windowing=win32
			;;
	esac
done

AC_MSG_CHECKING([for native Mac OS X])
if test "$gdk_windowing" = "quartz"; then
	os_osx=yes
else
	os_osx=no
fi
AC_MSG_RESULT([$os_osx])
AM_CONDITIONAL(OS_OSX, test "$os_osx" = "yes")

AM_CONDITIONAL(GDK_WINDOWING_X11, test "$gdk_windowing" = "x11")

if test "$os_osx" = "yes"; then
	AC_DEFINE([OS_OSX],[1],[Defined if os is Mac OSX])

	PKG_CHECK_MODULES(GTK_MAC, gtk-mac-integration-gtk3 >= 2.0.8)

	AC_SUBST(GTK_MAC_CFLAGS)
	AC_SUBST(GTK_MAC_LIBS)
fi

if test `uname` != "Darwin"; then
  if (test "$platform_win32" = "yes" || test "os_win32" = "yes" ) ; then
    LDFLAGS="$LDFLAGS -Wl,--export-all-symbols -mwindows"
  else
	LDFLAGS="$LDFLAGS -rdynamic"
  fi
fi

AC_CONFIG_FILES([Makefile src/Makefile po/Makefile.in])
AC_OUTPUT

