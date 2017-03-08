# SIRIL
=================

Summary
-------
SIRIL is an astronomical image processing tool.

SIRIL is an image processing tool specially tailored for noise reduction and improving the
signal/noise ratio of an image from multiple captures, as required in astronomy.
SIRIL can align automatically or manually, stack and enhance pictures from various file formats,
even images sequences (movies and SER files).

Contributors are welcome. Programming language is C.
Main development is done with most recent versions of libraries.

Requirements
------------
 * Gtk 3, (>= 3.6) for the GUI toolkit
 * cfitsio for fits related stuff
 * fftw3 as a FFT library
 * The GNU Scientific Library (gsl) for FWHM implementation, histograms and background extraction
 * libconfig++ (>= 1.4) for structured configuration file management
 * g++ for opencv code and avi exporter

SIRIL works internally with FITS files, but other file formats can be used as
input and converted using the conversion tab of the control window. Some file
formats are handled internally, like BMP, PPM and SER, some required external
libraries or programs as listed below. Libraries need to be present at
compilation time, or their support won't be compiled.

 * libraw for DSLR RAW files import
 * libffms2 for films import
 * libtiff (>= 4)
 * libjpeg or compatible libraries like libjpeg-turbo
 * libpng
 * libopencv
 * libavformat
 * libavutil >= 55.20
 * libavcodec
 * libswscale
 * libswresample

All these libraries are available in most Linux distributions and free systems,
maybe with the exception of ffms2 that is not as popular as others and may need
to be compiled on some systems.

Since version [0.9.6](http://free-astro.org/index.php?title=Siril:0.9.6) a new optional dependency is required to plot photometry data.
The following package is not needed at compilation time:

 * gnuplot
 
Source download
---------------

You need to use the following commands to clone SIRIL.

Branch 0.9:

    svn co https://free-astro.org/svn/siril/branches/0.9 Siril-0.9
    
or you can use the mirror on github:

    git clone https://github.com/lock042/Siril-0.9.git 

Development version (highly unstable):

    svn co https://free-astro.org/svn/siril/trunk/ Siril-master
 

Building SIRIL for GNU/Linux and OS X
-------------------------------------
The install is similar to the usual GNU/Linux package build, except that the
configure script is not directly shipped and has to be created and run with the
following command:

    ./autogen.sh
    make
    sudo make install

Note that a binary package for stable version of SIRIL is maintained for Debian. 
PPA repositories for Ubuntu and Linux Mint and maintained by SIRIL's authors are
now available in ppa:lock042/siril.
See the download page for other packages that could be available.

Notes on SIRIL FITS image format
--------------------------------
Flexible Image Transport System (FITS) is an open standard defining a digital 
file format useful for storage, transmission and processing of scientific and other images.
FITS is the most commonly used digital file format in astronomy.

Since FITS doesn't specify the order and size of data, it's useful to fix it at
some point. Currently, SIRIL uses unsigned 16-bit per channel values (TUSHORT),
and images are stored channel after channel on a bottom-to-top, left-to-right
order.

All files imported and converted in SIRIL or files exported by SIRIL are in this
FITS format, except sequence files like SER and films, which are read from the
file and converted on-the-fly.

Notes on image sequence files
-----------------------------
SIRIL makes a strong case for the use SER sequences against the generic film
containers that are not well suited for astronomy data and that may not be read
the same way by different players. SIRIL can convert any film format supported
by FFMS2 (probably all ffmpeg formats, which is a lot) to SER, and even any
image sequence to SER.
SIRIL supports SER v3. See https://free-astro.org/index.php/SER for more details.

Useful links
------------
 * [Project Homepage](http://free-astro.org/index.php/Siril)
 * [Documentation](http://free-astro.org/siril_doc-en/#Reference_documentation_1)
 * [Releases and Downloads](http://free-astro.org/index.php?title=Siril:releases)
 * [Report a bug](https://free-astro.org/bugs/view_all_bug_page.php)
