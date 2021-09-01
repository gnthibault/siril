crossroad source msys2
mkdir _deps && cd _deps

# Install deps from crossroad
crossroad install fftw \
                  exiv2 \
                  gtk3 \
                  libconfig \
                  gsl \
                  opencv \
                  libheif \
                  ffms2 \
                  cfitsio

# Build LibRaw from github
crossroad install lcms2
git clone --depth 1 https://github.com/LibRaw/LibRaw.git
cd LibRaw
autoreconf -fi
crossroad ./configure --disable-examples --disable-static
make install
cd ..