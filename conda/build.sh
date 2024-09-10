#!/bin/bash
autoreconf --install

mkdir build-conda && cd build-conda
CFLAGS="-O3 -DNDEBUG" ../configure --prefix=$PREFIX
make
#make check
make install
