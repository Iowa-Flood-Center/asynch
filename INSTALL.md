## Updating the package

```
autoreconf --install
```


## Installing the package

For an out of source build (prefered method):

```
mkdir build && cd build
../configure
make
make check
make install
```

### Installing the package on NEON

First, `git clone` the repository or `tar xzf` a released packages. To install the software for the IFC group:

```
mkdir build && cd build
../configure --prefix=/Groups/IFC/.local
make
make check
make install
```

## Standard Makefile Targets


`make all` Build programs, libraries, documentation, etc.
(Same as `make`.)
`make install` Install what needs to be installed.
`make install-strip` Same as `make install`, then strip debugging
symbols.
`make uninstall` The opposite of `make install`.
`make clean` Erase what has been built (the opposite of `make
all`).
`make distclean` Additionally erase anything `./configure`
created.
`make check` Run the test suite, if any.
`make installcheck` Check the installed programs or libraries, if
supported.
`make dist` Create PACKAGE-VERSION.tar.gz.

