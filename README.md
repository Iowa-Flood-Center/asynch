[![Build Status](https://travis-ci.org/Iowa-Flood-Center/asynch.svg?branch=master)](https://travis-ci.org/Iowa-Flood-Center/asynch)

# ASYNCH

A numerical library for solving differential equations with a tree structure. Emphasis is given to hillslope-link river basin models.


## Requirements

### Redhat & Centos

The same as **Fedora** but make sure you have the `epel-release` and the `PowerTools` repositories enabled.

```shell
 sudo dnf install epel-release
 sudo dnf config-manager --set-enabled PowerTools
 sudo dnf update
```

If you get an error when enabling `PowerTools` because of `config-manager`, run this:
```shell
 sudo dnf install 'dnf-command(config-manager)'
```

### Fedora

#### Packages

- autoconf
- automake
- openmpi-devel
- hdf5-devel
- libpq-devel
- zlib
- gcc-gfortran
- gcc

#### Optional
- openblas-devel

If you want ***openblas*** support with ***openmpi***, you can install the optional package, but it is not necessary for the app to work.

Run this command after installing the dependencies, this will create ***symbolic links*** of all the binaries of openmpi inside the `/usr/bin/` directory so the system can detect them:
```shell
sudo ln -s /usr/lib64/openmpi/bin/* /usr/bin/
```

### Ubuntu & Debian

#### Packages

- autoconf
- automake
- gcc
- make
- openmpi-bin
- libopenmpi-dev
- zlib1g-dev
- hdf5-tools
- libhdf5-dev
- libhdf5-openmpi-dev
- libhdf5-cpp-103
- libpq-dev
- pkg-config

#### Optional
- libopenblas64-0-openmp
- libopenblas64-0-openmp-dev
- libopenblas64-openmp-dev
- libopenblas-dev

If you want ***openblas*** support with ***openmpi***, you can install the optional package, but it is not necessary for the app to work.

## Compiling

Please run the following comands to compile `asynch`:

This will generate all the `configure` files and the `makefiles`.
```shell
autoreconf --install
cd build
../configure CFLAGS="-O3 -DNDEBUG -Wno-format-security"
make
make install
```

The default installation location in `/usr/local/`, you can find the `binaries` and  `libs` in `/usr/local/bin/` and `/usr/local/lib`.

If you want a custom location, you can use the `--prefix` option when running configure:

```shell
../configure CFLAGS="-O3 -DNDEBUG -Wno-format-security" --prefix=/custom/location/
```

## Docker

### Note/Disclaimer
The decision to use **fedora:latest** as the default container for docker was mainly because `asynch` needs `GLIBC > 2.29`, **Fedora 34** uses `GLIBC 2.33`, while **Redhat 8** and **Centos 8** uses `GLIBC 2.28`, I was having a little trouble setting up ``GLIBC`` inside the container for **centos:latest**. The app did compile, but is wasn't reading the ``GLIBC`` compilation I did for version ``2.29``. I think the app should be deployed using **Centos**, so if anyone wants to work around that issue, know that it is possible to solve and close to be solved.

Also, ***openblas*** was not working for docker, so *the docker compilation **does not** have **openblas support**.*

#### Windows

**Windows** users need to have `wsl2` installed and use a linux distro to run docker inside it.
**Enable the Windows instructions** inside the `Dockerfile` to be able to use the program correctly.


### Running Docker

Make sure you have **docker** installed in your machine, then follow this instructions:

1.  Clone this repository
2. Build the docker image, **this will take a long time**.
	```shell
	sudo docker build -t asynch-image .
	```
3. Log in into the image:
	```shell
	sudo docker run -it asynch-image
	```
4. By default, the dockerfile leaves you inside the `examples/` directory, so just see the **Running Example** instructions.

5. If you want to share data between your machine and the docker, run the followiing command:
	```shell
	sudo docker run -it -v /your/machine/examples:/docker-image/location/examples/ asynch-image
	```
	Change `/your/machine/examples/` for the location where you want your data to be located inside your machine, and change `/docker-image/location/examples/` to the location where you want your files to be inside the docker image.
	
    **Note:** All data is saved inside your machine.

6. Go to the location where your *example files* are inside the docker image, and execute them. **See Running Example.**

## Running Example

If everything when correctly, go to the `examples/` directory and execute the following command:

```shell
mpirun -n 4 asynch clearcreek.gbl
```

If you are running inside a docker container, add the `--allow-run-as-root` flag:

```shell
mpirun -n 4 --allow-run-as-root asynch clearcreek.gbl
```

The output should be the following:

```shell
Computations complete. Total time for calculations: 1.320425

Results written to file clearcreek.h5.
Peakflows written to file clearcreek.pea.
```

Inside the `clearcreek.pea` file, you should see the beginning exactly like this:

```shell
6359
254
```



## Documentation

The documentation is available [here](http://asynch.readthedocs.io/). Thank you to the people running Read the Docs for such an excellent service.

The source for the documentation is in the `docs` folder. Here are the instructions to built and read it locally. The documentation is built with [Doxygen](http://www.doxygen.org/) and [Sphinx](http://www.sphinx-doc.org). The sphinx template is from [ReadtheDocs](https://docs.readthedocs.io). [Breathe](https://breathe.readthedocs.io) provides a bridge between the Sphinx and Doxygen documentation systems.

    pip install --user sphinx sphinx-autobuild sphinx_rtd_theme breathe recommonmark
    apt-get doxygen

    cd docs  
    doxygen api.dox
    doxygen devel.dox
    make html

The html documentation is generated in `docs/.build/html`.

## Testing

Asynch doesn't have a good test covergage at the moment but the unit test framework is in place.
