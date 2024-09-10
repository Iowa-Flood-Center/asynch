FROM registry.fedoraproject.org/fedora-minimal:34
#FROM fedora:latest

copy . /app

# Use this if you want fedora-minimal
RUN microdnf upgrade -y && microdnf install autoconf automake \
    openmpi-devel hdf5-devel \
    libpq-devel zlib glibc-devel \
    gcc-gfortran gcc make -y


# Use this if you want fedora
# RUN dnf upgrade -y && dnf install autoconf automake \
#     openmpi-devel hdf5-devel \
#     libpq-devel zlib glibc-devel \
#     gcc-gfortran gcc make -y


#Enable the BLOCK below if you are a WINDOWS user, choose the corresponding package manager for fedora-minimal or fedora
# RUN microdnf install dos2unix -y
# RUN dnf install dos2unix -y
# WORKDIR /app
# RUN dos2unix *
# RUN dos2unix docker_config_files/configure.sh
# RUN dos2unix examples/*


RUN ln -s /usr/lib64/openmpi/bin/* /usr/bin/

WORKDIR /app
RUN autoreconf --install
WORKDIR /app/docker_config_files
RUN chmod +x configure.sh
RUN ./configure.sh
WORKDIR /app/build
RUN make
RUN make install
WORKDIR /app/examples