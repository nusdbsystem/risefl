### the dependencies image
FROM ubuntu:20.04 as dependencies

LABEL maintainer="Yuncheng Wu <lemonwyc@gmail.com>"

# Install basic libraries
RUN apt-get update && apt-get upgrade -y && \
        apt-get install -y --no-install-recommends \
        build-essential \
        python \
        python3 \
        python3-dev \
        python3-pip \
        libgmp3-dev \
        pkg-config \
        wget \
        automake \
        swig \
        gdb \
        git \
        unzip \
        libtool \
        libssl-dev \
        libtbb-dev \
        libsodium-dev \
        libcurl4-openssl-dev \
        && \
        pip3 install requests && \
        apt-get clean && \
        rm -rf /var/lib/apt/lists/*

# Change default sh to bash
RUN dpkg-reconfigure dash

# Install libsodium 1.0.18
RUN mkdir /root/temp && \
    cd /root/temp && \
    wget https://download.libsodium.org/libsodium/releases/libsodium-1.0.18-stable.tar.gz && \
    tar -xzvf libsodium-1.0.18-stable.tar.gz && \
    cd libsodium-stable/ && \
    ./configure && \
    make && make check && \
    make install

# Install GMP 6.3.0
RUN cd /root/temp && \
    wget https://ftp.gnu.org/gnu/gmp/gmp-6.2.0.tar.xz && \
    tar -xf gmp-6.2.0.tar.xz && \
    cd gmp-6.2.0/ && \
    ./configure && \
    make && make check && \
    make install

# Install GMP 6.3.0
RUN cd /root/temp && \
    wget https://libntl.org/ntl-11.5.1.tar.gz && \
    tar -xzvf ntl-11.5.1.tar.gz && \
    cd ntl-11.5.1/src/ && \
    ./configure && \
    make && make check && \
    make install

# Upgrade cmake version to 3.29
RUN cd /root/temp && \
    git clone https://github.com/Kitware/CMake.git && \
    cd CMake && \
    ./bootstrap && \
    make -j$(nproc) && \
    make install && \
    hash -r && \
    cmake --version

# Download risefl and build
RUN cd /root && \
    git clone https://github.com/nusdbsystem/risefl.git && \
    cd risefl && \
    bash make.sh

# Copy the built library to python packages path \
# Here the python3.8 package path is /usr/local/lib/python3.8/dist-packages
# Can use python3 -m site to check the correct sys.path
RUN cd /root/risefl/build && \
    cp -t /usr/local/lib/python3.8/dist-packages librisefl_crypto.a _risefl_interface.so risefl_interface.py