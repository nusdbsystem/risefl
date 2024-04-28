# API for RiseFL, the Robust Secure FL system, built on C/C++ core

## Installation:

RiseFL mainly require the following dependencies, you may refer to the following steps to install them.

- libsodium 1.0.18 https://doc.libsodium.org/

    ```bash
    mkdir temp
    cd temp
    wget https://download.libsodium.org/libsodium/releases/libsodium-1.0.18-stable.tar.gz
    tar -xzvf libsodium-1.0.18-stable.tar.gz 
    cd libsodium-stable/
    ./configure
    make && make check     
    make install
    ```

- GMP 6.2.0 https://gmplib.org/ (this version should match NTL)

    ```bash
    cd temp
    wget https://ftp.gnu.org/gnu/gmp/gmp-6.2.0.tar.xz
    tar -xf gmp-6.2.0.tar.xz
    cd gmp-6.2.0/
    ./configure
    make && make check     
    make install
    ```        

- NTL 11.5.1 https://libntl.org/: install this library following the documentation https://libntl.org/doc/tour-unix.html

    ```bash
    cd temp
    wget https://libntl.org/ntl-11.5.1.tar.gz
    tar -xzvf ntl-11.5.1.tar.gz
    cd ntl-11.5.1/src/
    ./configure
    make && make check     
    make install
    ```   

- TBB (https://en.wikipedia.org/wiki/Threading_Building_Blocks)

    ```bash
    sudo apt-get install libtbb-dev
    ```

- cmake version 3.24 or later: https://cmake.org/install/

    ```bash
    cd temp
    git clone https://github.com/Kitware/CMake.git
    cd CMake
    ./bootstrap
    make -j$(nproc)
    make install
    hash -r
    cmake --version
    ```       

- swig (https://github.com/swig/swig):

    ```bash 
    sudo apt-get update && sudo apt-get -y install swig
    ```

After finish installing the dependencies, can run `bash make.sh` to build the library. After that, can 
test if the build is successful by `cd build` and `./test`.

If use the library for python to invoke, can check the python sys.path by `python3 -m site` to get the path of 
python packages, and then cp the libraries to that path, for example:

   ```bash 
   cp -t /usr/local/lib/python3.8/dist-packages librisefl_crypto.a _risefl_interface.so risefl_interface.py
   ```

Then, can test if the library installation by `python3` and `import risefl_interface`. 


## Development Manual

- [Server side](doc/server.md)
- [Client side](doc/client.md)

Test script with explanations:
- [python/test.py](python/test.py) for a simulation of one iteration step of the system

Simulation experiments with three datasets as in the paper https://arxiv.org/abs/2311.15310 uses [a fork of FLSim](https://github.com/zhuyizheng/FLSim).