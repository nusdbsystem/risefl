# API for RiseFL, the Robust Secure FL system, built on C/C++ core

Install requirements:

- libsodium 1.0.18 https://doc.libsodium.org/
- GMP https://gmplib.org/
- NTL https://libntl.org/: install this library following the documentation https://libntl.org/doc/tour-unix.html
- TBB (https://en.wikipedia.org/wiki/Threading_Building_Blocks): `sudo apt-get install libtbb-dev`
- cmake version 3.24 or later: https://cmake.org/install/
- swig (https://github.com/swig/swig): `sudo apt-get update && sudo apt-get -y install swig`

How to install:
- Run `./make.sh`

Manual:
- [Server side](doc/server.md)
- [Client side](doc/client.md)

Test script with explanations:
- [python/test.py](python/test.py) for a simulation of one iteration step of the system

Simulation experiments with three datasets as in the paper https://arxiv.org/abs/2311.15310 uses [a fork of FLSim](https://github.com/zhuyizheng/FLSim).