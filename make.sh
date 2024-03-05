mkdir -p build && \
cmake -Bbuild -H. && \
make -C build -j 4 && \
cp build/_risefl_interface.so build/risefl_interface.py build/librisefl_crypto.a python/