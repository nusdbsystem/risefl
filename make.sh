mkdir -p build && \
cmake -Bbuild -H. && \
make -C build -j 4 && \
cp build/_di_zkp_interface.so build/di_zkp_interface.py build/libdi_zkp_crypto.a python/