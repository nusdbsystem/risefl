cmake_minimum_required(VERSION 3.24)
project(di_zkp_crypto)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_C_FLAGS "-Ofast -pthread -lpthread -DHAVE_TI_MODE=1 -march=native")
set(CMAKE_CXX_FLAGS "-Ofast -pthread -lpthread -fPIC -march=native")

message(STATUS "Project name is ${PROJECT_NAME}")

add_library(di_zkp_crypto include/ristretto.h include/ristretto_vector.h shamir.cpp
        include/shamir.h zkp.cpp include/zkp.h ristretto.cpp di_zkp_interface_server.cpp
        include/di_zkp_interface_server.h di_zkp_interface_client.cpp
        include/di_zkp_interface_client.h di_zkp_interface_common.cpp include/di_zkp_interface_common.h
        ristretto_vector.cpp main.cpp include/sodium/ed25519_ref10.c rist_fast_computation.cpp
        exception.cpp include/exception.h bench/simulate.h bench/simulate.cpp
        mpz_conversion.cpp include/mpz_conversion.h utils.cpp include/utils.h include/base64.h base64.cpp
        include/bulletin.h bulletin.cpp
        include/random_avx/mersenne-twister-avx2.h include/random_avx/mersenne-twister-avx2.cpp)
target_link_libraries(di_zkp_crypto sodium pthread ntl tbb gmpxx gmp m)

# Add install target
install(TARGETS ${PROJECT_NAME} di_zkp_crypto
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static
        PUBLIC_HEADER DESTINATION include)

# For testing the executable micro-benchmark
add_executable(test include/ristretto.h include/ristretto_vector.h shamir.cpp
        include/shamir.h zkp.cpp include/zkp.h ristretto.cpp di_zkp_interface_server.cpp
        include/di_zkp_interface_server.h di_zkp_interface_client.cpp
        include/di_zkp_interface_client.h di_zkp_interface_common.cpp include/di_zkp_interface_common.h
        ristretto_vector.cpp main.cpp include/sodium/ed25519_ref10.c rist_fast_computation.cpp
        exception.cpp include/exception.h bench/simulate.h bench/simulate.cpp
        mpz_conversion.cpp include/mpz_conversion.h utils.cpp include/utils.h include/base64.h base64.cpp
        include/bulletin.h bulletin.cpp
        include/random_avx/mersenne-twister-avx2.h include/random_avx/mersenne-twister-avx2.cpp)

target_link_libraries(test sodium pthread ntl tbb gmpxx gmp m)


# Set up swig
find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

# Include python
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})

message(STATUS "PROJECT_SOURCE_DIR name is ${PROJECT_SOURCE_DIR}")

set(CMAKE_SWIG_FLAGS "")
set_source_files_properties(python/di_zkp_interface.i PROPERTIES CPLUSPLUS ON)
include_directories(${PROJECT_SOURCE_DIR}/include)

# generate
swig_add_library(di_zkp_interface LANGUAGE python SOURCES python/di_zkp_interface.i)

swig_link_libraries(di_zkp_interface di_zkp_crypto ${PYTHON_LIBRARIES})
