%module risefl_interface_common

%include "std_vector.i"
%include "std_string.i"
%include "../include/risefl_interface_common.h"
%{
#include "../include/ristretto.h"
#include "../include/ristretto_vector.h"
#include "../include/zkp.h"
#include "../include/rist_fast_computation.h"
#include "../include/exception.h"
#include "../include/utils.h"
#include "../include/risefl_interface_common.h"
#include "../include/shamir.h"
#include "../include/random_avx/mersenne-twister-avx2.h"
%}

