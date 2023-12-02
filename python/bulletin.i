%module di_zkp_interface_client

%include "std_vector.i"
%include "std_string.i"
%include "../include/bulletin.h"

namespace std {
        %template(VecFloat) vector<float>;
        %template(VecVecFloat) vector< vector<float> >;
        %template(PairSignKeys) pair<SignPubKey, SignPrvKey>;
        %template(VecSignPubKeys) vector<SignPubKey>;
        %template(VecSignPrvKeys) vector<SignPrvKey>;
};

%{
#include "../include/ristretto.h"
#include "../include/ristretto_vector.h"
#include "../include/zkp.h"
#include "../include/rist_fast_computation.h"
#include "../include/exception.h"
#include "../include/utils.h"
#include "../include/di_zkp_interface_common.h"
#include "../include/shamir.h"
#include "../include/di_zkp_interface_client.h"
#include "../include/random_avx/mersenne-twister-avx2.h"
#include "../include/bulletin.h"
%}


