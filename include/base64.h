//
// Created by wuyuncheng on 4/24/23.
//

#ifndef FL_PROB_ZKP_BASE64_H
#define FL_PROB_ZKP_BASE64_H

#include <string>
#include <vector>

typedef unsigned char BYTE;

std::string base64_encode(BYTE const *buf, unsigned int bufLen);

inline std::string base64_encode(std::vector<BYTE> r) {
    return base64_encode(r.data(), r.size());
}

std::vector<BYTE> base64_decode(std::string const &);


#endif //FL_PROB_ZKP_BASE64_H
