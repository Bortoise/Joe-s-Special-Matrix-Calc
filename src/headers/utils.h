//
// Created by cole on 15/05/23.
//

#include "globals.h"
#include "lie_algebra.h"

#ifndef  JOEMATRIXCALC_UTILS_H
#define JOEMATRIXCALC_UTILS_H

namespace g = GiNaC;

namespace utils{
    void print_matrix(g::matrix &m);
    void print_exvectors(std::vector< g::exvector > &v);
}

#endif // JOEMATRIXCALC_UTILS_H
