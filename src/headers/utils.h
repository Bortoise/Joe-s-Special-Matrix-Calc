//
// Created by cole on 15/05/23.
//

#include "globals.h"

#ifndef  JOEMATRIXCALC_UTILS_H
#define JOEMATRIXCALC_UTILS_H

namespace g = GiNaC;

namespace utils{
    void print_matrix(g::matrix &m);
    void print_exvectors(std::vector< g::exvector > &v);
    void print_matrices(std::vector< g::matrix > m);
}

#endif // JOEMATRIXCALC_UTILS_H
