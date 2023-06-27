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
    void print_exvector(g::exvector &v);
    void print_column_matrices(mat_vec &v);
    void print_matrices(mat_vec m);

    
    /** Tests if two matrices are equal. */
    bool matrix_eq(g::matrix &m, g::matrix &n);

    
    /** Tests if two exvectors are equal. */
    bool exvector_eq(g::exvector &v, g::exvector &w);
}

#endif // JOEMATRIXCALC_UTILS_H
