#include "globals.h"

#ifndef JOEMATRIXCALC_LIN_ALG_H
#define JOEMATRIXCALC_LIN_ALG_H


namespace lin_alg{
    g::matrix bracket(g::matrix &m, g::matrix &n);
    std::vector< g::matrix* > spanning_subsequence(std::vector< g::matrix* > &matrices);
    std::vector< g::matrix* > is_linearly_ind(std::vector< g::matrix* > &matrices);
    std::vector< g::matrix* > intersect(std::vector< g::matrix* > &matrices1, std::vector< g::matrix* > &matrices2);
    std::vector< g::matrix* > symbolic_matrix_to_basis();
    
    g::exvector vectorize(g::matrix &m);

};

#endif //JOEMATRIXCALC_LIN_ALG_H