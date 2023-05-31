#include "globals.h"

#ifndef JOEMATRIXCALC_LIN_ALG_H
#define JOEMATRIXCALC_LIN_ALG_H


namespace lin_alg{
    /** Computes the lie bracket mn-nm. */
    g::matrix bracket(g::matrix &m, g::matrix &n);
    std::vector< g::matrix > spanning_subsequence(std::vector< g::matrix > &matrices);
    std::vector< g::matrix > is_linearly_ind(std::vector< g::matrix > &matrices);
    std::vector< g::matrix > intersect(std::vector< g::matrix > &matrices1, std::vector< g::matrix > &matrices2);
    std::vector< g::matrix > symbolic_matrix_to_basis();
    /** Swaps the i-th and j-th rows of m in place.
     *  
     *  If m has k rows then we need 0<= i,j < k.
     */
    void swap_rows(g::matrix* m, int i, int j);  
    g::matrix gaussian_elimination_columnless(g::matrix const &matrices);
    g::exvector vectorize(g::matrix &m);

};

#endif //JOEMATRIXCALC_LIN_ALG_H