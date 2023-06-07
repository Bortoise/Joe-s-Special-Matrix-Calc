#include "globals.h"
#include <algorithm>

#ifndef JOEMATRIXCALC_LIN_ALG_H
#define JOEMATRIXCALC_LIN_ALG_H


namespace lin_alg{
    /** Computes the lie bracket mn-nm. */
    g::matrix bracket(g::matrix &m, g::matrix &n);
    std::vector< g::matrix > spanning_subsequence(std::vector< g::matrix > &matrices);
    std::vector< g::matrix > is_linearly_ind(std::vector< g::matrix > &matrices);
    std::vector< g::matrix > intersect(std::vector< g::matrix > &matrices1, std::vector< g::matrix > &matrices2);
    std::vector< g::matrix > symbolic_matrix_to_basis();

    /** 
     * @brief Swaps the i-th and j-th rows of m in place.
     *  
     * @param[out] m A pointer to a matrix.
     * @param[in] i The index of the first row to swap.
     * @param[in] j The index of the second row to swap.
     * 
     * @warning If m has k rows then we need 0<= i,j < k.
     * @note This function modifies the matrix m.
     */
    void swap_rows(g::matrix* m, int i, int j);
   
    /**
     * @brief Construct a the gaussian elimination of the given matrix.
     * 
     * @param matrix The matrix to be put into row echelon form.
     * 
     * @return The row echelon form of matrix.
     * 
     * @note This function does not modify the matrix matrix. This is division free
     *       gaussian elimination, performing only multiplications and addition of elements.
     */
    g::matrix gaussian_elimination(g::matrix const &matrix);

    /**
     * @brief Compute the rank of the matrix m by constructing the gaussian elimination
     *        of m and counting the number of non-zero rows.
     * 
     * @param m A matrix.
     * @return The rank of the matrix m.
     */
    int rank(g::matrix &m);

    /**
     * @brief Returns the nullspace of the given matrix.
     * 
     * @param m A matrix.
     * @return A basis of the nullspace of m.
     */
    std::vector< g::exvector > nullspace(g::matrix &m);

    g::exvector vectorize(g::matrix &m);
    g::matrix basis_to_vectorized_matrix(std::vector< g::matrix > &basis);

};

#endif //JOEMATRIXCALC_LIN_ALG_H