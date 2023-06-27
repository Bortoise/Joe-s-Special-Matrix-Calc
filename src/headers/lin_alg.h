#include "globals.h"
#include <algorithm>

#ifndef JOEMATRIXCALC_LIN_ALG_H
#define JOEMATRIXCALC_LIN_ALG_H


namespace lin_alg{
    /** Computes the lie bracket mn-nm. */
    g::matrix bracket(g::matrix &m, g::matrix &n);

     /**
     * @brief Returns a vector representing a subsequence of matrices 
     * 
     * @param matrices A list of matrices
     * @return basis A basis for span(matrices) which is also a subsequence of (matrices)
     */
    mat_vec spanning_subsequence(mat_vec &matrices);
    
    mat_vec is_linearly_ind(mat_vec &matrices);
    mat_vec intersect(mat_vec &matrices1, mat_vec &matrices2);
    mat_vec symbolic_matrix_to_basis();

    /**
     * @brief Converts a matrix to a vector in a basis whose elements are matrices
     * 
     * @param m The matrix to be converted to a vector.
     * @param basis A basis of matrices.
     * 
     * @return A vector representing the matrix with respect to the given basis
     */
    g::exvector matrix_to_vector_in_basis(g::matrix &m, mat_vec &basis);

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
    void swap_rows(g::matrix &m, int i, int j);
   
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
     * 
     * @note currently we use .normal() at all times this may or may not be slow, we need to test this.
     * @warning division = true does not work. It does not simplify rational expressions to 0. We know how to fix this (call GiNaC .normal() on all expression that may contain rational expressions), but have not gotten around to it yet.
     */
    mat_vec nullspace(g::matrix &m, bool division = true);

    g::exvector vectorize(g::matrix &m);

    /**
     * @brief Returns the matrix representation of the given vector.
     * 
     * @param v Exvector to be converted into the matrix.
     * @param r Number of rows of the matrix
     * @param c Number of columns of the matrix
     * @param rowcol If false, then inserts the vector row by row. If true, then inserts the vector column by column.
     * @return Matrix representation of the vector v.
     * 
     * @note The matrix is constructed column by column.
     */
    g::matrix matricize(g::exvector &v, unsigned int r, unsigned int c, bool rowcol = false);
    g::matrix matricize(g::matrix &v, unsigned int r, unsigned int c, bool rowcol = false);
    

    /**
    * @brief Returns a vector of the coefficients of m with respect to the standard basis of sl(n).
    *
    * @param m An element of sl(n)
    * @return A vector of the coefficients of m with respect to the standard basis of sl(n).
    */
    g::exvector sl_ize(g::matrix m, int n);

    /**
     * @brief Returns the matrix represented by the given vector
     * 
     * @param v The vector of expressions to be converted to a matrix
     * @param basis The list of matrices to multiply the vector against
     * @return The matrix represented by the vector v in the given basis
     */
    g::matrix vector_to_matrix(g::exvector &v, mat_vec &basis);
    g::matrix vector_to_matrix(g::matrix &v, mat_vec &basis);

    /** Returns the matrix with the given matrices vectorized then put in as rows.
     * //TODO: finish docstring
     * @param basis
     * @return
     */
    g::matrix basis_to_vectorized_matrix(mat_vec &basis);

    /** Computes the tr(ab) assuming that the product ab makes sense and produces a square matrix*/
    g::ex prod_trace(g::matrix a, g::matrix b);


    /** Tests if the subspace spanned by A contains the subspace spanned by B, assuming A is linearly independant
    *
    */
    bool contains(std::vector< g::exvector >* A, std::vector< g::exvector >* B);

    /**
     * @brief Checks if A = B
     * 
     * @param A basis of A
     * @param B basis of B
     * @return true if A = B as subspaces
     * @return false if A != B as subspaces
     */
    bool equals(std::vector< g::exvector >* A, std::vector< g::exvector >* B);
    bool equals(mat_vec* A, mat_vec* B);
    bool equals(mat_vec* A, std::vector< g::exvector >* B);
};

#endif //JOEMATRIXCALC_LIN_ALG_H