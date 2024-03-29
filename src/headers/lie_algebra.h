#include "globals.h"
#include "lin_alg.h"
#include "utils.h"


#ifndef JOEMATRIXCALC_LIEALGEBRA_H
#define JOEMATRIXCALC_LIEALGEBRA_H



/** Lie subalgebra of sl(n) */
class lie_algebra { // We always refer the lie algebra we are working with L.

    private:
        int sl_size;
        int dim; /** Dimension of the lie algebra, L */
        mat_vec basis; /** Basis for sub algebra of sl(n).*/
        stdx::optional< std::vector< lie_algebra* > > derived_series;
        stdx::optional< std::vector< lie_algebra* > > lower_central_series;
        stdx::optional< lie_algebra* > normalizer;
        stdx::optional< lie_algebra* > centralizer;
        inline static stdx::optional< lie_algebra* > sl = stdx::optional< lie_algebra* >(); ;

    public:
        /** Constructor for lie_algebra.
         *
         * @param generators generators of your lie_algebra OR basis if you set basis=true
         * @param basis set basis=true if you are giving a basis in generators, instead of a generating set to save yourself computation time
         *
         * @note ensure all elements of generators are square matrices of the same size
         */
        lie_algebra(mat_vec generators, bool _basis=false ); //algebra from generators (set basis to false if going from set of generators)
        /** DESTROYS the lie algebra after.*/
        ~lie_algebra();

        /** Returns a (shallow) copy of the basis. */
        mat_vec get_basis();
        /** Returns sl_size. */
        [[nodiscard]] int get_sl_size() const;
        /** Returns dim = basis.size() */
        [[nodiscard]] int get_dim() const;

        static lie_algebra* get_sl(int n);

        /** Computes the normalizer in sl(n), N_{sl(n)}(L),  stores it in this->normalizer. */
        lie_algebra* compute_normalizer();

        /** Computes subset of M whose adjoint send x into L
         *
         * @param x A sl_size by sl_size matrix.
         * @param M A lie algebra contains L.
         */
        mat_vec compute_normalizer_element(g::matrix x, mat_vec M);

        /** Computes the centralizer in sl(n), C_{sl(n)}(L), and stores it in this->centralizer. */
        lie_algebra* compute_centralizer();
        /** Computes the bracket [L, sl(n)] */
        lie_algebra* bracket_with_sl();

        /** Returns the derived subalgebra [L, L]. */
        lie_algebra* compute_derived_subalgebra();
        /** Returns the derived series L^n=[L^{n-1}, L]. 
         *
         *  @note To save on effort, space, and a bit of computation time we do not fill 
         *  in the derived series of elements of the derived series. Do not bother to try and
         *  investigate the derived series of L^n, it is for your own good. 
         */
        std::vector< lie_algebra* > compute_derived_series(); // Returns the derived series of the algebra of the whole Lie algebra until it stabilizes
        /** Returns the lower central series series L^(n)=[L^(n-1), L^(n-1)] and stores it in this->centralizer. */
        std::vector< lie_algebra* > compute_lower_central_series(); // Returns the lower central series of the algebra of the whole Lie algebra until it stabilizes

        std::pair<int, int> nilpotent_indices(); // out[0] is length of the derived series, out[1] is the length of the lower central series

        bool is_abelian(); // Checks if the algebra is abelian.
        bool is_solvable(); // Checks if the algebra is solvable.
        bool is_nilpotent(); // Checks if the algebra is nilpotent.

        /** Returns true if N and L are equal subalgebras of sl(n). */
        bool equals(lie_algebra* N);
        /** Returns true if N is contained in L. */
        bool contains(lie_algebra* N); // Checks if this is contained in the other.
        /** Returns true if x is an element of L. */
        bool contains_element(g::matrix x);

        /** Extends the default basis for L to a basis for all of M, assuming L is contained in M.
         *
         * @param M A lie algebra containing L.
         * @return A basis of M extending the basis for L.
         */
        mat_vec extend_basis(lie_algebra* M);

        /** Computes the minimum rank of elements of L
        * @warning This is not currently working whatsoever
        */
        int min_rank();
        int max_rank();
};


/**
 * Computes the centralizer of x in N, C_N(x).
 *
 * @param x A sl_size by sl_size matrix.
 * @param N A lie-subalgebra of sl(sl_size).
 * @return Returns C_N(x)
 */
mat_vec compute_centralizer_element(g::matrix x, mat_vec N);
/** Returns the bracket of two lie algebras of sl(n). */
lie_algebra* bracket_lie_algebras(lie_algebra* algebra1, lie_algebra* algebra2);
//lie_algebra sl_basis(int dim);
//lie_algebra intersect(lie_algebra &lie_algebra1, lie_algebra &lie_algebra2);
// (index, sampleMatrices) nilpotent_index(); //TODO: figure out what this means
// subspace? intersect(A,B); //TODO: figure out what this means
// (int, int) minimum_rank_algorithm(, lie_algebra)


#endif // JOEMATRIXCALC_LIEALGEBRA_H