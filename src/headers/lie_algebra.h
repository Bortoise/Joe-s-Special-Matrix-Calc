#include "globals.h"
#include "lin_alg.h"

#ifndef JOEMATRIXCALC_LIEALGEBRA_H
#define JOEMATRIXCALC_LIEALGEBRA_H



/** Lie subalgebra of sl(n) */
class lie_algebra {

    private:
        int sl_size; /** Which sl(n) you're working in */
        int dim; /** Dimension of the Lie algebra */
        std::vector< g::matrix > basis; /** Basis for sub algebra of sl(n).*/
        std::experimental::optional< std::vector< lie_algebra* > > derived_series;
        std::experimental::optional< std::vector< lie_algebra* > > lower_central_series;
        std::experimental::optional< lie_algebra* > normalizer;
        std::experimental::optional< lie_algebra* > centralizer;

    public:
        /** Constructor for lie_algebra.
         *
         * @param generators generators of your lie_algebra OR basis if you set basis=true
         * @param basis set basis=true if you are giving a basis in generators, instead of a generating set to save yourself computation time
         *
         * @note ensure all elements of generators are square matrices of the same size
         */
        lie_algebra(std::vector< g::matrix > generators, bool _basis=false ); //algebra from generators (set basis to false if going from set of generators)

        /** Returns a (shallow) copy of the basis. */
        std::vector<g::matrix> get_basis();
        /** Returns sl_size. */
        int get_sl_size();
        /** Returns dim = basis.size() */
        int get_dim();

        /** Computes the normalizer in sl(n), N_{sl(n)}(L),  stores it in this.normalizer. */
        lie_algebra compute_normalizer();
        /** Computes the centralizer in sl(n), C_{sl(n)}(L), and stores it in this.centralizer. */
        lie_algebra compute_centralizer();
        /** Computes the bracket [L, sl(n)] */
        lie_algebra bracket_with_sl();

        /** Returns the derived subalgebra [L, L]. */
        lie_algebra compute_derived_subalgebra();
        /** Returns the derived series L^n=[L^{n-1}, L]. */
        std::vector< lie_algebra* > compute_derived_series(); // Returns the derived series of the algebra of the whole Lie algebra until it stabilizes
        /** Returns the derived series L^(n)=[L^(n-1), L^(n-1)]. */
        std::vector< lie_algebra* > compute_lower_central_series(); // Returns the lower central series of the algebra of the whole Lie algebra until it stabilizes

        bool is_abelian(); // Checks if the algebra is abelian.
        bool is_solvable(); // Checks if the algebra is solvable.
        bool is_nilpotent(); // Checks if the algebra is nilpotent.

        bool equals(lie_algebra &other); // Checks if they are the same vector space. First checks the dimension. 
        bool contained_in(lie_algebra &other); // Checks if this is contained in the other.

        std::vector< g::matrix > extend_basis(std::vector< g::matrix > sl_basis); // Extends basis to all of sl(n) using sl_basis
};

//lie_algebra bracket_lie_algebras( lie_algebra &algebra1, lie_algebra &algebra2, int dim); // < this seems studpid, why do we give it dim, it should be able to tell that from the matrices it's given
//lie_algebra sl_basis(int dim);
//lie_algebra intersect(lie_algebra &lie_algebra1, lie_algebra &lie_algebra2);
// (index, sampleMatrices) nilpotent_index(); //TODO: figure out what this means
// subspace? intersect(A,B); //TODO: figure out what this means
// (int, int) minimum_rank_algorithm(, lie_algebra)


#endif // JOEMATRIXCALC_LIEALGEBRA_H