#include "globals.h"

#ifndef JOEMATRIXCALC_LIEALGEBRA_H
#define JOEMATRIXCALC_LIEALGEBRA_H

class lie_algebra {

    private:
        int sl_size; /** Which sl(n) you're working in */
        int dim; /** Dimension of the Lie algebra */
        std::vector< g::matrix* > basis; /** Basis for sub algebra of sl(n).*/
        std::optional< std::vector< lie_algebra* > > derived_series;
        std::optional< std::vector< lie_algebra* > > lower_central_series;
        std::optional< lie_algebra* > normalizer;
        std::optional< lie_algebra* > centralizer;

    public:
        /** Lie subalgebra of sl(n). */
        lie_algebra(std::vector< g::matrix > &generators, bool basis=true); //algebra from generators (set basis to false if going from set of generators)

        std::vector<g::matrix> get_basis();
        int get_sl_size();
        int get_dim();
        
        lie_algebra compute_normalizer(); // Computes the normalizer of the algebra in sl(n)
        lie_algebra compute_centralizer(); // Computes the centralizer of the algebra in sl(n)
        lie_algebra bracket_with_sl( std::vector< g::matrix > &generators); // Returns [L, sl(n)]
        
        lie_algebra compute_derived_subalgebra(); // Returns the derived subalgebra [L, L] 
        std::vector< lie_algebra* > compute_derived_series(); // Returns the derived series of the algebra of the whole Lie algebra until it stabilizes
        std::vector< lie_algebra* > compute_lower_central_series(); // Returns the lower central series of the algebra of the whole Lie algebra until it stabilizes

        bool is_abelian(); // Checks if the algebra is abelian.
        bool is_solvable(); // Checks if the algebra is solvable.
        bool is_nilpotent(); // Checks if the algebra is nilpotent.

        bool equals(lie_algebra &other); // Checks if they are the same vector space. First checks the dimension. 
        bool contained_in(lie_algebra &other); // Checks if this is contained in the other.

        std::vector< g::matrix > extend_basis(std::vector< g::matrix* > sl_basis); // Extends basis to all of sl(n) using sl_basis
};

//g::matrix bracket(g::matrix &m, g::matrix &n);
lie_algebra bracket_lie_algebras( lie_algebra &algebra1, lie_algebra &algebra2, int dim); // < this seems studpid, why do we give it dim, it should be able to tell that from the matrices it's given
lie_algebra sl_basis(int dim);
std::vector< g::matrix* > is_linearly_ind(std::vector< g::matrix* > &matrices);
lie_algebra symbolic_matrix_to_basis();
std::vector< g::matrix* > intersect(std::vector< g::matrix* > &matrices1, std::vector< g::matrix* > &matrices2);
lie_algebra intersect(lie_algebra &lie_algebra1, lie_algebra &lie_algebra2);
// (index, sampleMatrices) nilpotent_index(); //TODO: figure out what this means
// subspace? intersect(A,B); //TODO: figure out what this means
// (int, int) minimum_rank_algorithm(, lie_algebra)


#endif // JOEMATRIXCALC_LIEALGEBRA_H