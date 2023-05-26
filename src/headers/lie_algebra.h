#include "globals.h"

class lie_algebra {
    private:
        std::vector<g::matrix> basis;
        std::optional< std::vector< lie_algebra > > derived_series;
        std::optional< std::vector< lie_aglebra > > lower_central_series;
        std::optional< lie_algebra > normalizer;
        std::optional< lie_algebra > centralizer;

    public:
        int n; // sl(n) you're working in
        int dim; // Dimension of the Lie algebra
        
        lie_algebra(std::vector< g::matrix > &generators, bool basis=true); //algebra from generators (set basis to false if going from set of generators)

        std::vector<g::matrix> get_basis();
        
        lie_algebra normalizer();
        lie_algebra centralizer();
        lie_algebra bracket_with_sl( std::vector< g::matrix > &generators); // Returns [L, sl(n)]
        
        lie_algebra derived_subalgebra(); // Returns the derived subalgebra [L, L] 
        std::vector< lie_algebra > derived_series(); // Returns the derived series of the algebra of the whole Lie algebra until it stabilizes
        std::vector< lie_algebra > lower_central_series(); // Returns the lower central series of the algebra of the whole Lie algebra until it stabilizes
        
        bool is_abelian(); // Checks if the algebra is abelian.
        bool is_solvable(); // Checks if the algebra is solvable.
        bool is_nilpotent(); // Checks if the algebra is nilpotent.

        bool equals(lie_algebra &other); // Checks if they are the same vector space.
        bool contained_in(lie_algebra &other); // Checks if this is contained in the other.

        std::vector< g::matrix > extend_basis(std::vector< g::matrix > sl_basis); // Extends basis to all of sl(n) using sl_basis
}

g::matrix bracket(g::matrix &m, g::matrix &n);
lie_algebra algebra_from_generators( std::vector< g::matrix > &generators); //TODO: handle matrix sizes
lie_algebra bracket_lie_algebras( lie_algebra &algebra1, lie_algebra &algebra2, int dim); // < this seems studpid, why do we give it dim, it should be able to tell that from the matrices it's given
lie_algebra sl_basis(int dim);
std::vector< g::matrix > is_linearly_ind(std::vector< g::matrix > &matrices);