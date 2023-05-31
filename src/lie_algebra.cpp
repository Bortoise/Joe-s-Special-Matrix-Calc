#include "headers/lie_algebra.h"

// using namespace L;
namespace se = std::experimental;

lie_algebra::lie_algebra(std::vector< g::matrix > generators, bool _basis){
    if (generators[0].rows() != generators[0].cols()) {
        throw "Input matrices are not square"; // make exception for
    }
    sl_size = generators[0].rows();
    std::vector< g::matrix > new_basis(generators);
    if(!_basis) {
        int old_dim = 0;
        while (old_dim != new_basis.size()) {
            old_dim = static_cast<int>(new_basis.size());
            for (int i = 0; i < old_dim - 1; i++) {
                for (int j = i + 1; j < old_dim; j++) {
                    g::matrix& m = *new g::matrix();
                    m = lin_alg::bracket(new_basis[i], new_basis[j]); // Do NOT try to make this one line, if you do it doesn't work for some reason
                    new_basis.push_back(m);
                }
            }
            new_basis = lin_alg::spanning_subsequence(new_basis);
        }
    }

    // Sets the this.dim to the number of basis elements
    dim = new_basis.size();
    // Sets this.basis to the new_basis we either generated or were given
    basis = new_basis; //TODO: maybe this can be replaced by just working with this.basis in every step

    // Instantiates the fields to None
    derived_series = se::optional< std::vector< lie_algebra* > >();
    lower_central_series = se::optional< std::vector< lie_algebra* > >();
    normalizer = se::optional< lie_algebra* >();
    centralizer = se::optional< lie_algebra* >();

}
