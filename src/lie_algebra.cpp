#include "headers/lie_algebra.h"

// using namespace L;
namespace se = std::experimental;

lie_algebra::lie_algebra(std::vector< g::matrix > generators, bool _basis){
    //TODO: figure out how to make this->sl / if we want to
    if (generators[0].rows() != generators[0].cols()) {
        throw "Input matrices are not square"; // make exception for
    }
    sl_size = generators[0].rows();
    std::vector< g::matrix > new_basis = lin_alg::spanning_subsequence(generators);
    if(!_basis) {
        int old_dim = 0;
        while (old_dim != static_cast<int>(new_basis.size())) {
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

lie_algebra::~lie_algebra() { //TODO: figure out funky warnings in the destructor
    delete(&basis);
    try {
        std::vector< lie_algebra* > D = derived_series.value();
        for(int i = 0; i < D.size(); i++) {
            delete(D[i]);
        }
        delete(&D);
    } catch(se::bad_optional_access) {}
    delete(&derived_series);
    try {
        std::vector< lie_algebra* > L = lower_central_series.value();
        for(int i = 0; i < L.size(); i++) {
            delete(L[i]);
        }
    } catch(se::bad_optional_access) {}
    delete(&lower_central_series);
    try {
        lie_algebra* N = normalizer.value();
        delete(N);
    } catch(se::bad_optional_access) {}
    delete(&normalizer);
    try {
        lie_algebra* C = centralizer.value();
        delete(C);
    } catch(se::bad_optional_access) {}
    delete(&centralizer);
}

int lie_algebra::get_sl_size() { return this->sl_size; }

int lie_algebra::get_dim() { return this->dim; }

std::vector<g::matrix> lie_algebra::get_basis() { 
    std::vector< g::matrix > v (this->basis);
    return v;
}

lie_algebra lie_algebra::compute_centralizer() {
    // Let L have basis e_i.
    lie_algebra out = compute_centralizer_element(this->basis[0], this->sl); // Sets out=C_{sl(n)}(e_1)
    for (int i = 1; i < this->dim; i++) { //TODO, only check on generating set? I.e. does this work mathematically and how to change the code to allow this
        out = compute_centralizer_element(this->basis[i], out); // Updates C_{sl(n)}(e_1,...,e_{i+1}) -> C_{C_{sl(n)}(e_1,...,e_i)}(e_{i+1})
    }
    this->centralizer = se::optional< lie_algebra* >(&out); // Records the centralizer
    return out;
}

lie_algebra lie_algebra::bracket_with_sl() { //Needs to be changed if we decide to go the other way about making sl a static member
    return bracket_lie_algebras(*this, this->sl);
}

lie_algebra lie_algebra::compute_derived_subalgebra() {
    lie_algebra m = bracket_lie_algebras(*this, *this);
    return m;
}

std::vector< lie_algebra* > lie_algebra::compute_derived_series() {
    std::vector< lie_algebra* > derived_series_out = std::vector< lie_algebra* >();
    int old_dim = this->dim;
    lie_algebra D = this->compute_derived_subalgebra(); // Sets L^1 = [L,L]
    
    while(D.get_dim() != old_dim) { // Terminates when L^n=L^{n+1}
        derived_series_out.push_back(&D); 
        old_dim = D.get_dim(); // Updates old_dim -> dim L^n
        D = D.compute_derived_subalgebra(); // Updates D -> L^{n+1}
    }
    
    if(derived_series_out.empty()) { // This deals with the case [L,L]=L
        derived_series_out.push_back(&D);
    }

    this->derived_series = se::optional<std::vector< lie_algebra* >>(derived_series_out); // Records the derived series 
    return std::vector(derived_series_out);
}

std::vector< lie_algebra* > lie_algebra::compute_lower_central_series() {
    std::vector< lie_algebra* > lower_central_series_out = std::vector< lie_algebra* >();
    int old_dim = this->dim;
    lie_algebra C = this->compute_derived_subalgebra(); // Sets L^(1) = [L,L]
    
    while(C.get_dim() != old_dim) { // Terminates when L^(n)=L^(n+1)
        lower_central_series_out.push_back(&C); 
        old_dim = C.get_dim(); // Updates old_dim -> dim L^(n)
        C = bracket_lie_algebras(C, *this); // Updates C -> L^(n+1)
    }

    if(lower_central_series_out.empty()) { // This deals with the case [L,L]=L
        lower_central_series_out.push_back(&C);
    }
    
    this->derived_series = se::optional<std::vector< lie_algebra* >>(lower_central_series_out); // Records the lower central series 
    return std::vector(lower_central_series_out);
}


bool lie_algebra::is_abelian() {
    lie_algebra N = this->compute_derived_subalgebra();
    int dim = N.get_dim();
    delete(&N);
    return dim;
}




lie_algebra compute_centralizer_element(g::matrix x, lie_algebra N) {
    std::vector< g::matrix > basis = std::vector< g::matrix >();
    std::vector< g::matrix > old_basis = N.get_basis();
    for (std::vector< g::matrix >::iterator iter = old_basis.begin(); iter < old_basis.end(); iter++) {
        if (lin_alg::bracket(*iter, x).is_zero_matrix()) {
            basis.push_back(*iter);
        }
    }
    return lie_algebra(basis, true);
}