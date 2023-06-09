#include "headers/lie_algebra.h"

// using namespace L;
// namespace se = std::experimental;
lie_algebra::lie_algebra(std::vector<g::matrix> generators, bool _basis) {
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
                    new_basis.push_back(lin_alg::bracket(new_basis[i], new_basis[j]));
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
    derived_series = stdx::optional< std::vector< lie_algebra* > >();
    lower_central_series = stdx::optional< std::vector< lie_algebra* > >();
    normalizer = stdx::optional< lie_algebra* >();
    centralizer = stdx::optional< lie_algebra* >();

}

lie_algebra::~lie_algebra() { //TODO: figure out funky warnings in the destructor
    delete (&basis);
    try {
        std::vector<lie_algebra *> D = derived_series.value();
        for (auto & i : D) {
            delete i;
        }
        delete (&D);
    } catch (stdx::bad_optional_access) {}
    delete (&derived_series);
    try {
        std::vector<lie_algebra *> L = lower_central_series.value();
        for (auto & i : L) {
            delete i;
        }
    } catch (stdx::bad_optional_access) {}
    delete (&lower_central_series);
    try {
        lie_algebra* N = normalizer.value();
        delete(N);
    } catch(stdx::bad_optional_access) {}
    delete(&normalizer);
    try {
        lie_algebra* C = centralizer.value();
        delete(C);
    } catch(stdx::bad_optional_access) {}
    delete(&centralizer);
}
/*
std::vector<g::matrix> lie_algebra::get_basis() {
    std::vector< g::matrix > v (this->basis);
    return v;
}

int lie_algebra::get_sl_size() const { return this->sl_size; }

int lie_algebra::get_dim() const { return this->dim; }

lie_algebra* lie_algebra::get_sl(int n) {
    if (lie_algebra::sl.has_value()) {
        return lie_algebra::sl.value();
    }
    std::vector<g::matrix> basis = std::vector<g::matrix>();
    g::matrix m = g::matrix(n,n);
    for (int i = 0; i < n-1; i++) {
        for (int j = i+1; j < n; j++) {
            m(i,j) = 1;
            basis.push_back(m);
            m(i,j) = 0;
            m(j,i) = 1;
            basis.push_back(m);
            m(j,i) = 0;
        }
    }
    m(0,0) = 1;
    for (int i = 1; i < n-1; i++) {
        m(i,i) = -1;
        basis.push_back(m);
        m(i,i) = 0;
    }
    lie_algebra* out = new lie_algebra(basis, true);
    lie_algebra::sl = stdx::optional(out);
    return out;
}

lie_algebra* lie_algebra::compute_centralizer() {
    // Let L have basis e_i.
    lie_algebra* out = compute_centralizer_element(this->basis[0], get_sl(this->get_sl_size())); // Sets out=C_{sl(n)}(e_1)
    for (int i = 1; i < this->dim; i++) {
        out = compute_centralizer_element(this->basis[i], out); // Updates C_{sl(n)}(e_1,...,e_{i+1}) -> C_{C_{sl(n)}(e_1,...,e_i)}(e_{i+1})
    }
    this->centralizer = stdx::optional< lie_algebra* >(out); // Records the centralizer
    return out;
}

lie_algebra* lie_algebra::compute_normalizer() {
    // Let L have basis {e_i, i<=r}. Set M_0 = sl(n) and M_{i+1}=N(x,L,M_i), where N(x,L,M) is the elements y of M such that ad(y) x in L.
    // It is clear that N(L)=\bigcap_{j<= r} N(x,L,sl(n))=M_r
    lie_algebra* out = this->compute_normalizer_element(this->basis[0], get_sl(this->get_sl_size())); // Sets out=N_{sl(n)}(e_1,L)=M_1
    for (int i = 1; i < this->dim; i++) {
        out = this->compute_normalizer_element(this->basis[i], out); // Updates out -> M_{i+1}
    }
    this->normalizer = stdx::optional< lie_algebra* >(out); // Records the normalizer
    return out;
}

lie_algebra* lie_algebra::bracket_with_sl() { //Needs to be changed if we decide to go the other way about making sl a static member
    lie_algebra* sl = get_sl(this->get_sl_size());
    return bracket_lie_algebras(this, sl);
}

lie_algebra* lie_algebra::compute_derived_subalgebra() {
    if(this->derived_series.has_value()) {
        return this->derived_series.value()[0];
    }
    lie_algebra* m = bracket_lie_algebras(this, this);
    return m;
}

std::vector< lie_algebra* > lie_algebra::compute_derived_series() {
    std::vector< lie_algebra* > derived_series_out = std::vector< lie_algebra* >();
    int old_dim = this->dim;
    lie_algebra* D = this->compute_derived_subalgebra(); // Sets L^1 = [L,L]
    
    while(D->get_dim() != old_dim) { // Terminates when L^n=L^{n+1}
        derived_series_out.push_back(D);
        old_dim = D->get_dim(); // Updates old_dim -> dim L^n
        D = D->compute_derived_subalgebra(); // Updates D -> L^{n+1}
    }
    
    if(derived_series_out.empty()) { // This deals with the case [L,L]=L
        derived_series_out.push_back(D);
    }

    this->derived_series = stdx::optional<std::vector< lie_algebra* >>(derived_series_out); // Records the derived series 
    return std::vector(derived_series_out);
}

std::vector< lie_algebra* > lie_algebra::compute_lower_central_series() {
    std::vector< lie_algebra* > lower_central_series_out = std::vector< lie_algebra* >();
    int old_dim = this->dim;
    lie_algebra* C = this->compute_derived_subalgebra(); // Sets L^(1) = [L,L]
    
    while(C->get_dim() != old_dim) { // Terminates when L^(n)=L^(n+1)
        lower_central_series_out.push_back(C);
        old_dim = C->get_dim(); // Updates old_dim -> dim L^(n)
        C = bracket_lie_algebras(C, this); // Updates C -> L^(n+1)
    }

    if(lower_central_series_out.empty()) { // This deals with the case [L,L]=L
        lower_central_series_out.push_back(C);
    }
    
    this->derived_series = stdx::optional<std::vector< lie_algebra* >>(lower_central_series_out); // Records the lower central series 
    return std::vector(lower_central_series_out);
}

bool lie_algebra::is_abelian() { return this->compute_derived_subalgebra()->get_dim() == 0; }


bool lie_algebra::equals(lie_algebra* N) {
    if (N->get_dim() != this->get_dim()) {
        return false;
    }
    return N->contains(this) && this->contains(N);
}

bool lie_algebra::contains(lie_algebra* N) {
    if (N->get_dim() > this->get_dim()) {
        return false;
    }
    std::vector< g::matrix > vectors = std::vector< g::matrix >();

    // Let e_1,...,e_n be a basis of L
    std::vector< g::matrix > basis_1 = this->get_basis();
    // Let f_1,...,f_m be a basis of N
    std::vector< g::matrix > basis_2 = N->get_basis();
    // Makes the list S=(e_1,...,e_n,f_1,...,f_m)
    vectors.insert(vectors.end(), basis_1.begin(), basis_1.end());
    vectors.insert(vectors.end(), basis_2.begin(), basis_2.end());

    // Makes a matrix whose rows are the basis elements of L and N.
    g::matrix M = lin_alg::basis_to_vectorized_matrix(vectors);

    // N is a subset of L iff dim(L)=dim(L+N), but that is to say rank(e_1,...,e_n)
    // = rank(S). Taking transposes this is what we actually verify.
    return this->get_dim() == lin_alg::rank(M);
}


bool lie_algebra::contains_element(g::matrix x) {
    std::vector< g::matrix > N_basis = std::vector< g::matrix >();
    N_basis.push_back(x);
    lie_algebra N = {N_basis, true};
    return this->contains(&N);
}

std::vector< g::matrix > lie_algebra::extend_basis(lie_algebra* M) {
    std::vector< g::matrix > vec_list = std::vector< g::matrix >(this->basis);
    for(g::matrix v : M->basis) {
        vec_list.push_back(v);
    }
    return lin_alg::spanning_subsequence(vec_list);
}
*/

/**
 * TODO: FIGURE OUT THE PROPER INPUT REQUIREMENTS OR MAKE THIS A PRIVATE MEMBER OF lie_algebra
 */
// TODO: Fix this algorithm (need to take the kernel of ad(x):N->sl.
//lie_algebra compute_centralizer_element(g::matrix x, lie_algebra N) {
//    std::vector< g::matrix > basis = std::vector< g::matrix >();
//    std::vector< g::matrix > old_basis = N.get_basis();
//    for (int i = 0; i < old_basis.size(); i++) { // Forms the matrix of ad(x) acting on N
//        old_basis[i] = lin_alg::bracket(x, old_basis[i]);
//    }
//    g::matrix m = lin_alg::basis_to_vectorized_matrix(old_basis);
//    std::vector< g::exvector > vec_null_basis = lin_alg::nullspace(m);
//    for (g::exvector v : vec_null_basis) {
//        basis.push_back(g::matrix(N.get_sl_size(), N.get_sl_size(), v)); //TODO: why you like this ginac
//    }
//    return lie_algebra(basis, true);
//}
/*

lie_algebra* bracket_lie_algebras(lie_algebra &algebra1, lie_algebra  &algebra2) {
    std::vector< g::matrix > basis = std::vector< g::matrix >();
    for (g::matrix v1 : algebra1.get_basis()) {
        for (g::matrix v2 : algebra2.get_basis()) {
            basis.push_back(lin_alg::bracket(v1, v2));
        }
    }
    basis = lin_alg::spanning_subsequence(basis);
    lie_algebra* out = new lie_algebra(basis, true);
    return out;
}
*/