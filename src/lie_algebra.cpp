#include "headers/lie_algebra.h"

// using namespace L;
// namespace se = std::experimental;
lie_algebra::lie_algebra(std::vector<g::matrix> generators, bool _basis) {
    sl_size = 0;
    std::vector<g::matrix> new_basis = {};
    if(!generators.empty()) {
        if (generators[0].rows() != generators[0].cols()) {
            throw "Input matrices are not square"; // make exception for
        }
        sl_size = generators[0].rows();
        new_basis = lin_alg::spanning_subsequence(generators);
        if (!_basis) {
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
    try {
        std::vector<lie_algebra *> D = derived_series.value();
        for (auto & i : D) {
            delete i;
        }
        delete (&D);
    } catch (stdx::bad_optional_access) {}
    try {
        std::vector<lie_algebra *> L = lower_central_series.value();
        for (auto & i : L) {
            delete i;
        }
    } catch (stdx::bad_optional_access) {}
    try {
        lie_algebra* N = normalizer.value();
        delete(N);
    } catch(stdx::bad_optional_access) {}
    try {
        lie_algebra* C = centralizer.value();
        delete(C);
    } catch(stdx::bad_optional_access) {}
}

std::vector<g::matrix> lie_algebra::get_basis() {
    std::vector< g::matrix > v (this->basis);
    return v;
}

int lie_algebra::get_sl_size() const { return this->sl_size; }

int lie_algebra::get_dim() const { return this->dim; }

lie_algebra* lie_algebra::get_sl(int n) {
    if (lie_algebra::sl.has_value()) {
        if (lie_algebra::sl.value()->get_sl_size() == n) {
            return lie_algebra::sl.value();
        }
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
    for (int i = 1; i < n; i++) {
        m(i,i) = -1;
        basis.push_back(m);
        m(i,i) = 0;
    }
    lie_algebra* out = new lie_algebra(basis, true);
    lie_algebra::sl = stdx::optional(out);
    return out;
}

lie_algebra* lie_algebra::compute_centralizer() { //TODO: fix, broken
    if (this->centralizer.has_value()) {
        return this->centralizer.value();
    }
    // Let L have basis e_i.
    lie_algebra* out = compute_centralizer_element(this->basis[0], get_sl(this->get_sl_size())); // Sets out=C_{sl(n)}(e_1)
    for (int i = 1; i < this->dim; i++) {
        // Updates C_{sl(n)}(e_1,...,e_{i+1}) -> C_{C_{sl(n)}(e_1,...,e_i)}(e_{i+1})
        lie_algebra* out_new = compute_centralizer_element(this->basis[i], out);
        delete(out);
        out = out_new;
    }
    this->centralizer = stdx::optional< lie_algebra* >(out); // Records the centralizer
    return out;
}

lie_algebra* lie_algebra::compute_normalizer() {
    if (this->normalizer.has_value()) {
        return this->normalizer.value();
    }
    // Let L have basis {e_i, i<=r}. Set M_0 = sl(n) and M_{i+1}=N(x,L,M_i), where N(x,L,M) is the elements y of M such that ad(y) x in L.
    // It is clear that N(L)=\bigcap_{j<= r} N(x,L,sl(n))=M_r

    // Sets out=N_{sl(n)}(e_1,L)=M_1
    std::vector< g::matrix > out = this->compute_normalizer_element(this->basis[0], get_sl(this->get_sl_size())->get_basis());
    for (int i = 1; i < this->dim; i++) {
        // Updates out -> M_{i+1}
        std::vector< g::matrix > out_new = this->compute_normalizer_element(this->basis[i], out);
        delete(&out);
        out = out_new;
    }
    lie_algebra* out_alg = new lie_algebra(out, true); 
    this->normalizer = stdx::optional< lie_algebra* >(out_alg); // Records the normalizer
    return out_alg;
}

std::vector< g::matrix > lie_algebra::compute_normalizer_element(g::matrix x, std::vector< g::matrix > M) {


    lie_algebra* sl_alg = lie_algebra::get_sl(this->get_sl_size());

    // Computes the matrix of ad(x) with respect to the standard basis of M.
    g::exvector ad_x_on_basis = {};
    for (g::matrix v : M) {
        g::exvector temp = lin_alg::sl_ize(lin_alg::bracket(x,v), sl_alg->get_sl_size());
        ad_x_on_basis.insert(ad_x_on_basis.end(), temp.begin(), temp.end());
    }
    g::matrix adx = lin_alg::matricize(ad_x_on_basis, sl_alg->get_dim(), M.size());
    // Let sl = span(alpha) where alpha is the extension of the basis of L to sl.
    std::vector< g::matrix > alpha = this->extend_basis(sl_alg);
    // Let H be the span of the remaining basis elements of sl not in L. Let P (proj) be the projection onto H, in the basis alpha.
    g::matrix proj = {sl_alg->get_dim(), sl_alg->get_dim()};
    for (int i = this->get_dim(); i < sl_alg->get_dim(); i++) {
        proj(i,i) = 1;
    }
    
    // Vectorize the matrices in alpha, put them as columns of a matrix and then invert to get the change of basis from std of sl to alpha
    g::exvector alpha_to_sl = {};
    for (g::matrix v : alpha) {
        g::exvector temp_insert = lin_alg::sl_ize(v, this->get_sl_size());
        alpha_to_sl.insert(alpha_to_sl.end(), temp_insert.begin(), temp_insert.end());
    }
    g::matrix alpha_to_sl_matrix = lin_alg::matricize(alpha_to_sl, sl_alg->get_dim(), sl_alg->get_dim());
    g::matrix sl_to_alpha = alpha_to_sl_matrix.inverse();
 
    // TODO: to save time we can precompute Id_{std_sl}^{alpha}, alpha, in compute_normalizer.
    // N(x,L,M) = ker(P * ad(x)), which we can compute as null(P_{alpha}^{alpha} * Id_{std_sl}^{alpha} * ad(x)_{std_M}^{std_sl}),
    // since y is in ker(P) iff y is in L
    g::matrix padx = {sl_alg->get_dim(), M.size()};
    padx = proj.mul(sl_to_alpha).mul(adx);
    std::vector< g::exvector > null_basis_M_alpha = lin_alg::nullspace(padx);
    // Puts the nullspace in the standard basis of sl. // TODO: Don't need helper function
    std::vector< g::matrix > null_basis_sl = {};
    std::vector< g::matrix > sl_basis = sl_alg->get_basis();
    for (g::exvector v : null_basis_M_alpha) {
        g::matrix v_sl_basis = alpha_to_sl_matrix.mul(lin_alg::matricize(v,sl_alg->get_dim(), 1));        null_basis_sl.push_back(lin_alg::vector_to_matrix(v_sl_basis, sl_basis));
    }
    return null_basis_sl;
}

lie_algebra* lie_algebra::bracket_with_sl() { // Needs to be changed if we decide to go the other way about making sl a static member
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
    if (this->derived_series.has_value()) {
        return std::vector(this->derived_series.value());
    }
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
    if (this->lower_central_series.has_value()) {
        return std::vector(this->lower_central_series.value());
    }
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

bool lie_algebra::is_solvable() {
    if (derived_series.has_value()) { //If we computed the derived series we just check if the last term is 0
        lie_algebra* ds = derived_series.value()[derived_series.value().size() - 1];
        return ds->dim == 0;
    }
    // Otherwise we use Cartan's criterion and check if tr(L,[L,L]), and it suffices
    //  to check this on a basis of L and of [L,L]
    lie_algebra* da = compute_derived_subalgebra();
    for (g::matrix a : basis) {
        for (g::matrix b : da->basis) {
            if (lin_alg::prod_trace(a,b).is_zero()) {
                return false;
            }
        }
    }
    return true;
}

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

lie_algebra* compute_centralizer_element(g::matrix x, lie_algebra *N) {
    std::vector< g::matrix > basis = std::vector< g::matrix >();
    std::vector< g::matrix > old_basis = N->get_basis();
    for (int i = 0; i < old_basis.size(); i++) { // Forms the matrix of ad(x) acting on N
        old_basis[i] = lin_alg::bracket(x, old_basis[i]);
    }
    g::matrix m = lin_alg::basis_to_vectorized_matrix(old_basis);
    std::vector< g::exvector > vec_null_basis = lin_alg::nullspace(m);
    for (g::exvector v : vec_null_basis) {
        basis.push_back(lin_alg::matricize(v, N->get_sl_size(), N->get_sl_size())); //TODO: why you like this ginac
    }
    lie_algebra* l = new lie_algebra(basis, true);
    return l;
}


lie_algebra* bracket_lie_algebras(lie_algebra *algebra1, lie_algebra  *algebra2) {
    std::vector< g::matrix > basis = std::vector< g::matrix >();
    for (g::matrix v1 : algebra1->get_basis()) {
        for (g::matrix v2 : algebra2->get_basis()) {
            basis.push_back(lin_alg::bracket(v1, v2));
        }
    }
    basis = lin_alg::spanning_subsequence(basis);
    lie_algebra* out = new lie_algebra(basis, true);
    return out;
}