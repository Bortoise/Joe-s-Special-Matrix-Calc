//
// Created by cole on 24/05/23.
//
#include "headers/utils.h"
#include "headers/lin_alg.h"
#include "headers/lie_algebra.h"
#include <random>

/*
void test_spanning_subsequence(){
    g::symbol x("x");
    g::matrix a = {{0,0},{0,1}};
    g::matrix b = {{0,1},{0,0}};
    g::matrix c = {{0,0},{1,0}};
    //g::matrix d = {{0,1}, {1,1}};
    //g::matrix d(2,2);
    std::vector< g::matrix > v;
    v.push_back(a);
    v.push_back(b);
    //v.push_back(d);
    v.push_back(c);

        
    std::vector <g::matrix> v2 = lin_alg::spanning_subsequence(v);

    for(g::matrix i : v2){
        utils::print_matrix(i);
        std::cout<<std::endl;
    }
}
*/
void test_lie_alg_equals() {
    g::matrix a = {{1,0},{0,-1}};
    g::matrix b = {{0,1},{0,0}};
    g::matrix c = {{0,0},{1,0}};
    lie_algebra alg_1 = {{a,b,c}};
    lie_algebra alg_2 = {{b,c}};
    if(! alg_1.equals(&alg_2) ) {
        throw std::invalid_argument("I goofed");
    }
    lie_algebra alg_3 = {{a}};
    if ( alg_1.equals(&alg_3)) {
        throw std::invalid_argument("I goofed");
    }
}

void test_get_sl(int n) {
    lie_algebra* sl = lie_algebra::get_sl(n);
    if (sl->get_dim() != n*n-1) {
        throw std::invalid_argument("I goofed");
    }
}

void test_bracket_algebra_sl_sl(int n) {
    lie_algebra* sl = lie_algebra::get_sl(n);
    if (! sl->equals(bracket_lie_algebras(sl, sl))) {
        throw std::invalid_argument("I goofed");
    }
}

lie_algebra* get_L5_1() {
    g::matrix a = {{1,4,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,-3}};
    g::matrix b = {{0,1,0,0},{0,0,1,0},{0,0,0,0},{0,0,0,0}};
    g::matrix c = {{0,0,1,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    std::vector< g::matrix> matrices = {a,b,c};
    lie_algebra* out = new lie_algebra(matrices);
    return out;
}

lie_algebra* get_L5_2() {
    g::matrix a = {{1,0,0,0},{0,1,0,4},{0,0,-3,0},{0,0,0,1}};
    g::matrix b = {{0,1,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    g::matrix c = {{0,0,1,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    std::vector< g::matrix> matrices = {a,b,c};
    lie_algebra* out = new lie_algebra(matrices);
    return out;
}
/*
void benchmark_spanning_subsequence(){

    std::srand(1984);

    std::vector<g::matrix > v;
    //g::symbol x("x");
    


    for(int i = 0; i < 100; i++){
        g::matrix m(4,4);
        for(int i = 0; i < 4; i++){
            for(int j = 0; j < 4; j++){
                m.set(i,j,std::rand()%10);

            }
        }
        v.push_back(m);
    }
    //v[0].set(2,3,x);

    v = lin_alg::spanning_subsequence(v);
    std::cout<<v.size();



}
*/

int main() {
//    test_lie_alg_equals();
//    test_get_sl(6);
//    test_bracket_algebra_sl_sl(6);
    lie_algebra* alg = get_L5_1();
//    alg->compute_derived_series();
//    alg->compute_lower_central_series();
//    std::vector< lie_algebra* > ds = alg->compute_derived_series();
//    std::vector< lie_algebra* > lcs = alg->compute_lower_central_series();
//    std::cout << ds.size() << std::endl;
//    for (lie_algebra* a : ds) {
//        std::cout << a->get_dim() << std::endl;
////        utils::print_matrices(a->get_basis());
//    }
//    std::cout << std::endl;
//    std::cout << lcs.size() << std::endl;
//    for (lie_algebra* a : lcs) {
//        std::cout << a->get_dim() << std::endl;
////        utils::print_matrices(a->get_basis());
//    }
    std::cout << alg->compute_normalizer()->get_dim() << std::endl;
    utils::print_matrices(alg->compute_normalizer()->get_basis());
//    for (int n = 1; n < 15; n++){ // This causes a segfault cause only 1 static sl is being allowed at a time
//        test_get_sl(n);
//    }
    return 0;

//    g::symbol x("x");
//    g::symbol y("y");
//
//    g::matrix a = {{0,x,0,1},{0,0,y,1},{0,0,0,0}};
//    g::exvector v = lin_alg::vectorize(a);
//    g::matrix b = lin_alg::matricize(v, 3, 4);
//    utils::print_matrix(a);
//    std::cout<<std::endl;
//    utils::print_matrix(b);
//
//    /*
//    g::matrix a = {{0,x,0,1},{0,0,y,1},{0,0,0,0}};
//    std::cout << "Printing basis of nullspace of:" << std::endl;
//    utils::print_matrix(a);
//    std::cout << "Gaussian elimination:" << std::endl;
//    g::matrix b = lin_alg::gaussian_elimination(a);
//    utils::print_matrix(b);
//    std::cout << "Basis of nullspace:" << std::endl;
//    std::vector< g::exvector > v = lin_alg::nullspace(a);
//    utils::print_exvectors(v);
//    std::cout << "Dimension of nullspace: " << v.size() << std::endl;
//
//    g::matrix c = {{0,0,1,1},{0,0,x,0},{1,2*y,x,0}};
//    g::matrix d = lin_alg::gaussian_elimination(c);
//    std::cout << "Printing gaussian elimination of:" << std::endl;
//    utils::print_matrix(c);
//    std::cout << "Gaussian elimination:" << std::endl;
//    utils::print_matrix(d);
//    */
//
//    g::matrix e = {{0,0},{0,1}};
//    g::matrix f = {{lin_alg::prod_trace(e,e)}};
//    utils::print_matrix(f);
//    return 0;
}
