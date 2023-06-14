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

void test_get_normalizer_element() {
    g::matrix a = {{0,1,0,0},{0,0,1,0},{0,0,0,0},{0,0,0,0}};
    g::matrix b = {{1,0,0,0}, {0,1,0,0}, {0,0,1,0}, {0,0,0,-3}};
    std::vector<g::matrix> L_basis = {b};
    lie_algebra* sl = lie_algebra::get_sl(4);
    lie_algebra* L = new lie_algebra(L_basis);

    std::vector< g::matrix > M = sl->get_basis();
    std::vector< g::matrix > normalizer = L->compute_normalizer_element(a, M);
    std::cout << "Printing a" << std::endl;
    utils::print_matrix(a);
    std::cout << "Printing basis of L" << std::endl;
    utils::print_matrices(L->get_basis());
    std::cout << "Printing normalizer of a in M" << std::endl;
    utils::print_matrices(normalizer);

    std::vector<g::matrix> brackets;
    for (g::matrix v : normalizer){
        brackets.push_back(lin_alg::bracket(a,v));
    }
    std::cout << "Printing brackets" << std::endl;
    utils::print_matrices(brackets);

    std::cout << "Computing normalizer of b in M" << std::endl;
    normalizer = L->compute_normalizer_element(b, M);
    utils::print_matrices(normalizer);
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
    test_get_normalizer_element();

    // lie_algebra* alg = get_L5_1();

    // std::cout << alg->compute_normalizer()->get_dim() << std::endl;
    // utils::print_matrices(alg->compute_normalizer()->get_basis());

    return 0;
}
