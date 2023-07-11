//
// Created by cole on 24/05/23.
//
#include "headers/utils.h"
#include "headers/lin_alg.h"
#include "headers/lie_algebra.h"
#include <random>

//TODO: Add tests for derived / central series, max_rank

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
    g::matrix c = {{0,0,0,1},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    std::vector< g::matrix> matrices = {a,b,c};
    lie_algebra* out = new lie_algebra(matrices);
    return out;
}

lie_algebra* get_L5_4(g::symbol x) {
    g::matrix a = {{0,1,0,0},{0,0,1,0},{0,0,0,1},{0,0,0,0}};
    g::matrix b = {{0,0,x,0},{0,0,0,x+1},{0,0,0,0},{0,0,0,0}};
    g::matrix c = {{0,0,0,1},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    std::vector< g::matrix> matrices = {a,b,c};
    lie_algebra* out = new lie_algebra(matrices);

    return out;
}

lie_algebra* get_L5_5(g::symbol x) {
    g::matrix a = {{0,1,0,0},{0,0,0,0},{0,0,0,1},{0,0,0,0}};
    g::matrix b = {{0,0,x,0},{0,0,0,x+1},{0,0,0,0},{0,0,0,0}};
    g::matrix c = {{0,0,0,1},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    std::vector< g::matrix> matrices = {a,b,c};
    lie_algebra* out = new lie_algebra(matrices);
    return out;
}

lie_algebra* get_L3_025_5() {
    g::matrix m1 = {{-1, 0, 0, 0, 0, 2, 0, 0, 0, 0},
                    {0,  1, 0,  0, 0, 0, 0, 0, 0, 0},
                    {0,  0, -1, 0, 0, 0, 0, 0, 0, 0},
                    {0,  0, 0,  1, 0, 0, 0, 0, 0, 0},
                    {0,  0, 0,  0, 0, 0, 0, 1, 0, 0},
                    {0,  0, 1,  0, 0, -1, 0, 0, 0, 0},
                    {0,  0, 0,  0, 0, 0, 0, 0, 0, 1},
                    {0,  0, 0,  0, 0, 0, 0, 0, 0, 0},
                    {0,  0, 0,  0, 0, 0, 0, 0, 1, 0},
                    {0,  0, 0,  0, 0, 0, 0, 0, 0, 0}};
    g::matrix m2 = {{0, 0, 0, 0, 2, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 2, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 1, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 1, 1, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}};
    g::matrix m3 = {{0, 0, 0, 0, 0, 0, 1, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 1},
                    {0, 0, 0, 1, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
    mat_vec basis = mat_vec {m1,m2,m3};
    mat_vec brackets = {lin_alg::bracket(m1,m2), lin_alg::bracket(m2,m3), lin_alg::bracket(m1,m3)};
    g::matrix M = lin_alg::bracket(m1,m2).add(m2).sub(m3);
    utils::print_matrix(M);
    //utils::print_matrices(brackets);
    lie_algebra* L3 = new lie_algebra(basis);
    return L3;
}

lie_algebra* get_L3_0_9a(g::symbol x){
    g::matrix m1 = {{0, 0, 0, 0, 0, 0, 2, 0, 0, 0}, {0, 2, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, -2, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, -1, 0, 0, 0, 1}, {0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, -1}};
    g::matrix m2 = {{2*x, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 2*x+2, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, -6*x-2, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 2*x, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 2*x+1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, -2*x-1, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 2*x, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, -2*x, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 2*x+1, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, -2*x-1}};
    g::matrix m3 = {{0, 0, 0, 0, 2, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0, 0, 0,}, {0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 1, 0},  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
    mat_vec basis = {m1,m2,m3};
    lie_algebra* L3_0_9a = new lie_algebra(basis);
    return L3_0_9a;
}

void test_spanning_subsequence(){
    g::symbol x("x");
    g::symbol y("y");

    g::matrix a = {{1,0}};
    g::matrix b = {{2,0}};
    g::matrix c = {{3,1}};
    mat_vec v = {a,b,c};

    mat_vec subsequence = lin_alg::spanning_subsequence(v);

    std::cout << "Printing matrices" << std::endl;
    utils::print_matrices(v);
    std::cout << "Printing spanning subsequence" << std::endl;
    utils::print_matrices(subsequence);

    a = {{1,0},{0,1}};
    b = {{2,0},{2,1}};
    c = {{3,1},{1,2}};
    g::matrix d = {{3,1},{3,1}};
    v = {a,b,c,d};

    subsequence = lin_alg::spanning_subsequence(v);

    std::cout << "Printing matrices" << std::endl;
    utils::print_matrices(v);
    std::cout << "Printing spanning subsequence" << std::endl;
    utils::print_matrices(subsequence);
}

void test_lie_alg_equals() {
    g::matrix a = {{1,0},{0,-1}};
    g::matrix b = {{0,1},{0,0}};
    g::matrix c = {{0,0},{1,0}};
    lie_algebra alg_1 = {{a,b,c}};
    lie_algebra alg_2 = {{b,c}};
    if(! alg_1.equals(&alg_2) ) {
        throw std::runtime_error("I goofed");
    }
    lie_algebra alg_3 = {{a}};
    if ( alg_1.equals(&alg_3)) {
        throw std::runtime_error("I goofed");
    }
}

void test_get_sl(int n) {
    lie_algebra* sl = lie_algebra::get_sl(n);
    if (sl->get_dim() != n*n-1) {
        throw std::runtime_error("I goofed");
    }
}

void test_bracket_algebra_sl_sl(int n) {
    lie_algebra* sl = lie_algebra::get_sl(n);
    if (! sl->equals(bracket_lie_algebras(sl, sl))) {
        throw std::runtime_error("I goofed");
    }
}

void test_get_normalizer_element() {
    g::matrix a = {{0,1,0}, {0,0,1}, {0,0,0}};
    g::matrix b = {{1,0,0}, {0,1,0}, {0,0,-2}};
    g::matrix c = {{1,0,0}, {0,1,0}, {0,0,1}};
    g::matrix d = {3,3};
    d(1,2) = 1;

    mat_vec L_basis = {d};
    lie_algebra* sl = lie_algebra::get_sl(3);
    lie_algebra* L = new lie_algebra(L_basis);

    mat_vec M = sl->get_basis();
    mat_vec normalizer = L->compute_normalizer_element(d, M);

    // std::cout << "Printing d" << std::endl;
    // utils::print_matrix(d);
    // std::cout << "Printing basis of L" << std::endl;
    // utils::print_matrices(L->get_basis());
    // std::cout << "Printing dim of normalizer of d in M" << std::endl;
    // // utils::print_matrices(normalizer);
    // std::cout << normalizer.size() << std::endl;

    normalizer = L->compute_normalizer_element(a, M);
    std::cout << "Printing a" << std::endl;
    utils::print_matrix(a);
    std::cout << "Printing basis of L" << std::endl;
    utils::print_matrices(L->get_basis());
    std::cout << "Printing normalizer of a in M" << std::endl;
    utils::print_matrices(normalizer);

    mat_vec a_brackets;
    for (g::matrix v : normalizer){
        a_brackets.push_back(lin_alg::bracket(a,v));
    }
    std::cout << "Printing brackets" << std::endl;
    utils::print_matrices(a_brackets);

    std::cout << "Printing b" << std::endl;
    utils::print_matrix(b);
    std::cout << "Computing normalizer of b in M" << std::endl;
    normalizer = L->compute_normalizer_element(b, M);
    utils::print_matrices(normalizer);

    mat_vec b_brackets;
    for (g::matrix v : normalizer){
        b_brackets.push_back(lin_alg::bracket(b,v));
    }
    std::cout << "Printing brackets" << std::endl;
    utils::print_matrices(b_brackets);
}

void test_gaussian_elimination(){
    g::symbol x("x");
    g::symbol y("y");
    g::matrix a = {3,3};
    g::matrix b = {3,3};
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            b(i,j) = 1;
        }
    }
    g::matrix c = {{x, 0, 0}, {0, 0, y}};
    g::matrix d = {3,2};
    g::matrix e = {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,-2,-1},{0,0,0,1,0,0,0,0},{-1,0,0,0,1,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,-1,0,0,0,0},{0,1,0,0,0,-1,0,0},{0,0,0,0,0,1,0,0}};

    g::matrix gaussian_a_computed = lin_alg::gaussian_elimination(a);
    g::matrix gaussian_b_computed = lin_alg::gaussian_elimination(b);
    g::matrix gaussian_c_computed = lin_alg::gaussian_elimination(c);
    g::matrix gaussian_d_computed = lin_alg::gaussian_elimination(d);
    g::matrix gaussian_e_computed = lin_alg::gaussian_elimination(e);

    g::matrix gaussian_a_actual = {{0,0,0},{0,0,0},{0,0,0}};
    g::matrix gaussian_b_actual = {{1,1,1},{0,0,0},{0,0,0}};
    g::matrix gaussian_c_actual = {{x,0,0},{0,0,y}};
    g::matrix gaussian_d_actual = {{0,0},{0,0},{0,0}};
    g::matrix gaussian_e_actual = {{-1,0,0,0,1,0,0,0},{0,1,0,0,0,0,0,0},{0,0,0,1,0,0,0,0},{0,0,0,0,0,1,0,0},
                                        {0,0,0,0,0,0,-2,-1},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}};
                                        

    if(!utils::matrix_eq(gaussian_a_computed, gaussian_a_actual)) {
        throw std::runtime_error("gaussian_elimination fails on a");
    }
    if(!utils::matrix_eq(gaussian_b_computed, gaussian_b_actual)) {
        throw std::runtime_error("gaussian_elimination fails on b");
    }
    if(!utils::matrix_eq(gaussian_c_computed, gaussian_c_actual)) {
        throw std::runtime_error("gaussian_elimination fails on c");
    }
    if(!utils::matrix_eq(gaussian_d_computed, gaussian_d_actual)) {
        throw std::runtime_error("gaussian_elimination fails on d");
    }
    if(!utils::matrix_eq(gaussian_e_computed, gaussian_e_actual)) {
        throw std::runtime_error("gaussian_elimination fails on e");
    }
}

void test_nullspace(){
    g::symbol x("x");
    g::symbol y("y");
    g::matrix a = {3,3};
    g::matrix b = {3,3};
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            b(i,j) = 1;
        }
    }
    g::matrix c = {{x, 0, 0}, {0, 0, y}};
    g::matrix d = {3,2};
    g::matrix e = {{x, 1, 0}, {0, 0, y}};
    g::matrix f = {{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,-2,-1},{0,0,0,1,0,0,0,0},{-1,0,0,0,1,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,-1,0,0,0,0},{0,1,0,0,0,-1,0,0},{0,0,0,0,0,1,0,0}};

    mat_vec null_a_computed = lin_alg::nullspace(a);
    mat_vec null_b_computed = lin_alg::nullspace(b);
    mat_vec null_c_computed = lin_alg::nullspace(c);
    mat_vec null_d_computed = lin_alg::nullspace(d);
    mat_vec null_e_computed = lin_alg::nullspace(e);
    mat_vec null_f_computed = lin_alg::nullspace(f);

    std::vector< g::exvector > null_a_actual = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    std::vector< g::exvector > null_b_actual = {{1, -1, 0}, {1, 0, -1}};
    std::vector< g::exvector > null_c_actual = {{0, 1, 0}};
    std::vector< g::exvector > null_d_actual = {{1, 0}, {0, 1}};
    std::vector< g::exvector > null_e_actual = {{1, -x, 0}};
    std::vector< g::exvector > null_f_actual = {{0,0,1,0,0,0,0,0},{1,0,0,0,1,0,0,0},{0,0,0,0,0,0,-1,2}};


    if(!lin_alg::equals(&null_a_computed, &null_a_actual)) {
        throw std::runtime_error("nullspace fails on a");
    }
    if(!lin_alg::equals(&null_b_computed, &null_b_actual)) {
        throw std::runtime_error("nullspace fails on b");
    }
    if(!lin_alg::equals(&null_c_computed, &null_c_actual)) {
        throw std::runtime_error("nullspace fails on c");
    }
    if(!lin_alg::equals(&null_d_computed, &null_d_actual)) {
        throw std::runtime_error("nullspace fails on d");
    }
    if(!lin_alg::equals(&null_e_computed, &null_e_actual)) {
        throw std::runtime_error("nullspace fails on e");
    }
    if(!lin_alg::equals(&null_f_computed, &null_f_actual)) {
        throw std::runtime_error("nullspace fails on f");
    }
}

void test_sl_ize(){
    lie_algebra* sl = lie_algebra::get_sl(3);
    mat_vec sl_basis = sl->get_basis();

    g::matrix a = {{-14,2,3},{4,5,6},{7,8,9}};
    g::matrix b = {3,3};
    g::matrix c = {{-2,0,0},{0,1,0},{0,0,1}};
    g::matrix d = {{1, 3, 0, -6}, {2, 1, 0, 4}, {0, 0, 1, 8}, {-2, 4, 0, -3}};

    g::exvector slize_a = lin_alg::sl_ize(a,3);
    g::exvector slize_b = lin_alg::sl_ize(b,3);
    g::exvector slize_c = lin_alg::sl_ize(c,3);
    g::exvector slize_d = lin_alg::sl_ize(d,4);

    g::exvector slize_a_real = {2,4,3,7,6,8,-5,-9};
    g::exvector slize_b_real = {0,0,0,0,0,0,0,0};
    g::exvector slize_c_real = {0,0,0,0,0,0,-1,-1};
    g::exvector slize_d_real = {3,2,0,0,-6,-2,0,0,4,4,8,0,-1,-1,3};

    g::matrix a_reconstruction = lin_alg::vector_to_matrix(slize_a, sl_basis);
    g::matrix b_reconstruction = lin_alg::vector_to_matrix(slize_b, sl_basis);
    g::matrix c_reconstruction = lin_alg::vector_to_matrix(slize_c, sl_basis);

    sl = lie_algebra::get_sl(4);
    sl_basis = sl->get_basis();
    g::matrix d_reconstruction = lin_alg::vector_to_matrix(slize_d, sl_basis);


    // test slize correct
    if (!utils::exvector_eq(slize_a, slize_a_real)) { 
        throw std::runtime_error("slize fails on a");
    }
    if (!utils::exvector_eq(slize_b, slize_b_real)) {
        throw std::runtime_error("slize fails on b");
    }
    if (!utils::exvector_eq(slize_c, slize_c_real)) {
        throw std::runtime_error("slize fails on c");
    }
    if (!utils::exvector_eq(slize_d, slize_d_real)) {
        throw std::runtime_error("slize fails on d");
    }

    // test slize is inverse of vector_to_matrix with sl_basis
    if (!utils::matrix_eq(a, a_reconstruction)) {
        throw std::runtime_error("slize fails to be an inverse to vector_to_matrix on a");
    }
    if (!utils::matrix_eq(b, b_reconstruction)) {
        throw std::runtime_error("slize fails to be an inverse to vector_to_matrix on b");
    }
    if (!utils::matrix_eq(c, c_reconstruction)) {
        throw std::runtime_error("slize fails to be an inverse to vector_to_matrix on c");
    }
    if (!utils::matrix_eq(d, d_reconstruction)) {
        throw std::runtime_error("slize fails to be an inverse to vector_to_matrix on d");
    }
}

void test_matricize(){
    g::symbol x("x");
    g::matrix a = {{1},{2},{3},{4},{5}};
    g::matrix b = {{x}, {1}, {0}, {2}};

    g::exvector c = {1,2,3,4,5,6,7,8,9};
    g::exvector d = {x,2,0,-x,1,2};
    
    g::matrix a1_matricize = lin_alg::matricize(a,1,5);
    g::matrix a2_matricize = lin_alg::matricize(a,5,1);
    g::matrix b1_matricize = lin_alg::matricize(b,2,2);
    g::matrix b2_matricize = lin_alg::matricize(b,2,2,true);
    
    g::matrix c1_matricize = lin_alg::matricize(c,3,3);
    g::matrix c2_matricize = lin_alg::matricize(c,3,3,true);
    g::matrix d1_matricize = lin_alg::matricize(d,2,3);
    g::matrix d2_matricize = lin_alg::matricize(d,2,3,true);

    g::matrix a1_matricize_true = {{1,2,3,4,5}};
    g::matrix a2_matricize_true = {{1},{2},{3},{4},{5}};
    g::matrix b1_matricize_true = {{x,1},{0,2}};
    g::matrix b2_matricize_true = {{x,0},{1,2}};
    
    g::matrix c1_matricize_true = {{1,2,3},{4,5,6},{7,8,9}};
    g::matrix c2_matricize_true = {{1,4,7},{2,5,8},{3,6,9}};
    g::matrix d1_matricize_true = {{x,2,0},{-x,1,2}};
    g::matrix d2_matricize_true = {{x,0,1},{2,-x,2}};
    
    if(!utils::matrix_eq(a1_matricize, a1_matricize_true)) {
        throw std::runtime_error("matricize(a,1,5) failed");
    }
    if(!utils::matrix_eq(a2_matricize, a2_matricize_true)) {
        throw std::runtime_error("matricize(a,5,1) failed");
    }
    if(!utils::matrix_eq(b1_matricize, b1_matricize_true)) {
        throw std::runtime_error("matricize(b,2,2) by rows failed");
    }
    if(!utils::matrix_eq(b2_matricize,b2_matricize_true)) {
        throw std::runtime_error("matricize(b,2,2) by cols failed");
    }
    if(!utils::matrix_eq(c1_matricize,c1_matricize_true)) {
        throw std::runtime_error("matricize(c,3,3) by rows failed");
    }
    if(!utils::matrix_eq(c2_matricize,c2_matricize_true)) {
        throw std::runtime_error("matricize(c,3,3) by cols failed");
    }
    if(!utils::matrix_eq(d1_matricize,d1_matricize_true)) {
        throw std::runtime_error("matricize(d,2,3) by rows failed");
    }
    if(!utils::matrix_eq(d2_matricize,d2_matricize_true)) {
        throw std::runtime_error("matricize(d,2,3) by cols failed");
    }
}

void test_normalizer() {
    g::symbol x("x");
    lie_algebra* L5_1 = get_L5_1();
    lie_algebra* L5_2 = get_L5_2();
    lie_algebra* L5_4 = get_L5_4(x);
    lie_algebra* L5_5 = get_L5_5(x);
    lie_algebra* L3_0_9a = get_L3_0_9a(x);
    // lie_algebra* test_1 = get_test_alg_1();

    lie_algebra* normalizer_L5_1 = L5_1->compute_normalizer();
    lie_algebra* normalizer_L5_2 = L5_2->compute_normalizer();
    lie_algebra* normalizer_L5_4 = L5_4->compute_normalizer();
    lie_algebra* normalizer_L5_5 = L5_5->compute_normalizer();
    lie_algebra* normalizer_L3_0_9a = L3_0_9a->compute_normalizer();
    // lie_algebra* normalizer_test_1 = test_1->compute_normalizer();
    
    mat_vec L5_1_normalizer_basis = {
        {{0,  1,  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}}, 
        {{0,  0,  1,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}},
        {{0,  0,  0,  0}, {0,  0,  1,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}},
        {{1,  0,  0,  0}, {0,  1,  0,  0}, {0,  0,  1,  0}, {0,  0,  0,  -3}}
    };

    mat_vec L5_2_normalizer_basis = {
        {{1,  0,  0,  0}, {0,  0,  0,  0}, {0,  0,  -1,  0}, {0,  0,  0,  0}},
        {{0,  1,  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}},
        {{0,  0,  0,  1}, {0,  0,  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}},
        {{0,  0,  0,  0}, {0,  0,  0,  1}, {0,  0,  0,  0}, {0,  0,  0,  0}},
        {{0,  0,  0,  0}, {0,  1,  0,  0}, {0,  0,  -2,  0}, {0,  0,  0,  1}}
    };

    mat_vec L5_4_normalizer_basis = {
        {{0, 0, 1, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}},                               
        {{0, (2*x + 1)/(x+1),  0,  0}, {0,  0,  1,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}},      
        {{-3,  0,  0,  0}, {0, -1,  0,  0}, {0,  0,  1,  0}, {0,  0,  0,  3}},                  
        {{0,  0,  0,  1}, {0,  0,  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}},                   
        {{0,  0,  0,  0}, {0,  0,  0,  1}, {0,  0,  0,  0}, {0,  0,  0,  0}},                   
        {{0,  -x/(x+1), 0,  0}, {0,  0,  0,  0}, {0,  0,  0,  1}, {0,  0,  0, 0}}               
    };

    mat_vec L5_5_normalizer_basis = {
        {{0 , 1 , 0 , 0}, {0 , 0 , 0 , 0}, {0,  0,  0,  0}, {0,  0,  0,  0}},
        {{0 , 0 , 0 , 0}, {0 , 1 , 0 , 0}, {0,  0,  -1,  0}, {0,  0,  0,  0}},
        {{0 , 0 , 1 , 0}, {0 , 0 , 0 , 0}, {0,  0,  0,  0}, {0,  0,  0,  0}},
        {{0 , 0 , 0 , 1}, {0 , 0 , 0 , 0}, {0,  0,  0,  0}, {0,  0,  0,  0}},
        {{0 , 0 , 0 , 0}, {0 , 0 , 0 , 1}, {0,  0,  0,  0}, {0,  0,  0,  0}}, 
        {{0 , 0 , 0 , 0}, {0 , 0 , 0 , 0}, {0,  0,  0,  1}, {0,  0,  0,  0}},
        {{-1,  0,  0, 0},{ 0,  0,  0,  0}, { 0,  0,  0,  0}, {0,  0,  0,  1}}
    };

    lie_algebra normalizer_actual_L5_1 = {L5_1_normalizer_basis, true};
    lie_algebra normalizer_actual_L5_2 = {L5_2_normalizer_basis, true};
    lie_algebra normalizer_actual_L5_4 = {L5_4_normalizer_basis, true};
    lie_algebra normalizer_actual_L5_5 = {L5_5_normalizer_basis, true};

    if(!normalizer_actual_L5_1.equals(normalizer_L5_1)) {
        throw std::runtime_error("I normalizer goofed with L5_1");
    }
    if(!normalizer_actual_L5_2.equals(normalizer_L5_2)) {
        throw std::runtime_error("I normalizer goofed with L5_2");
    }
    if(!normalizer_actual_L5_4.equals(normalizer_L5_4)) {
        throw std::runtime_error("I normalizer goofed with L5_4");
    }
    if(!normalizer_actual_L5_5.equals(normalizer_L5_5)) {
        throw std::runtime_error("I normalizer goofed with L5_5");
    }
}

void test_centralizer() {
    g::symbol x("x");
    lie_algebra* L5_1 = get_L5_1();
    lie_algebra* L5_2 = get_L5_2();
    lie_algebra* L5_4 = get_L5_4(x);
    lie_algebra* L5_5 = get_L5_5(x);
    
    lie_algebra* sl = lie_algebra::get_sl(4);

    lie_algebra* centralizer_L5_1 = L5_1->compute_centralizer();
    lie_algebra* centralizer_L5_2 = L5_2->compute_centralizer();
    lie_algebra* centralizer_L5_4 = L5_4->compute_centralizer();
    lie_algebra* centralizer_L5_5 = L5_5->compute_centralizer();
    
    mat_vec L5_1_cent_basis = {
        {{1,  0,  0,  0}, {0,  1,  0,  0}, {0,  0,  1,  0}, {0,  0,  0,  -3}}, 
        {{0,  0,  1,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}}
    };

    mat_vec L5_2_cent_basis = {
        {{1,  0,  0,  0}, {0,  1,  0,  0}, {0,  0,  -3,  0}, {0,  0,  0,  1}},
        {{0,  0,  0,  1}, {0,  0,  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}}
    };
    
    mat_vec L5_4_cent_basis = {
        {{0,  0,  1,  0}, {0,  0,  0,  1}, {0,  0,  0,  0}, {0,  0,  0,  0}},
        {{0,  0,  0,  1}, {0,  0,  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}}
    };

    mat_vec L5_5_cent_basis = {
        {{0,  x/(x+1),  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  1}, {0,  0,  0,  0}},
        {{0,  0,  0,  1}, {0,  0,  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}},
        {{0,  0,  1,  0}, {0,  0,  0,  1}, {0,  0,  0,  0}, {0,  0,  0,  0}},
    };

    lie_algebra centralizer_actual_L5_1 = {L5_1_cent_basis, true};
    lie_algebra centralizer_actual_L5_2 = {L5_2_cent_basis, true};
    lie_algebra centralizer_actual_L5_4 = {L5_4_cent_basis, true};
    lie_algebra centralizer_actual_L5_5 = {L5_5_cent_basis, true};
    
    if(! centralizer_actual_L5_1.equals(centralizer_L5_1)) {
        throw std::runtime_error("I centralizer goofed with L5_1");
    }
    if(! centralizer_actual_L5_2.equals(centralizer_L5_2)) {
        throw std::runtime_error("I centralizer goofed with L5_2");
    }

    if(! centralizer_actual_L5_4.equals(centralizer_L5_4)) {
        throw std::runtime_error("I centralizer goofed with L5_4");
    }
    if(! centralizer_actual_L5_5.equals(centralizer_L5_5)) {
        throw std::runtime_error("I centralizer goofed with L5_5");
    }
}



int main() {
    // test_spanning_subsequence();
    // test_lie_alg_equals();
    // test_get_sl(6);
    // test_bracket_algebra_sl_sl(6);
    // test_get_normalizer_element();
    test_gaussian_elimination();
    test_nullspace();
    test_sl_ize();
    test_matricize();
    test_normalizer();
    test_centralizer();

    // lie_algebra* alg = get_L5_1();

    // std::cout << alg->compute_normalizer()->get_dim() << std::endl;
    // utils::print_matrices(alg->compute_normalizer()->get_basis());

    return 0;
}
