//
// Created by cole on 24/05/23.
//
#include "headers/utils.h"
#include "headers/lin_alg.h"
#include "headers/lie_algebra.h"
#include <random>


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

    std::cout << "Printing a" << std::endl;
    utils::print_matrix(a);
    std::cout << "Printing gaussian(a)" << std::endl;
    utils::print_matrix(gaussian_a_computed);

    std::cout << "Printing b" << std::endl;
    utils::print_matrix(b);
    std::cout << "Printing gaussian(b)" << std::endl;
    utils::print_matrix(gaussian_b_computed);

    std::cout << "Printing c" << std::endl;
    utils::print_matrix(c);
    std::cout << "Printing gaussian(c)" << std::endl;
    utils::print_matrix(gaussian_c_computed);

    std::cout << "Printing d" << std::endl;
    utils::print_matrix(d);
    std::cout << "Printing gaussian(d)" << std::endl;
    utils::print_matrix(gaussian_d_computed);

    std::cout << "Printing e" << std::endl;
    utils::print_matrix(e);
    std::cout << "Printing gaussian(e)" << std::endl;
    utils::print_matrix(gaussian_e_computed);
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

    std::cout << "Printing a" << std::endl;
    utils::print_matrix(a);
    std::cout << "Printing null(a)" << std::endl;
    utils::print_column_matrices(null_a_computed);

    std::cout << "Printing b" << std::endl;
    utils::print_matrix(b);
    std::cout << "Printing null(b)" << std::endl;
    utils::print_column_matrices(null_b_computed);

    std::cout << "Printing c" << std::endl;
    utils::print_matrix(c);
    std::cout << "Printing null(c)" << std::endl;
    utils::print_column_matrices(null_c_computed);

    std::cout << "Printing d" << std::endl;
    utils::print_matrix(d);
    std::cout << "Printing null(d)" << std::endl;
    utils::print_column_matrices(null_d_computed);

    std::cout << "Printing e" << std::endl;
    utils::print_matrix(e);
    std::cout << "Printing null(e)" << std::endl;
    utils::print_column_matrices(null_e_computed);

    std::cout << "Printing f" << std::endl;
    utils::print_matrix(f);
    std::cout << "Printing null(f)" << std::endl;
    utils::print_column_matrices(null_f_computed);

    std::vector< g::exvector > null_a_actual = {{1,-1,0},{1,0,-1}};    
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

    g::matrix a_hopefully = lin_alg::vector_to_matrix(slize_a, sl_basis);
    g::matrix b_hopefully = lin_alg::vector_to_matrix(slize_b, sl_basis);
    g::matrix c_hopefully = lin_alg::vector_to_matrix(slize_c, sl_basis);

    sl = lie_algebra::get_sl(4);
    sl_basis = sl->get_basis();
    g::matrix d_hopefully = lin_alg::vector_to_matrix(slize_d, sl_basis);

    std::cout << "Printing a" << std::endl;
    utils::print_matrix(a);
    std::cout << "Printing sl_ize(a)" << std::endl;
    utils::print_exvector(slize_a);
    std::cout << "Printing a?" << std::endl;
    utils::print_matrix(a_hopefully);

    std::cout << "Printing b" << std::endl;
    utils::print_matrix(b);
    std::cout << "Printing sl_ize(b)" << std::endl;
    utils::print_exvector(slize_b);
    std::cout << "Printing b?" << std::endl;
    utils::print_matrix(b_hopefully);

    std::cout << "Printing c" << std::endl;
    utils::print_matrix(c);
    std::cout << "Printing sl_ize(c)" << std::endl;
    utils::print_exvector(slize_c);      
    std::cout << "Printing c?" << std::endl;
    utils::print_matrix(c_hopefully);  

    std::cout << "Printing d" << std::endl;
    utils::print_matrix(d);
    std::cout << "Printing sl_ize(d)" << std::endl;
    utils::print_exvector(slize_d);      
    std::cout << "Printing d?" << std::endl;
    utils::print_matrix(d_hopefully);  
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
    
    std::cout << "Printing a" << std::endl;
    utils::print_matrix(a);
    std::cout << "Printing matricize(a,1,5)" << std::endl;
    utils::print_matrix(a1_matricize);
    std::cout << "Printing matricize(a,5,1)" << std::endl;
    utils::print_matrix(a2_matricize);

    std::cout << "Printing b" << std::endl;
    utils::print_matrix(b);
    std::cout << "Printing matricize(b,2,2)" << std::endl;
    utils::print_matrix(b1_matricize);
    std::cout << "Printing matricize(b,2,2,true)" << std::endl;
    utils::print_matrix(b2_matricize);

    std::cout << "Printing c" << std::endl;
    utils::print_exvector(c);
    std::cout << "Printing matricize(c,3,3)" << std::endl;
    utils::print_matrix(c1_matricize);
    std::cout << "Printing matricize(c,3,3,true)" << std::endl;
    utils::print_matrix(c2_matricize);

    std::cout << "Printing d" << std::endl;
    utils::print_exvector(d);
    std::cout << "Printing matricize(d,2,3)" << std::endl;
    utils::print_matrix(d1_matricize);
    std::cout << "Printing matricize(d,2,3,true)" << std::endl;
    utils::print_matrix(d2_matricize);
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
    
    
    // mat_vec out_basis = out->get_basis();
    // utils::print_matrices(out_basis);
    // g::matrix t = {{-3,0,0,0},{0,-1,0,0},{0,0,1,0},{0,0,0,3}};
    // g::matrix ta = lin_alg::bracket(t,a);
    // g::matrix tb = lin_alg::bracket(t,b);
    // g::matrix tc = lin_alg::bracket(t,c);
    // std::cout << "Printing brackets:" << std::endl;
    // std::cout << std::endl;
    // utils::print_matrix(ta);
    // std::cout << std::endl;
    // std::cout << std::endl;
    // utils::print_matrix(tb);
    // std::cout << std::endl;
    // std::cout << std::endl;
    // utils::print_matrix(tc);

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

lie_algebra* get_test_alg_1() {
    g::matrix a = {{0,1,0,0},{0,0,1,0},{0,0,0,1},{0,0,0,0}};
    g::matrix b = {{0,0,176,0},{0,0,0,176+1},{0,0,0,0},{0,0,0,0}};
    g::matrix c = {{0,0,0,1},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    std::vector< g::matrix> matrices = {a,b,c};
    lie_algebra* out = new lie_algebra(matrices);

    return out;
}

void test_normalizer() {
    g::symbol x("x");
    lie_algebra* L5_1 = get_L5_1();
    lie_algebra* L5_2 = get_L5_2();
    lie_algebra* L5_4 = get_L5_4(x);
    lie_algebra* L5_5 = get_L5_5(x);
    lie_algebra* test_1 = get_test_alg_1();

    lie_algebra* normalizer_L5_1 = L5_1->compute_normalizer();
    lie_algebra* normalizer_L5_2 = L5_2->compute_normalizer();
    lie_algebra* normalizer_L5_4 = L5_4->compute_normalizer();
    lie_algebra* normalizer_L5_5 = L5_5->compute_normalizer();
    lie_algebra* normalizer_test_1 = test_1->compute_normalizer();
    
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
        {{0, 0, 1, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}, //x
        {{0, (2*x + 1)/(x+1),  0,  0}, {0,  0,  1,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}},//x
        {{-3,  0,  0,  0}, {0, -1,  0,  0}, {0,  0,  1,  0}, {0,  0,  0,  3}},
        {{0,  0,  0,  1}, {0,  0,  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}},//x
        {{0,  0,  0,  0}, {0,  0,  0,  1}, {0,  0,  0,  0}, {0,  0,  0,  0}},//x
        {{0,  -x/(x+1), 0,  0}, {0,  0,  0,  0}, {0,  0,  0,  1}, {0,  0,  0, 0}}//x
    };

    mat_vec test_1_normalizer_basis = {
        {{0, 0, 1, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}}, //x
        {{0, (2*176 + 1)/(176+1),  0,  0}, {0,  0,  1,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}},//x
        {{-3,  0,  0,  0}, {0, -1,  0,  0}, {0,  0,  1,  0}, {0,  0,  0,  3}},
        {{0,  0,  0,  1}, {0,  0,  0,  0}, {0,  0,  0,  0}, {0,  0,  0,  0}},//x
        {{0,  0,  0,  0}, {0,  0,  0,  1}, {0,  0,  0,  0}, {0,  0,  0,  0}},//x
        {{0,  -176/(176+1), 0,  0}, {0,  0,  0,  0}, {0,  0,  0,  1}, {0,  0,  0, 0}}//x
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
    lie_algebra normalizer_actual_test_1 = {test_1_normalizer_basis, true};

    /*
    mat_vec L5_4_basis = L5_4->get_basis();
    mat_vec out = lie_algebra::get_sl(L5_4->get_sl_size())->get_basis();
    for (int i = 0; i < L5_4->get_dim(); i++) {              
        g::matrix adx = L5_4->compute_normalizer_element1(L5_4_basis[i], out);
        mat_vec extension = L5_4->compute_normalizer_element2(L5_4_basis[i], out);
        g::matrix alpga2asl = L5_4->compute_normalizer_element3(L5_4_basis[i], out);
        g::matrix alpga2aslminusone = L5_4->compute_normalizer_element4(L5_4_basis[i], out);
        g::matrix padx = L5_4->compute_normalizer_element5(L5_4_basis[i], out);
        mat_vec nullspasis = L5_4->compute_normalizer_element6(L5_4_basis[i], out);

        std::cout << "adx" << std::endl;
        utils::print_matrix(adx);
        std::cout << "extension" << std::endl;
        utils::print_matrices(extension);
        std::cout << "alpga2asl" << std::endl;
        utils::print_matrix(alpga2asl);
        std::cout << "alpga2aslminusone" << std::endl;
        utils::print_matrix(alpga2aslminusone);
        std::cout << "padx" << std::endl;
        utils::print_matrix(padx);
        std::cout << "nullspasis" << std::endl;
        utils::print_matrices(nullspasis);


        // Updates out -> M_{i+1}  
        out = L5_4->compute_normalizer_element(L5_4_basis[i], out);
        
        std::cout <<"nullspace as matrices" << std::endl;
        utils::print_matrices(out);
    }
    */

    mat_vec normalizer_test_basis = normalizer_test_1->get_basis();
    utils::print_matrices(normalizer_test_basis);
    
    if(! normalizer_actual_test_1.equals(normalizer_test_1)) {
        throw std::invalid_argument("I normalizer goofed with test 1");
    }
    if(! normalizer_actual_L5_1.equals(normalizer_L5_1)) {
        throw std::invalid_argument("I normalizer goofed with L5_1");
    }
    if(! normalizer_actual_L5_2.equals(normalizer_L5_2)) {
        throw std::invalid_argument("I normalizer goofed with L5_2");
    }
    if(! normalizer_actual_L5_4.equals(normalizer_L5_4)) {
        // mat_vec normalizer_basis = normalizer_L5_4->get_basis();
        // utils::print_matrices(normalizer_basis);
        // std::cout << normalizer_basis.size() << std::endl;
        throw std::invalid_argument("I normalizer goofed with L5_4");
    }
    if(! normalizer_actual_L5_5.equals(normalizer_L5_5)) {
        throw std::invalid_argument("I normalizer goofed with L5_5");
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
        {{0,  0,  1,  0}, {0,  1,  0,  1}, {0,  0,  0,  0}, {0,  0,  0,  0}},
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
        mat_vec v = centralizer_L5_1->get_basis();
        utils::print_matrices(v);
        throw std::invalid_argument("I centralizer goofed with L5_1");
    }
    if(! centralizer_actual_L5_2.equals(centralizer_L5_2)) {
        throw std::invalid_argument("I centralizer goofed with L5_2");
    }

    if(! centralizer_actual_L5_4.equals(centralizer_L5_4)) {
        throw std::invalid_argument("I centralizer goofed with L5_4");
    }
    if(! centralizer_actual_L5_5.equals(centralizer_L5_5)) {
        throw std::invalid_argument("I centralizer goofed with L5_5");
    }
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
    // test_lie_alg_equals();
    // test_get_sl(6);
    // test_bracket_algebra_sl_sl(6);
    // test_get_normalizer_element();
    // test_gaussian_elimination();
    // test_nullspace();
    // test_spanning_subsequence();
    // test_sl_ize();
    // test_matricize();
    test_normalizer();
    // test_centralizer();

    // lie_algebra* alg = get_L5_1();

    // std::cout << alg->compute_normalizer()->get_dim() << std::endl;
    // utils::print_matrices(alg->compute_normalizer()->get_basis());

    return 0;
}
