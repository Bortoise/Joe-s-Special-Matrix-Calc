//
// Created by cole on 24/05/23.
//
#include "headers/utils.h"
#include "headers/lin_alg.h"




int main() {
    g::symbol x("x");
    g::symbol y("y");
    g::matrix a = {{1,0},{0,0}};
    g::matrix b = {{0,1},{0,0}};
    g::exvector v = lin_alg::vectorize(a);


    //g::matrix ab = bracket(a, b);
    // print_matrix(ab);
    // for(int i = 0; i < v.size(); i++){
    //     v[i].print(g::print_python(std::cout));
    //     std::cout<<std::endl;
    // }

    g::matrix c = {{0,1, 3},{1,1, -2}, {-1, -1, 0}};
    // std::cout << c.determinant() << std::endl;
    g::matrix m = lin_alg::gaussian_elimination(c);
    
        
    //utils::print_matrix(c);
    utils::print_matrix(m);
    //lin_alg::swap_rows(&m, 0, 1);
    //utils::print_matrix(m);


    //a.add(b).print(g::print_latex(std::cout));
    //a.print(g::print_latex(std::cout));
    return 0;
}
