//
// Created by cole on 24/05/23.
//
#include "headers/utils.h"
#include "headers/lin_alg.h"


std::vector < g::matrix  > hi(){
    g::matrix a = {{0,0},{0,1}};
    g::matrix b = {{0,1},{0,0}};
    g::matrix c = {{0,0},{1,0}};
    g::matrix d = {{0,1}, {1,1}};
    //g::matrix d(2,2);
    std::vector< g::matrix > v;
    v.push_back(a);
    v.push_back(b);
    v.push_back(d);
    v.push_back(c);

    return v;
}

int main() {


    //g::matrix ab = bracket(a, b);
    // print_matrix(ab);
    // for(int i = 0; i < v.size(); i++){
    //     v[i].print(g::print_python(std::cout));
    //     std::cout<<std::endl;
    // }

    g::matrix c = 
    {{0, 1, 1, 0},
     {1, 1, 0, 0},
     {0, 0, 0, 0},
     {0, 0, 0, 1}};
    // std::cout << c.determinant() << std::endl;
    std::vector < g::matrix > v = hi();
    //(v[3]).print(g::print_python(std::cout));
    //utils::print_matrix(v[3]);


    
    std::vector <g::matrix> v2 = lin_alg::spanning_subsequence(v);

    for(g::matrix i : v2){
        utils::print_matrix(i);
        std::cout<<std::endl;
    }





    
        
    //utils::print_matrix(c);
    //utils::print_matrix(m);
    //lin_alg::swap_rows(&m, 0, 1);
    //utils::print_matrix(m);


    //a.add(b).print(g::print_latex(std::cout));
    //a.print(g::print_latex(std::cout));
    return 0;
}
