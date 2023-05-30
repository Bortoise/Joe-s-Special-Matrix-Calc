//
// Created by cole on 24/05/23.
//
#include "headers/utils.h"
#include "headers/lin_alg.h"


   g::exvector vectorize(g::matrix &m){
        g::exvector v;
        int r,c;
        r = m.rows();
        c = m.cols();
        for(int i = 0; i < r; i++){
            for(int j = 0; j < r; j++){
                v.push_back(m[i][j]);
            }
        }
        return v;
    }


int main() {
    g::symbol x("x");
    g::symbol y("y");
    g::matrix a = {{1,0},{0,0}};
    g::matrix b = {{0,1},{0,0}};
    g::exvector v = lin_alg::vectorize(a);


    //g::matrix ab = bracket(a, b);
    // print_matrix(ab);
    for(auto i : v){
        i.print(g::print_latex(std::cout));
        std::cout<<std::endl;
    }


    //a.add(b).print(g::print_latex(std::cout));
    //a.print(g::print_latex(std::cout));
}
