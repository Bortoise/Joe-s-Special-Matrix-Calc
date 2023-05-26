//
// Created by cole on 24/05/23.
//
#include "headers/utils.h"


int main() {
    g::symbol x("x");
    g::symbol y("y");
    g::matrix a = {{1,0},{0,0}};
    g::matrix b = {{0,1},{0,0}};
    //g::matrix ab = bracket(a, b);
    // print_matrix(ab);
    a.print(g::print_latex(std::cout));
}
