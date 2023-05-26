//
// Created by cole on 24/05/23.
//
#include "iostream"
//#include "LieAlgebra.h"
//
//using namespace GiNaC;
#include <ginac/matrix.h>

GiNaC::ex joe() {
    GiNaC::matrix m = {{-1,1},{2,0}};
    return m(0,0);
}

int main() {
    GiNaC::ex e = joe();
    std::cout << e;
}
