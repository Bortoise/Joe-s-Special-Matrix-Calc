//
// Created by cole on 15/05/23.
//

#include "LieAlgebra.h"
#include "iostream"
#include <ginac/matrix.h>


GiNaC::ex joe() {
    GiNaC::matrix m = {{-1,1},{2,0}};
    return m(0,0);
}

