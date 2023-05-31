//
// Created by cole on 15/05/23.
//

#include "headers/utils.h"

namespace utils{
    void print_matrix(g::matrix &m){
        for (int i = 0; i < m.rows(); i++){
            for (int j = 0; j < m.cols(); j++){
                m(i,j).print(g::print_latex(std::cout));
                std::cout<<" ";
            }
            std::cout << std::endl;
        }
    }
};