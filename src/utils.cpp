//
// Created by cole on 15/05/23.
//

#include "headers/utils.h"


g::matrix bracket(g::matrix &m, g::matrix &n){
    return m.mul(n).sub(n.mul(m));
}
 	
// static void my_print(const g::ex & e)
// {   
//     if (g::is_a<g::function>(e))
//         std::cout << g::ex_to<g::function>(e).get_name();
//     else
//         std::cout << g::ex_to<g::basic>(e).class_name();
//     std::cout << "(";
//     size_t n = e.nops();
//     if (n)
//         for (size_t i=0; i<n; i++) {
//             my_print(e.op(i));
//             if (i != n-1)
//                 std::cout << ",";
//         }
//     else
//         std::cout << e;
//     std::cout << ")";
// }

void print_matrix(g::matrix &m){
    for (int i = 0; i < m.rows(); i++){
        for (int j = 0; j < m.cols(); j++){
            m(i,j).print(g::print_latex(std::cout));
            std::cout<<" ";
        }
        std::cout << std::endl;
    }
}

g::ex joe() {
    g::matrix m = {{-1,1},{2,3}};
    return m(0,0);
}