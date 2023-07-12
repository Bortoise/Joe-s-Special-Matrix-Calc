//
// Created by cole on 15/05/23.
//

#include "headers/utils.h"

namespace utils{
    void print_matrix(g::matrix &m){
        for (int i = 0; i < m.rows(); i++){
            for (int j = 0; j < m.cols(); j++){
                m(i,j).print(g::print_python(std::cout));
                std::cout<<"       ";
            }
            std::cout << std::endl;
        }
    }

    void print_exvectors(std::vector< g::exvector > &v){
        for (g::exvector i : v){
            std::cout<<"(";
            for (std::vector<g::ex>::iterator j = i.begin(); j != i.end(); ++j){
                j->print(g::print_python(std::cout));
                if (j != i.end()-1) std::cout<<", ";
            }
            std::cout<<")"<<std::endl;
        }
    }

    void print_exvector(g::exvector &i){
        std::cout<<"(";
        for (std::vector<g::ex>::iterator j = i.begin(); j != i.end(); ++j){
            j->print(g::print_python(std::cout));
            if (j != i.end()-1) std::cout<<", ";
        }
        std::cout<<")"<<std::endl;
    }

    void print_column_matrices(mat_vec &v){
        for (g::matrix m : v){
            std::cout<<"(";
            for (int i = 0; i < m.rows(); i++){
                m(i,0).print(g::print_python(std::cout));
                if (i != m.rows()-1) std::cout<<", ";
            }
            std::cout<<")"<<std::endl;
        }
    }

    void print_matrices(mat_vec m){
        for (g::matrix i : m){
            print_matrix(i);
            std::cout << std::endl;
        }
    }

    bool matrix_eq(g::matrix &m, g::matrix &n) {
        return m.sub(n).is_zero_matrix();
    }

    
    bool exvector_eq(g::exvector &v, g::exvector &w) {
        if (v.size() != w.size()) {
            return false;
        }
        for(int i = 0; i < v.size(); i++) {
            if (!(v[i]-w[i]).is_zero()) {
                return false;
            }
        }
        return true;
    }
};