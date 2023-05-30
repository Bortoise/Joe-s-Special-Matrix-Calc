#include "headers/lin_alg.h"


namespace lin_alg{


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


    g::matrix bracket(g::matrix &x, g::matrix &y){
        // Compute [x,y] = xy - yx
        return x.mul(y).sub(y.mul(x)); 
    }

    std::vector< g::matrix* > spanning_subsequence(std::vector< g::matrix* > &matrices){
        // Remove linear dependences from a set of matrices
        std::vector< g::exvector* > vectorized;



    }

};

