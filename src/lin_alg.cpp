#include "headers/lin_alg.h"
#include "headers/utils.h"


namespace lin_alg{

    
    g::exvector vectorize(g::matrix &m){
        g::exvector v;
        int r,c;
        r = m.rows();
        c = m.cols();
        for(int i = 0; i < m.nops(); i++){
            v.push_back(m.op(i));
        }
        return v;
    }


    g::matrix bracket(g::matrix &x, g::matrix &y){
        // Compute [x,y] = xy - yx
        return x.mul(y).sub(y.mul(x)); 
    }

    std::vector< g::matrix* > spanning_subsequence(std::vector< g::matrix* > &matrices){
        // Remove linear dependences from a set of matrices
        std::vector< g::exvector > vectorized;
        unsigned int tr, tc;
        tr = vectorized[0].size();
        tc = vectorized.size();
        g::matrix temp(tr, tc);
        
        for(unsigned int i = 0; i < tr; i++){
            for(unsigned int j = 0; j < tr; j++){
                temp.set(i,j, vectorized[i][j]);
            }
        }
        temp = gaussian_elimination_columnless(temp);
        //find remaining rows

        //filter out remaining rows
    

    }

    void swap_rows(g::matrix *mat, int row1, int row2){
        // Swap two rows in a matrix
        for (int i = 0; i < mat->cols(); i++){
            g::ex temp = (*mat)(row1,i);
            (*mat)(row1,i) = (*mat)(row2,i);
            (*mat)(row2,i) = temp;
        }
    }   



    g::matrix gaussian_elimination_columnless(g::matrix const &mat){
        // Perform gaussian elimination on a matrix without swapping any columns

        // Make a deep copy of the matrix to perform operations on
        g::matrix mat_temp = g::matrix(mat.rows(), mat.cols());    
        g::ex pivot;
        int pivot_row;

        // Copies the matrix
        for (int i = 0; i < mat.rows(); i++){
            for (int j = 0; j < mat.cols(); j++){
                mat_temp(i,j) = mat(i,j);
            }
        }

        for (int i = 0; i < mat_temp.cols(); i++){
            // Find a pivot in the ith column

            pivot = 0;
            pivot_row = -1;
            for (int j = i; j < mat_temp.rows(); j++){
                // Check if the pivot is non-zero
                if (!mat_temp(j,i).is_zero()){
                    mat_temp(j,i).print(g::print_latex(std::cout));
                    std::cout << std::endl;
                    pivot = mat_temp(j,i);
                    pivot_row = j;
                    break;
                }
            }

            if (pivot_row != -1){
                swap_rows(&mat_temp, i, pivot_row);
                for (int j = i+1; j < mat_temp.rows(); j++){
                    // Subtract the ith row from the jth row
                    std::cout << "Before:" << std::endl;
                    utils::print_matrix(mat_temp);
                    for (int k = i; k < mat_temp.cols(); k++){ 
                        std::cout << "(" << i << "," << j << "," << k << ")" << std::endl;
                        std::cout << mat_temp(j,k) << mat_temp(j,k) << mat_temp(i,i) << mat_temp(i,k) << mat_temp(j,i) << std::endl;
                        mat_temp.set(j,k,(mat_temp(j,k) * mat_temp(i,i)) - mat_temp(i,k) * mat_temp(j,i));
                        std::cout << mat_temp(j,k) << std::endl;
                    }
                    std::cout << "After:" << std::endl;
                    utils::print_matrix(mat_temp);
                }
            }
        }

        return mat_temp;
    }
};