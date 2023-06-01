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
        temp = gaussian_elimination(temp);
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



    g::matrix gaussian_elimination(g::matrix const &mat){
        // Perform gaussian elimination on a matrix without swapping any columns

        // Make a deep copy of the matrix to perform operations on
        g::matrix mat_temp = g::matrix(mat.rows(), mat.cols());    

        // Initialize the pivot and pivot row
        g::ex pivot;
        int pivot_row;

        // Initialize the number of non-zero columns, where to search for the pivot
        int numNonZeroCols = 0;

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
            for (int j = numNonZeroCols; j < mat_temp.rows(); j++){
                // Check if the pivot is non-zero
                if (!mat_temp(j,i).is_zero()){
                    pivot = mat_temp(j,i);
                    pivot_row = j;
                    break;
                }
            }

            if (pivot_row != -1){
                numNonZeroCols += 1; // Increment the number of non-zero columns to use as a starting point for the next column

                swap_rows(&mat_temp, i, pivot_row);
                for (int j = i+1; j < mat_temp.rows(); j++){
                    // Subtract the ith row from the jth row
                    g::ex temp_first_ex = mat_temp(j,i);
                    for (int k = i; k < mat_temp.cols(); k++){ 
                        mat_temp.set(j,k,(mat_temp(j,k) * mat_temp(i,i)) - mat_temp(i,k) * temp_first_ex);
                    }
                }
            }
        }

        return mat_temp;
    }
};