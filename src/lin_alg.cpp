#include "headers/lin_alg.h"
#include "headers/utils.h"


namespace lin_alg{

    
    g::exvector vectorize(g::matrix &m){
        g::exvector v;
        for(int i = 0; i < m.nops(); i++){
            v.push_back(m.op(i));
        }

        return v;

    }

    g::matrix matricize(g::exvector &v, unsigned int r, unsigned int c){
        g::matrix x(r,c);
        for(int i  = 0; i < r; i++){
            for(int j = 0; j < c; j++){
                x.set(i,j, v[i*c + j]);
            }
        }
        return x;
    }



    g::matrix bracket(g::matrix &x, g::matrix &y){
        // Compute [x,y] = xy - yx
        return x.mul(y).sub(y.mul(x)); 
    }

    std::vector< g::matrix > spanning_subsequence(std::vector< g::matrix > &matrices){
        // Remove linear dependences from a set of matrices
        std::vector< g::exvector > vectorized;
        for(g::matrix i: matrices){
            vectorized.push_back(vectorize(i));
        }
        unsigned int tr, tc;
        tr = vectorized.size();
        tc = vectorized[0].size();
        g::matrix temp(tc, tr);
       
        for(unsigned int i = 0; i < vectorized.size(); i++){
            for(unsigned int j = 0; j < vectorized[0].size(); j++){
                temp.set(j,i, vectorized[i][j]);
            
            }
        }

        utils::print_matrix(temp);
        std::cout<<std::endl;
        temp = gaussian_elimination(temp);
        utils::print_matrix(temp);
        std::cout<<std::endl;

        std::vector< g::matrix > rtn; 
        int csd = 0;
        for(int i = 0; i < temp.cols(); i++){
            for(int j = csd; j < temp.rows(); j++){
                if(!temp(j,i).is_zero()){
                    csd++;
                    rtn.push_back(matrices[i]);
                    break;
                }
             
            }
            
        }
        return rtn;

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
                swap_rows(&mat_temp, numNonZeroCols, pivot_row);
                for (int j = numNonZeroCols+1; j < mat_temp.rows(); j++){
                    // Subtract the ith row from the jth row
                    g::ex temp_first_ex = mat_temp(j,i);
                    for (int k = i; k < mat_temp.cols(); k++){ 
                        mat_temp.set(j,k,(mat_temp(j,k) * mat_temp(numNonZeroCols,i)) - mat_temp(numNonZeroCols,k) * temp_first_ex);
                    }
                }

                numNonZeroCols += 1; // Increment the number of non-zero columns to use as a starting point for the next column
            }
        }

        return mat_temp;
    }

    int rank(g::matrix &matrix){
        // Compute the gaussian elimination of a matrix and count the number of non-zero rows
        g::matrix temp = gaussian_elimination(matrix);

        // Initialize the number of non-zero rows and the column value to start searching for non-zero values
        int rk = 0;
        int colVal = 0;
        
        // Go through each row starting from the column the previous row ended at and move until the next non-zero value is found.
        for(int i = 0; i < temp.rows(); i++){
            for(int j = colVal; j < temp.cols(); j++){
                if(!temp(i,j).is_zero()){
                    rk++;
                    colVal = j+1;
                    break;
                }
            }
        }

        return rk;
    }

    std::vector< g::exvector > nullspace(g::matrix &matrix){
        // Return a basis of the nullspace of the given matrix

        // Compute the gaussian elimination of the matrix
        g::matrix temp = gaussian_elimination(matrix);

        // Initialize the nullspace basis        
        std::vector< g::exvector > nullspace_basis;

        // The nullspace of the matrix is the columns after the pivots
        // minus e_i where i the i is the numNonZeroRows plus the number
        // of columns past the last pivot. We then multiply the result
        // such that there is no division.

        // Find the number of columns before the nullspace starts
        // and store the pivots and the columns they are in
        std::vector< g::ex > pivots;
        std::vector< int > pivot_cols;
        int colVal = 0;
        int rk = 0;
        for(int i = 0; i < temp.rows(); i++){
            for(int j = colVal; j < temp.cols(); j++){
                if(!temp(i,j).is_zero()){
                    rk++;
                    colVal = j+1;
                    pivot_cols.push_back(j);
                    pivots.push_back(temp(i,j));
                    break;
                }
            }
        }

        // Find the columns that are not pivots before the nullspace starts
        // and add them to the nullspace basis
        colVal = 0;
        for(int i = 0; i < rk; i++){
            for(int j = colVal; j < temp.cols(); j++){
                if(!temp(i,j).is_zero()){
                    colVal = j+1;
                    break;
                } else {
                    // Generate the basis vector corresponding to the column
                    g::exvector nullspace_vector;
                    for (int k = 0; k < temp.cols(); k++){
                        if (k == j){
                            nullspace_vector.push_back(1);
                        } else {
                            nullspace_vector.push_back(0);
                        }
                    }
                    // Add the basis vector to the nullspace basis
                    nullspace_basis.push_back(nullspace_vector);
                }
            }
        }

        // Compute the product of the pivots
        g::ex totalProduct = 1;
        for (g::ex pivot : pivots){
            totalProduct *= pivot;
        }

        for (int i = colVal; i < temp.cols(); i++){
            // Initialize the nullspace vector
            g::exvector nullspace_vector;
            std::vector< int >::iterator pivot_it = pivot_cols.begin();
            int pivot_ind = 0;
            for (int j = 0; j < temp.cols(); j++){
                // Add the correct multiplying factor times the element in the column
                if (j == *pivot_it && pivot_it != pivot_cols.end()){
                    nullspace_vector.push_back(temp(pivot_ind,i)*totalProduct/pivots[pivot_ind]);
                    pivot_it++; pivot_ind++;
                }
                // Add the negative of the total product
                else if (j == colVal){
                    nullspace_vector.push_back(-totalProduct);
                }
                // Add 0's
                else{
                    nullspace_vector.push_back(0);
                }
            }
            nullspace_basis.push_back(nullspace_vector);
        }

        return nullspace_basis;
    }

    g::ex prod_trace(g::matrix a, g::matrix b) {
        unsigned int num_rows = a.rows();
        unsigned int num_cols = a.cols();
        g::ex out = g::ex();
        for (int i = 0; i < num_rows; i++) {
            for (int j = 0; j < num_cols; j++) {
                out = g::add(out, g::mul(a(i,j),b(j,i)));
            }
        }
        return out;
    }
};