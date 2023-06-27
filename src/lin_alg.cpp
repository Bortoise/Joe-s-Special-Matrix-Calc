#include "headers/lin_alg.h"
#include "headers/utils.h"


namespace lin_alg{

    g::matrix basis_to_vectorized_matrix(mat_vec &basis) {
        if (basis.empty()) {
            return {0,0};
        }
        g::matrix out = {static_cast<unsigned int>(basis.size()),
                         static_cast<unsigned int>(basis[0].nops())};
        int i = 0;
        for (g::matrix m : basis) {
            int j = 0;
            for (const g::ex& e : vectorize(m)) {
                out(i,j) = e;
                j++;
            }
            i++;
        }
        return out;
    }

    g::exvector vectorize(g::matrix &m){
        g::exvector v;
        for(int i = 0; i < m.nops(); i++){
            v.push_back(m.op(i));
        }

        return v;

    }

    g::matrix matricize(g::exvector &v, unsigned int r, unsigned int c, bool rowcol){
        g::matrix x(r,c);
        if (rowcol == false){
            for(int i  = 0; i < r; i++){
                for(int j = 0; j < c; j++){
                    x.set(i,j, v[i*c + j]);
                }
            }
        } else {
            for(int i  = 0; i < c; i++){
                for(int j = 0; j < r; j++){
                    x.set(j,i, v[i*r + j]);
                }
            }
        }
        return x;
    }

    g::matrix matricize(g::matrix &v, unsigned int r, unsigned int c, bool rowcol){
        g::matrix x(r,c);
        if (rowcol == false){
            for(int i  = 0; i < r; i++){
                for(int j = 0; j < c; j++){
                    x.set(i,j, v(i*c + j,0));
                }
            }
        } else {
            for(int i  = 0; i < c; i++){
                for(int j = 0; j < r; j++){
                    x.set(j,i, v(i*r + j,0));
                }
            }
        }
        return x;
    }

    g::matrix vector_to_matrix(g::exvector &v, mat_vec &basis){
        // Convert a vector to a matrix in a given basis
        if(v.size() != basis.size()){
            throw("Error: vector size does not match basis size");
        }

        g::matrix rtn(basis[0].rows(), basis[0].cols());
        for(int i = 0; i < v.size(); i++){
            rtn = rtn.add(basis[i].mul_scalar(v[i]));
        }
        return rtn;
    }

    g::matrix vector_to_matrix(g::matrix &v, mat_vec &basis){
        // Convert a vector to a matrix in a given basis
        if(v.nops() != basis.size()){
            throw("Error: vector size does not match basis size");
        }

        g::matrix rtn(basis[0].rows(), basis[0].cols());
        for(int i = 0; i < v.nops(); i++){
            rtn = rtn.add(basis[i].mul_scalar(v(i,0)));
        }
        return rtn;
    }

    g::matrix bracket(g::matrix &x, g::matrix &y){
        // Compute [x,y] = xy - yx
        return x.mul(y).sub(y.mul(x)); 
    }

    mat_vec spanning_subsequence(mat_vec &matrices){
        if(matrices.empty()) {
            return mat_vec();
        }
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

        temp = gaussian_elimination(temp);
        mat_vec rtn; 
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

    void swap_rows(g::matrix &mat, int row1, int row2){
        // Swap two rows in a matrix
        for (int i = 0; i < mat.cols(); i++){
            g::ex temp = mat(row1,i);
            mat(row1,i) = mat(row2,i);
            mat(row2,i) = temp;
        }
    }

    g::matrix gaussian_elimination(g::matrix const &mat){
        // Perform gaussian elimination on a matrix without swapping any columns

        // Make a deep copy of the matrix to perform operations on
        g::matrix mat_temp = g::matrix(mat.rows(), mat.cols());    

        // Copies the matrix
        for (int i = 0; i < mat.rows(); i++){
            for (int j = 0; j < mat.cols(); j++){
                mat_temp(i,j) = mat(i,j);
            }
        }

        // Initialize the pivot and pivot row
        g::ex pivot;
        int pivot_row;

        // Initialize the number of non-zero columns, where to search for the pivot
        int numNonZeroCols = 0;

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
                swap_rows(mat_temp, numNonZeroCols, pivot_row); // SUSPICIOUS 
                for (int j = 0; j < mat_temp.rows(); j++){
                    if (j == numNonZeroCols) {
                        continue;
                    }
                    g::ex temp_first_ex = mat_temp(j,i);
                    // Subtracts an appropriate multiple of the numNonZeroCols-th row from an appropriate multiple of the jth row
                    if (!temp_first_ex.is_zero()){
                        for (int k = 0; k < mat_temp.cols(); k++){ 
                            mat_temp.set(j,k,(mat_temp(j,k) * pivot) - (mat_temp(numNonZeroCols,k) * temp_first_ex));
                        }
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
        int col = 0;
        int row = 0;

        while (col < temp.cols() && row < temp.rows()) {
            if (!temp(row, col).is_zero()) row++;
            col++;
        }

        return row;
    }

    mat_vec nullspace(g::matrix &matrix, bool division){
        // Return a basis of the nullspace of the given matrix

        // Compute the gaussian elimination of the matrix
        g::matrix temp = gaussian_elimination(matrix);

        // Initialize the nullspace basis        
        mat_vec nullspace_basis;

        // Go through the matrix, store the pivot columns and for each non-pivot column, add the corresponding vector to the nullspace basis
        std::vector< int > pivot_columns;
        g::exvector pivots;

        int col = 0;
        int row = 0;

        while (col < temp.cols() && row < temp.rows()) {
            if(!temp(row, col).normal().is_zero()){
                pivot_columns.push_back(col);
                pivots.push_back(temp(row,col));
                row++;
            }
            col++;
        }
        
        g::exvector total_product = {};
        for (int i = 0; i < pivots.size(); i++){
            g::ex temp_total_product = 1;
            if (!division) {
                for (int j = 0; j < pivots.size(); j++){
                    if (i != j){
                        temp_total_product = temp_total_product * pivots[j];
                    }
                }
            } else {
                temp_total_product = 1/pivots[i];
            }
            total_product.push_back(temp_total_product);
        }
        
        int num_current_pivots = 0;
        for (int i = 0; i < temp.cols(); i++){
            if (num_current_pivots < pivot_columns.size() && i == pivot_columns[num_current_pivots]){
                num_current_pivots++;
                continue;
            }

            g::matrix nullspace_vector = g::matrix(temp.cols(), 1);
            for (int j = 0; j < num_current_pivots; j++){
                nullspace_vector(pivot_columns[j],0) = -(temp(j,i) * total_product[j]).normal();
            }
            nullspace_vector(i,0) = 1;

            nullspace_basis.push_back(nullspace_vector);
        }
        
        return nullspace_basis;
    }
    
    g::exvector sl_ize(g::matrix m, int n) {
        g::exvector out = {}; 
        for (int i = 0; i < n-1; i++) {
            for (int j = i+1; j < n; j++) {
                out.push_back(m(i,j));
                out.push_back(m(j,i));
            }
        }
        for (int i = 1; i < n; i++) {
            out.push_back(-1*m(i,i));
        }
        return out;
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

    
    bool contains(std::vector< g::exvector >* A, std::vector< g::exvector >* B) {
        if (A->size() > B->size() ) {
            return false;
        }
        std::vector< g::exvector > vectors = {};

        vectors.insert(vectors.end(), A->begin(), A->end());
        vectors.insert(vectors.end(), B->begin(), B->end());

        mat_vec vectors2 = {};
        for (g::exvector v : vectors) {
            vectors2.push_back(matricize(v,1,v.size()));
        }
        // Makes a matrix whose rows are the elements of A and B.
        g::matrix M = lin_alg::basis_to_vectorized_matrix(vectors2);

        // span(B) is a subset of span(A) iff dim(span(A))=dim(span(L)+span(N)), but that is to say rank(e_1,...,e_n)
        // = rank(S). Taking transposes this is what we actually verify.
        return A->size() == lin_alg::rank(M);
    }

    bool equals(std::vector< g::exvector >* A, std::vector< g::exvector >* B){
        if (A->size() != B->size() ) {
            return false;
        }

        return contains(A, B) && contains(B, A);
    }
};