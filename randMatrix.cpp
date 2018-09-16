/* MIT License
 
 Copyright (c) 2018 Max Kapsecker
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

==============================================================================*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <ctime>
#include <random>
#include <map>

#define MAX_DISPLAY_SIZE 10
#define MATRIX_ENTRY_SPACE 12
#define SAMPLING_STRATEGY 2

namespace matrixOperations{
    
    inline void sumMatrix(float **A, float **B, float **C, int size_l, int size_m){
        
        /*
         Description:   - Sum of two matrices of the same size
         Runtime:       - O(l*m)
         Input:         - (float **A)   Reference on input matrix A
                        - (float **B)   Refrence on input matrix B
                        - (float **C)   Reference on output matrix C
                        - (int size_l)  row count of A resp. B
                        - (int size_m)  col count of A resp. B
         Output:        -  void
         */
        
        for(int i = 0; i < size_l; i++){
            for(int j = 0; j < size_m; j++){
                C[i][j] = A[i][j] + B[i][j];
            }
        }
    }

    inline void differenceMatrix(float **A, float **B, float **C, int size_l, int size_m){
        
        /*
         Description:   - Difference of two matrices of the same size
         Runtime:       - O(l*m)
         Input:         - (float **A)   Reference on input matrix A
                        - (float **B)   Refrence on input matrix B
                        - (float **C)   Reference on output matrix C
                        - (int size_l)  row count of A resp. B
                        - (int size_m)  col count of A resp. B
         Output:        - void
         */
        
        for(int i = 0; i < size_l; i++){
            for(int j = 0; j < size_m; j++){
                C[i][j] = A[i][j] - B[i][j];
            }
        }
    }


    inline float frobeniusNorm(float **A, int size_l, int size_m){
        
        /*
         Description:   - Compute Frobenius norm of input matrix
         Runtime:       - O(l*m)
         Input:         - (float **A)   Reference on input matrix A
                        - (int size_l)  row length of A
                        - (int size_m)  col length of B
         Output:        - (float)       Frobenius norm of input matrix
         */
        
        float norm = 0.0;
        for(int i = 0; i < size_l; i++){
            for(int j = 0; j < size_m; j++){
                norm += A[i][j] * A[i][j];
            }
        }
        
        return sqrt(norm);
    }


    float euclideanNormRows(float **A, int index, int size_l){
        
        /*
         Description:   - Compute euclidean norm of input vector
         Runtime:       - O(l)
         Input:         - (float **A)   Reference on input matrix
                        - (int index)   column index to choose from matrix
                        - (int size_l)  row length of A
         Output:        - (float)       Euclidean norm of input vector
         */
        
        float norm = 0.0;
        for(int i = 0; i < size_l; i++){
            norm += A[i][index] * A[i][index];
        }
        return sqrt(norm);
    }


    float euclideanNormCols(float **array, int index, int size_l){
        
        /*
         Description:   - Compute euclidean norm of input vector
         Runtime:       - O(l)
         Input:         - (float **A)   Reference on input matrix
                        - (int index)   row index to choose from matrix
                        - (int size_l)  col length of A
         Output:        - (float) Euclidean norm of input vector
         */
        
        float norm = 0.0;
        for(int i = 0; i < size_l; i++){
            norm += array[index][i] * array[index][i];
        }
        return sqrt(norm);
    }


    
    
    void subsampleMatrix(float **inputMatrix1, float **inputMatrix2, float **outputMatrix1, float **outputMatrix2, int l, int m, int n, int red_m, int sampling){
        
        /*
         Description:   - Subsample a matrix by rows resp. columns by either a uniform or a custom distribution
         Input:         - (float **inputMatrix1)    reference on original input matrix 1
                        - (float **inputMatrix1)    reference on original input matrix 2
                        - (float **inputMatrix1)    reference on subsampled outcome matrix
                        - (float **inputMatrix1)    reference on subsampled outcome matrix
                        - (int l)                   row count of input matrix 1
                        - (int m)                   col / row count of input matrix 1 / input matrix 2
                        - (int n)                   col count of input matrix 2
                        - (int red_m)               col / row count subsampled matrices
                        - (int sampling)            type of subsampling (uniform or custom (1 = uniform, 2 = see paper))
         Output:        - void
         */
        
        // TODO: Check if B is a transpose of A. If so, we can do much more efficient in space-complexity
        if(sampling == SAMPLING_STRATEGY){
            
            // Runtime: - O(l*red_m)
            float sum = 0.0;
            float *A_row = new float[m];
            float *B_col = new float[m];
            int *value = new int[red_m];
            
            // compute probability vector
            for(int i = 0; i < m; i++){
                A_row[i] = euclideanNormRows(inputMatrix1, i, l);
                B_col[i] = euclideanNormCols(inputMatrix2, i, n);
                sum += (A_row[i] * B_col[i]);
            }
            for(int i = 0; i < m; i++){
                A_row[i] = (A_row[i] * B_col[i]) / sum;
            }
            
            std::random_device rd;
            std::mt19937 eng(rd());
            std::discrete_distribution<> distr(A_row[0], A_row[m]);
            
            for(int i = 0; i < red_m; i++){
                value[i] = distr(eng);
                for(int j = 0; j < n; j++){
                    outputMatrix2[value[i]][j] = inputMatrix2[value[i]][j] / (A_row[value[i]] * red_m);
                }
            }
            for(int i = 0; i < l; i++){
                for(int j = 0; j < red_m; j++){
                    outputMatrix1[i][value[j]] = inputMatrix1[i][value[j]];
                }
            }
            
            delete[] A_row;
            delete[] B_col;
            delete[] value;

        }
        else {
            
            // Runtime: - O(l*red_m)
            int *value = new int[red_m];
            
            // uniform index sampling
            std::random_device rd;
            std::mt19937 eng(rd());
            std::uniform_int_distribution<> distr(0, red_m - 1);
            
            // downsample matrices
            for(int i = 0; i < red_m; i++){
                value[i] = int(distr(eng));
                for(int j = 0; j < n; j++){
                    outputMatrix2[value[i]][j] = inputMatrix2[value[i]][j];
                }
            }
            for(int i = 0; i < l; i++){
                for(int j = 0; j < red_m; j++){
                    outputMatrix1[i][value[j]] = inputMatrix1[i][value[j]];
                }
            }
            delete[] value;
        }
    }


    void outerMatrixProduct(float **A, float **B, float **C, int size_l, int size_m, int size_n){
        
        /*
         Description:   - Matrix multiplication AB by outer products
         Runtime:       - O(l*m*n)
         Input:         - (float **A)   Reference on input matrix A
                        - (float **B)   Refrence on input matrix B
                        - (float **C)   Reference on output matrix C
                        - (int size_l)  row length of A
                        - (int size_m)  row and col length of A and B
                        - (int size_n)  col length of B
         Output:        - void
         */
        
        for(int k = 0; k < size_m; k++){
            for(int i = 0; i < size_l; i++){
                for(int j = 0; j < size_n; j++){
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
    }


    void innerMatrixProduct(float **A, float **B, float **C, int size_l, int size_m, int size_n){
        
        /*
         Description:   - Classical matrix multiplication AB by inner products
         Runtime:       - O(l*m*n)
         Input:         - (float **A)   Reference on input matrix A
                        - (float **B)   Refrence on input matrix B
                        - (float **C)   Reference on output matrix C
                        - (int size_l)  row length of A
                        - (int size_m)  row and col length of A and B
                        - (int size_n)  col length of B
         Output:        - void
         */
        
        for(int i = 0; i < size_l; i++){
            for(int j = 0; j < size_n; j++){
                C[i][j] = 0.0;
                for(int k = 0; k < size_m; k++){
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }
    }


    void printMatrix(float **matrix, int size_l, int size_m){
        
        /*
         Description:   - Prints a matrix
         Input:         - (float **matrix)  Reference on matrix
                        - (int size_l)      Row count of input matrix A
                        - (int size_m)      Col count of input matrix A
         Output:        - void
         */
        
        std::cout << std::endl;
        for(int i = 0; i < size_l; i++){
            std::cout << std::setw(MATRIX_ENTRY_SPACE);
            for (int j = 0; j < size_m; j++){
                std::cout << matrix[i][j] << std::setw(MATRIX_ENTRY_SPACE);
            }
            std::cout << std::endl;
        }
    }


    void fillMatrix(float **matrix, int size_l, int size_m){
        
        /*
         Description:   - Fill input matrix with random numbers
         Runtime:       - O(l*m)
         Input:         - (**matrix)    Reference on input matrix A
                        - (int size_l)  Row count of input matrix A
                        - (int size_m)  Col count of input matrix A
         Output:        - Matrix with random numbers (range: [0,1]) as entry
         */
        
        // gaussian
        /*std::random_device rd{};
        std::mt19937 gen{rd()};
        std::normal_distribution<> d{10,1};
        std::map<int, int> hist{};
        
        for(int i = 0; i < size_l; i++){
            for (int j = 0; j < size_m; j++){
                matrix[i][j] = double(std::round(d(gen)));
            }
        }*/
        
        // uniform in [0,1]
        for(int i = 0; i < size_l; i++){
            for (int j = 0; j < size_m; j++){
                matrix[i][j] = double(rand()) / (RAND_MAX);
            }
        }
    }
}

int main(){
    

    bool consoleMode;
    std::cout << "Row of experiments (type: 0) or singular experiment (type: 1)?" << std::endl;
    std::cin >> consoleMode;
    std::cout << "------------------------------" << std::endl;
    
    // key values
    double elapsed_secs_inner_exact, elapsed_secs_outer_exact, elapsed_secs_outer_approx;
    float error_outer_exact, error_outer_approx;
    
    // set seed
    srand(time(NULL));
    
    if(consoleMode){
        
        int l, m, n, red_m, subsampling;
        float subsample;
        std::cout << "Rowcount of matrix A: ";
        std::cin >> l;
        std::cout << "Row resp. columncount of matrix A resp. B: ";
        std::cin >> m;
        std::cout << "Columncount of matrix B: ";
        std::cin >> n;
        std::cout << "Factor of subsampling: ";
        std::cin >> subsample;
        red_m = ceil(m*subsample);
        std::cout << "Using " << red_m << " rows/cols in the probabilistic product instead of " << m << "." << std::endl;
        std::cout << "Use uniform sampling (1) or custom sampling (2): ";
        std::cin >> subsampling;
        std::cout << "------------------------------" << std::endl;
        
        
        // allocate
        float **A, **B, **A_red, **B_red, **C_exact_inner, **C_exact_outer, **C_approx_outer;
        
        // original matrices
        A = new float *[l];
        B = new float *[m];
        
        // subsampled matrices
        A_red = new float *[l];
        B_red = new float *[red_m];
        
        // outcome matrices
        C_exact_inner = new float *[l];
        C_exact_outer = new float *[l];
        C_approx_outer = new float *[l];
        
        for(int i = 0; i < l; i++){
            A[i] = new float[m];
            A_red[i] = new float[red_m];
            C_exact_inner[i] = new float[n];
            C_exact_outer[i] = new float[n];
            C_approx_outer[i] = new float[n];
        }
        for(int i = 0; i <m; i++){
            B[i] = new float[n];
        }
        for(int i = 0; i < red_m; i++){
            B_red[i] = new float[n];
        }
        
        // fill matrices with random number
        matrixOperations::fillMatrix(A, l, m);
        matrixOperations::fillMatrix(B, m, n);
        
        // print if suitable
        if(l <= MAX_DISPLAY_SIZE && m <= MAX_DISPLAY_SIZE && n <= MAX_DISPLAY_SIZE){
            std::cout << "A = ";
            matrixOperations::printMatrix(A, l, m);
            std::cout << std::endl;
            std::cout << "B = ";
            matrixOperations::printMatrix(B, m, n);
            std::cout << std::endl;
        }
        else{
            std::cout << "A = " << std::endl;
            std::cout << "To large to display." << std::endl;
            std::cout << "B = " << std::endl;
            std::cout << "To large to display." << std::endl;
        }
        std::cout << "------------------------------" << std::endl;
        
        // EXACT INNER MATRIX PRODUCT
        clock_t begin_inner_exact = clock();
        matrixOperations::innerMatrixProduct(A, B, C_exact_inner, l, m, n);
        clock_t end_inner_exact = clock();
        elapsed_secs_inner_exact = (double(end_inner_exact - begin_inner_exact) / CLOCKS_PER_SEC) * 1000.0;
        std::cout << "INNER MATRIX PRODUCT" << std::endl << std::endl;
        std::cout << "A*B = ";
        if(l <= MAX_DISPLAY_SIZE && n <= MAX_DISPLAY_SIZE){
            matrixOperations::printMatrix(C_exact_inner, l, n);
            std::cout << std::endl;
        }
        else{
            std::cout << "To large to display." << std::endl;
        }
        std::cout << "Computed the exact inner matrix product in " << elapsed_secs_inner_exact << " milliseconds" << std::endl;
        std::cout << "------------------------------" << std::endl;
        
        // EXACT OUTER MATRIX PRODUCT
        clock_t begin_outer_exact = clock();
        matrixOperations::outerMatrixProduct(A, B, C_exact_outer, l, m, n);
        clock_t end_outer_exact = clock();
        elapsed_secs_outer_exact = (double(end_outer_exact - begin_outer_exact) / CLOCKS_PER_SEC) * 1000.0;
        std::cout << "OUTER MATRIX PRODUCT" << std::endl << std::endl;
        std::cout << "A*B = ";
        if(l <= MAX_DISPLAY_SIZE && n <= MAX_DISPLAY_SIZE){
            matrixOperations::printMatrix(C_exact_outer, l, n);
            std::cout << std::endl;
        }
        else{
            std::cout << "To large to display." << std::endl;
        }
        std::cout << "Computed the exact product in " << elapsed_secs_outer_exact << " milliseconds" << std::endl;
        matrixOperations::differenceMatrix(C_exact_inner, C_exact_outer, C_exact_outer, l, n);
        error_outer_exact = matrixOperations::frobeniusNorm(C_exact_outer, l, n);
        std::cout << "Error of the exact outer matrix product: " << error_outer_exact << ". Surprise ;-)" << std::endl;
        std::cout << "------------------------------" << std::endl;
        
        // APPROXIMATE OUTER MATRIX PRODUCT - UNIFORM SAMPLING
        clock_t begin_outer_approx = clock();
        matrixOperations::subsampleMatrix(A, B, A_red, B_red, l, m, n, red_m, subsampling);
        matrixOperations::outerMatrixProduct(A_red, B_red, C_approx_outer, l, red_m, n);
        clock_t end_outer_approx = clock();
        elapsed_secs_outer_approx = (double(end_outer_approx - begin_outer_approx) / CLOCKS_PER_SEC) * 1000.0;
        std::cout << "APPROXIMATE OUTER MATRIX PRODUCT" << std::endl << std::endl;
        std::cout << "C*R = ";
        if(l <= MAX_DISPLAY_SIZE && n <= MAX_DISPLAY_SIZE){
            matrixOperations::printMatrix(C_approx_outer, l, n);
            std::cout << std::endl;
        }
        else{
            std::cout << "To large to display." << std::endl;
        }
        std::cout << "Computed the approximate outer matrix product in " << elapsed_secs_outer_approx << " milliseconds" << std::endl;
        matrixOperations::differenceMatrix(C_exact_inner, C_approx_outer, C_approx_outer, l, n);
        error_outer_approx = matrixOperations::frobeniusNorm(C_approx_outer, l, n);
        std::cout << "Error of the approximate outer product: " << error_outer_approx << "." << std::endl;
        std::cout << "------------------------------" << std::endl;
        
        // CLEAN DYNAMIC MEMORY
        // FREE SUBARRAYS
        for(int i = 0; i <l; i++){
            delete[] A[i];
            delete[] A_red[i];
            delete[] C_exact_inner[i];
            delete[] C_exact_outer[i];
            delete[] C_approx_outer[i];
        }
        for(int i = 0; i < m; i++){
            delete[] B[i];
        }
        for(int i = 0; i < red_m; i++){
            delete[] B_red[i];
        }
        // FREE ARRAYS
        delete[] A;
        delete[] B;
        delete[] A_red;
        delete[] B_red;
        delete[] C_exact_inner;
        delete[] C_exact_outer;
        delete[] C_approx_outer;
    }
    
    else{
    
        int a, b, c, l, m, n, red_m, subsampling;
        float subsample;
        
        std::cout << "We consider squared matrices." << std::endl;
        std::cout << "Set start dimension: ";
        std::cin >> a;
        std::cout << "Set end dimension: ";
        std::cin >> b;
        std::cout << "Set increment: ";
        std::cin >> c;
        std::cout << "Factor of subsampling: ";
        std::cin >> subsample;
        std::cout << std::endl;
        std::cout << "Use uniform sampling (1) or custom sampling (2): ";
        std::cin >> subsampling;
        std::cout << "------------------------------" << std::endl;
        
        std::ofstream result;
        result.open ("result.txt");
        result << "matrixDimension,timeInnerExact,errorInnerExact,timeOuterExact,errorOuterExact,timeOuterApprox,errorOuterApprox,subsample,subsampling\n";
        
        for(int k = a; k <= b; k += c){
            
            l = k;
            m = k;
            n = k;
            
            red_m = ceil(m*subsample);
            
            // allocate
            float **A, **B, **A_red, **B_red, **C_exact_inner, **C_exact_outer, **C_approx_outer;
            
            // original matrices
            A = new float *[l];
            B = new float *[m];
            
            // subsampled matrices
            A_red = new float *[l];
            B_red = new float *[red_m];
            
            // outcome matrices
            C_exact_inner = new float *[l];
            C_exact_outer = new float *[l];
            C_approx_outer = new float *[l];
            
            for(int i_1 = 0; i_1 < l; i_1++){
                A[i_1] = new float[m];
                A_red[i_1] = new float[red_m];
                C_exact_inner[i_1] = new float[n];
                C_exact_outer[i_1] = new float[n];
                C_approx_outer[i_1] = new float[n];
            }
            for(int i_2 = 0; i_2 <m; i_2++){
                B[i_2] = new float[n];
            }
            for(int i_3 = 0; i_3 < red_m; i_3++){
                B_red[i_3] = new float[n];
            }
            
            // fill matrix randomly
            matrixOperations::fillMatrix(A, l, m);
            matrixOperations::fillMatrix(B, m, n);
            
            // inner product exact
            clock_t begin_inner_exact = clock();
            matrixOperations::innerMatrixProduct(A, B, C_exact_inner, l, m, n);
            clock_t end_inner_exact = clock();
            elapsed_secs_inner_exact = (double(end_inner_exact - begin_inner_exact) / CLOCKS_PER_SEC) * 1000.0;
            
            // outer product exact
            clock_t begin_outer_exact = clock();
            matrixOperations::outerMatrixProduct(A, B, C_exact_outer, l, m, n);
            clock_t end_outer_exact = clock();
            elapsed_secs_outer_exact = (double(end_outer_exact - begin_outer_exact) / CLOCKS_PER_SEC) * 1000.0;
            matrixOperations::differenceMatrix(C_exact_inner, C_exact_outer, C_exact_outer, l, n);
            error_outer_exact = matrixOperations::frobeniusNorm(C_exact_outer, l, n);
            
            // outer product approx - uniform sampling
            clock_t begin_outer_approx = clock();
            matrixOperations::subsampleMatrix(A, B, A_red, B_red, l, m, n, red_m, subsampling);
            matrixOperations::outerMatrixProduct(A_red, B_red, C_approx_outer, l, red_m, n);
            clock_t end_outer_approx = clock();
            elapsed_secs_outer_approx = (double(end_outer_approx - begin_outer_approx) / CLOCKS_PER_SEC) * 1000.0;
            matrixOperations::differenceMatrix(C_exact_inner, C_approx_outer, C_approx_outer, l, n);
            error_outer_approx = matrixOperations::frobeniusNorm(C_approx_outer, l, n);
            
            
            // write results to file
            result << std::fixed << l << "," << elapsed_secs_inner_exact << "," << 0.0 << "," << elapsed_secs_outer_exact << "," << error_outer_exact << "," << elapsed_secs_outer_approx << "," << double(error_outer_approx) << "," << subsample << "," << subsampling << "\n";
            
            // clean dynamic memory
            // free sub-arrays
            // TODO: Manage to free memory within loop! IMPORTANT!
            /*for(int i1 = 0; i1 < l; i1++){
                delete[] A[i1];
                delete[] A_red[i1];
                delete[] C_exact_inner[i1];
                delete[] C_exact_outer[i1];
                delete[] C_approx_outer[i1];
            }
            for(int i2 = 0; i2 < m; i2++){
                delete[] B[i2];
            }
            for(int i3 = 0; i3 < red_m; i3++){
                delete[] B_red[i3];
            }
            // free arrays
            delete[] A;
            delete[] B;
            delete[] A_red;
            delete[] B_red;
            delete[] C_exact_inner;
            delete[] C_exact_outer;
            delete[] C_approx_outer;
            
            A = nullptr;
            B = nullptr;
            A_red = nullptr;
            B_red = nullptr;
            C_exact_inner = nullptr;
            C_exact_outer = nullptr;
            C_approx_outer = nullptr;*/
            
        }
        
        // close file
        result.close();
        
    }
	return 0;

}
