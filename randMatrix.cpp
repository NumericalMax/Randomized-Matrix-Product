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
#include <ctime>
#include <math.h>
#include <fstream>

using namespace std;

void sumMatrix(float **A, float **B, float **C, int size_l, int size_m){
    
    /*
     Description:   - Sum of two matrices of the same size
     Runtime:       - O(l*m)
     Input:         - Reference on input matrix A
                    - Refrence on input matrix B
                    - Reference on output matrix C
     Output:        - Matrix difference
     */
    
    for(int i = 0; i < size_l; i++){
        for(int j = 0; j < size_m; j++){
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}

void differenceMatrix(float **A, float **B, float **C, int size_l, int size_m){
    
    /*
     Description:   - Subtraction of two matrices of the same size
     Runtime:       - O(l*m)
     Input:         - Reference on input matrix A
                    - Refrence on input matrix B
                    - Reference on output matrix C
     Output:        - Matrix difference
     */
    
    for(int i = 0; i < size_l; i++){
        for(int j = 0; j < size_m; j++){
            C[i][j] = A[i][j] - B[i][j];
        }
    }
}


float frobeniusNorm(float **A, int size_l, int size_m){
    
    /*
     Description:   - Compute Frobenius norm of input matrix
     Runtime:       - O(l*m)
     Input:         - Reference on input matrix A
                    - Refrence on input matrix B
                    - Reference on output matrix C
                    - row length of A
                    - row and col length of A and B
                    - col length of B
     Output:        - Frobenius norm of input matrix (float)
     */
    
    float norm = 0.0;
	for(int i = 0; i < size_l; i++){
		for(int j = 0; j < size_m; j++){
			norm += A[i][j] * A[i][j];
        }
	}
	return sqrt(norm);
}


void subsampleIndex(int *index, float ratio, int size_m){
    
    /*
     Description:   - Uniformly subsample an array of indices
     Runtime:       - O(m)
     Input:         - Index set
                    - Size of the index set
                    - Amount of indices which shall be sampled in percent
     Output:        - Subsampled index set
     */
    
    

}


void outerMatrixProduct(float **A, float **B, float **C, int size_l, int size_m, int size_n){
    
    /*
     Description:   - Matrix multiplication by outer products
     Runtime:       - O(l*m*n)
     Input:         - Reference on input matrix A
                    - Refrence on input matrix B
                    - Reference on output matrix C
                    - row length of A
                    - row and col length of A and B
                    - col length of B
     Output:        - Matrix Product
     */
    
    float **C_temp;
    C_temp = new float *[size_l];
    
    //for(int i = 0; i < size_l; i++){
    //    C_temp[i] = new float[size_n];
    //}
    
    for(int k = 0; k < size_m; k++){
        for(int i = 0; i < size_l; i++){
            C_temp[i] = new float[size_n];
            for(int j = 0; j < size_n; j++){
                C_temp[i][j] = A[i][k] * B[k][j];
            }
        }
        sumMatrix(C, C_temp, C, size_l, size_n);
    }
}


void innerMatrixProduct(float **A, float **B, float **C, int size_l, int size_m, int size_n){
    
    /*
     Description:   - Classical matrix multiplication by inner products
     Runtime:       - O(l*m*n)
     Input:         - Reference on input matrix A
                    - Refrence on input matrix B
                    - Reference on output matrix C
                    - row length of A
                    - row and col length of A and B
                    - col length of B
     Output:        - Matrix product
     */
    
    for(int i = 0; i < size_l; i++){
        for(int j = 0; j < size_n; j++){
            for(int k = 0; k < size_m; k++){
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}


void printMatrix(float **matrix, int size_l, int size_m){
	
    /*
     Description:   - Prints a matrix
     Input:         - Reference on input matrix A
                    - Row count of input matrix A
                    - Col count of input matrix A
     Output:        - Prints out matrix values
     */
    
    cout << endl;
    for(int i = 0; i < size_l; i++){
        for (int j = 0; j < size_m; j++){
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}


void fillMatrix(float **matrix, int size_l, int size_m){
    
    /*
     Description:   - Fill input matrix with random numbers
     Runtime:       - O(l*m)
     Input:         - Reference on input matrix A
     Output:        - Matrix with random numbers (range: [0,1]) as entry
     */
    
	for(int i = 0; i < size_l; i++){
		for (int j = 0; j < size_m; j++){
			
            matrix[i][j] = double(rand()) / (RAND_MAX);
        
        }
	}
}


int main(){
    
    bool consoleMode;
    cout << "Row of experiments (type: 0) or singular experiment (type: 1)?" << endl;
    cin >> consoleMode;
    cout << "------------------------------" << endl;
    
    // enable randomness
    srand(time(NULL));
    
    if(consoleMode){
    
        int l, m, n;
        float subsample;
        cout << "Rowcount of matrix A: ";
        cin >> l;
        cout << "Row resp. columncount of matrix A resp. B: ";
        cin >> m;
        cout << "Columncount of matrix B: ";
        cin >> n;
        cout << "Factor of subsampling: ";
        cin >> subsample;
        cout << endl;
        cout << "------------------------------" << endl;
        

        // allocate
        float **A, **B, **C_exact_inner, **C_exact_outer, **C_approx;
        A = new float *[l];
        B = new float *[m];
        C_exact_inner = new float *[l];
        C_exact_outer = new float *[l];
        C_approx = new float *[l];
        for(int i = 0; i <l; i++){
            A[i] = new float[m];
            C_exact_inner[i] = new float[n];
            C_exact_outer[i] = new float[n];
            C_approx[i] = new float[n];
        }
        for(int i = 0; i <m; i++){
            B[i] = new float[n];
        }
        
        
        // fill matrices with random number
        fillMatrix(A, l, m);
        fillMatrix(B, m, n);
        
        // print if suitable
        if(l < 50 && m < 50 && n < 50){
            cout << "A = ";
            printMatrix(A, l, m);
            cout << endl;
            cout << "B = ";
            printMatrix(B, m, n);
            cout << endl;
        }
        else{
            cout << "A = " << endl;
            cout << "To large to display." << endl;
            cout << "B = " << endl;
            cout << "To large to display." << endl;
        }
        cout << "------------------------------" << endl;
        
        // EXACT INNER MATRIX PRODUCT
        clock_t begin_inner_exact = clock();
        innerMatrixProduct(A, B, C_exact_inner, l, m, n);
        clock_t end_inner_exact = clock();
        double elapsed_secs_inner_exact = (double(end_inner_exact - begin_inner_exact) / CLOCKS_PER_SEC) * 1.0;
        cout << "INNER MATRIX PRODUCT" << endl << endl;
        cout << "A*B = ";
        if(l < 50 && n < 50){
            printMatrix(C_exact_inner, l, n);
            cout << endl;
        }
        else{
            cout << "To large to display." << endl;
        }
        cout << "Computed the exact product in " << elapsed_secs_inner_exact << " seconds" << endl;
        cout << "------------------------------" << endl;
        
        // EXACT OUTER MATRIX PRODUCT
        clock_t begin_outer_exact = clock();
        outerMatrixProduct(A, B, C_exact_outer, l, m, n);
        clock_t end_outer_exact = clock();
        double elapsed_secs_outer_exact = (double(end_outer_exact - begin_outer_exact) / CLOCKS_PER_SEC) * 1.0;
        cout << "OUTER MATRIX PRODUCT" << endl << endl;
        cout << "A*B = ";
        if(l < 50 && n < 50){
            printMatrix(C_exact_outer, l, n);
            cout << endl;
        }
        else{
            cout << "To large to display." << endl;
        }
        cout << "Computed the exact product in " << elapsed_secs_outer_exact << " seconds" << endl;
        differenceMatrix(C_exact_inner, C_exact_outer, C_exact_inner, l, n);
        float error = frobeniusNorm(C_exact_inner, l, n);
        cout << "Error of the exact outer product: " << error << ". Surprise ;-)" << endl;
        cout << "------------------------------" << endl;
        
        // APPROXIMATE OUTER MATRIX PRODUCT
        // ...
        
        
        //delete A, B, C_exact_inner, C_exact_outer;
    }
    
    else{
    
        int a, b, c, l, m, n;
        float subsample;
        cout << "We consider squared matrices." << endl;
        cout << "Set start dimension: ";
        cin >> a;
        cout << "Set end dimension: ";
        cin >> b;
        cout << "Set increment: ";
        cin >> c;
        cout << "Factor of subsampling: ";
        cin >> subsample;
        cout << endl;
        cout << "------------------------------" << endl;
        
        double elapsed_secs_inner_exact, elapsed_secs_outer_exact;
        float error_inner_exact, error_outer_exact, error_outer_approx;
        
        ofstream result;
        result.open ("result.txt");
        result << "matrixDimension,timeInnerExact,errorInnerExact,timeOuterExact,errorOuterExact,timeOuterApprox,errorOuterApprox,subsample\n";
        
        for(int k = a; k <= b; k = k + c){
        
            l = k;
            m = k;
            n = k;
            
            // allocate
            float **A, **B, **C_exact_inner, **C_exact_outer, **C_approx_outer;
            A = new float *[l];
            B = new float *[m];
            C_exact_inner = new float *[l];
            C_exact_outer = new float *[l];
            C_approx_outer = new float *[l];
            for(int i = 0; i <l; i++){
                A[i] = new float[m];
                C_exact_inner[i] = new float[n];
                C_exact_outer[i] = new float[n];
                C_approx_outer[i] = new float[n];
            }
            for(int i = 0; i <m; i++){
                B[i] = new float[n];
            }
            
            // fill
            fillMatrix(A, l, m);
            fillMatrix(B, m, n);
            
            // inner product exact
            clock_t begin_inner_exact = clock();
            innerMatrixProduct(A, B, C_exact_inner, l, m, n);
            clock_t end_inner_exact = clock();
            elapsed_secs_inner_exact = (double(end_inner_exact - begin_inner_exact) / CLOCKS_PER_SEC);
            error_inner_exact = 0;
            
            // outer product exact
            clock_t begin_outer_exact = clock();
            outerMatrixProduct(A, B, C_exact_outer, l, m, n);
            clock_t end_outer_exact = clock();
            elapsed_secs_outer_exact = (double(end_outer_exact - begin_outer_exact) / CLOCKS_PER_SEC) * 1.0;
            differenceMatrix(C_exact_inner, C_exact_outer, C_exact_outer, l, n);
            error_outer_exact = frobeniusNorm(C_exact_outer, l, n);
            result << l << "," << elapsed_secs_inner_exact << "," << error_inner_exact << "," << elapsed_secs_outer_exact << "," << error_outer_exact << "," << elapsed_secs_outer_exact << "," << error_outer_exact << "," << subsample << "\n";
            
            // outer product approx
            
            //delete A, B, C_exact_inner, C_exact_outer;
        }

        result.close();
        
    }
	return 0;

}
