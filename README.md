# Randomized Matrix Product
Let us consider the matrix product of A and B (Here: A has l rows and m columns, while B has m rows and n columns).
It is well known that the classical approach (3 for-loop solutions) for the matrix product can be estimated by a runtime of O(lmn). This runtime makes the three-loop-algorithm a computational inefficient approach for large-scale matrices A and B. Theoretical considerations have been accomplished, to reduce the number of FLOPs. One of these considerations is to use a probabilistic method. That is a random subsampling (according to a specified probability distribution) of the rows in A and according to columns in B to the reduced matrices A' and B'. By representing the multiplication concerning the dyadic product, it becomes apparent how the supposed dimensions (mxn) of the approximate outcome matrix A'B' is assured. Note: Even though the size/dimension of the resulting matrix is the same as for the exact product, the rank of A'B' is likely to be lower. The question remains how to subsample the rows in A and columns in B. Therefore we state a minimization problem to reduce the expected difference of the exact matrix product and the approximated matrix product measured regarding the Frobenius norm. Following link provides a more detailed and theoretical explanation (see Chapter 2: Pages 18 - 22):

https://arxiv.org/pdf/1608.04481v1.pdf

The original work can be found in:

https://www.stat.berkeley.edu/~mmahoney/pubs/matrix1_SICOMP.pdf

The chart below draws the outcome of the program. As expected we have a lower runtime by considering the randomized matrix multiplication. The average error on each cell represents the accompanying error for the probabilistic approach, that is the total error measured by the Frobenius norm divided by the number of total entries. For instance: Let's assume an average error of 0.1 measured for a matrix dimension of 25x25. The overall error is then reverse computed by 0.1*25*25 = 62.5.

![alt text](https://github.com/NumericalMax/RandomizedMatrixProduct/blob/master/images/all_avg_Uniform.png)
![alt text](https://github.com/NumericalMax/RandomizedMatrixProduct/blob/master/images/all_avg_Custom.png)

# Application
There are applications, in which an approximate solution can replace the exact answer. The critical part of such claims is to find the tradeoff between error size and runtime.

# Structure of the Repository
- images: Contains png-files drawing key-values of the considered algorithms
- evaluate.ipynb: Jupyter notebook for creating the png-files based on the outcome of the code (see: result.txt)
- makefile.make: Makefile to compile (gcc is used) to code.
- randMatrix.cpp: C++ code for performing the different approaches
- result.txt: Automatically generated output file from the C++ code.

# Machine
The Code was solely implemented and tested on a MacBook Pro i5 / 8GB RAM. Due to background processes, the resulting time graphs look volatile.

# Dependencies of the Code
iostream, iomanip, fstream, random, ctime, map, math.h

# Execution
The compilation / execution of the program might look like this:
![alt text](https://github.com/NumericalMax/RandomizedMatrixProduct/blob/master/images/run_0.png)
![alt text](https://github.com/NumericalMax/RandomizedMatrixProduct/blob/master/images/run_1.png)

# Compiler
Note: The g++ compiler with the c++11 standard is used.

# TODO:
- Free dynamic memory in for the loop! 

# Author
Max Kapsecker, 2018: Remark, that there is always room for optimization. The Code can indeed be optimized regarding time and space efficiency. Feel free to report any mistake to max.kapsecker@tum.de
