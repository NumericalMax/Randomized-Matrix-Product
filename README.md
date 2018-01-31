# Randomized Matrix Product
Let us consider the matrix product of A and B (Here: A has l rows and m columns, while B has m rows and n columns).
It is well known, that the classical approach (3 for-loop solution) for the matrix product can be estimated by a runtime of O(lmn), which makes it a computational inefficient approach for large scaled matrices A and B. However, theoretical considerations have been accomplished, in order to reduce the amount of FLOPs. One of these approaches is to use a probabilistic method, that is a random subsampling (according to a certain probability distribution) of the rows in A and according columns in B. By performing the multiplication in terms of the dyadic product it becomes apparent how the supposed dimensions (mxn) of the outcome matrix can be achieved. Note: Even though the size / dimension of the matrix is the same as for the correct produc, the rank is likely to be lower. The question remains how to subsample the rows and columns. Therefore we state a minimization problem in order to reduce the expected value of the difference between supposed matrix product and the approximated matrix product in Frobenius norm. You can find the theoretical consideration in more detailed here:

https://arxiv.org/pdf/1608.04481v1.pdf

Chapter 2: Pages 18 - 22

The chart below draws the outcome of the programm. As expected we have a lower runtime by considering the randomized matrix multiplication. The accompanied error is represented by the average error on each cell, that is the complete error for a particular matrix dimension divided by the number od entries. In order to compute the total error you have to multiply the error by the squared number of dimensions. For instance: Let's assume an average error of 0.1 measured for a matrix dimension of 25x25. The overall error is then computed by 0.1*25*25 = 62.5.

![alt text](https://github.com/NumericalMax/RandomizedMatrixProduct/blob/master/images/all.png)

There are applications, in which an approximate solution of large scaled matrix products is sufficient. The important part might be to find the tradeoff between error size and runtime.

# Dependencies
iostream, iomanip, fstream, random, ctime, map, math.h

# Author
Max Kapsecker, 2018: Remark, that there is always room for optimization. The Code can surely optimized in terms of time -and spaceefficiency. Feel free to report any mistake: max.kapsecker@tum.de
