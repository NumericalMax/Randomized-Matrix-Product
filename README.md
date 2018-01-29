# Randomized Matrix Product
Let us consider the matrix product of A \in \mathbb{R}{l \times m} and B \in \mathbb{R}{m \times n}.
It is well known, that the classical approach for matrix products can estimated by a runtime of O(lmn), which makes it a computational inefficient approach for large scaled matrices A and B. However, theoretical considerations have been accomplished, in order to reduce the amount of FLOPs. One of these approaches is to use a probabilistic method for the matrix product. That is a random subsampling (followed by a certain probability distribution) of the rows in A and according columns in B. By representing the result in terms of outer products one can ensure the supposed dimensions of the outcome matrix, but with a lower rank. Finally we like to minimize the distance of the exact answer AB and the approximate answer A'B'. Hence, the probability distribution for the process of subsampling has to be chose such that the norm ||AB - A'B'||_{F} is minimized. You can find the theoretical consideration on this more detailed in:

https://arxiv.org/pdf/1608.04481v1.pdf

Chapter 2: Pages 18 - 22
