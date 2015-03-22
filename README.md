# Description
Header only matrix library supporting basic arithmetical operations on matrices as well as more advanced mathematical 
functions such as transposition, inversion, calculation of determinant, row echelon form, qr_decomposition, 
householder matrix calculation, jordan form and finding eigenvalues and eigenvectors. 
Supports complex numbers from standard library.
Interface is similar to interfaces of containers from C++ standard library.

#Usage
To use this library, you have to copy all .h files except test.h to your source directory and include file "matrix.h" 
from your source. Also this library needs C++11, so maybe you will have to tell compiler to use it.
Then you can use class matrix<T> in your code. This class acts like container class, with iterators over rows, 
columns and all elements in matrix. Row and column iterators use special proxy classes matrix<T>::row and matrix<T>::col,
which also acts like containers. They share data with the matrix, so for example writing to row can also change state 
of matrix or some column. You can see examples of usage in file "test.cpp".


