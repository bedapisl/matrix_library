//Bedrich Pisl - Programming in C++, MFF, 2013/2014

#ifndef TEST_H
#define TEST_H

#include "matrix.h"
#include <iostream>
#include <stdexcept>
#include <memory>
#include <fstream>

int test1();

int test2();

int test3();

int test4();

int test_arithmetical_operations();

int template_iterator();

int inversion_test(); 

int typedef_test();

int matrix_of_matrices_test();

int checked_test();

int determinant_and_ref_test();

int qr_decomposition_test();

int eigenvalues_test();

int definite_and_jordan_form();

bool almost_equal(const matrix<double> & first, const matrix<double> & second);

template <typename T>
void print_vector_of_vectors(const std::vector<std::vector<T>> & v)
{
	std::cout << "printing vector of vectors" << std::endl;
	for(int i=0; i<v.size(); ++i)
	{
		for(int j=0; j<v[i].size(); ++j)
		{
			std::cout << v[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

#endif
