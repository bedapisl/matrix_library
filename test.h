//Zapoctovy program pro predmet Programovani v C++ v roce 2013/2014
//Bedrich Pisl

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

int test_aritmeticke_operace();

int template_iterator();

int inversion_test(); 

int typedef_test();

int matice_matic_test();

int checked_test();

int determinant_and_ref_test();

int qr_rozklad_test();

int vlastni_cisla();

int definitnost_a_jordanova_forma();

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
