//Bedrich Pisl - Programming in C++, MFF, 2013/2014

#include "test.h"

using namespace std;

int main()
{
	if(test1() != 0)
		return -1;

	if(test2() != 0)
		return -1;

	if(test3() != 0)
		return -1;

	if(test4() != 0)
		return -1;

	if(test_arithmetical_operations() != 0)
		return -1;
	
	if(template_iterator() != 0)
		return -1;
	
	if(inversion_test() != 0)
		return -1;
	
	if(typedef_test() != 0)
		return -1;
	
	if(matrix_of_matrices_test() != 0)
		return -1;
	
	//if(checked_test() != 0)
	//	return -1;
	
	if(determinant_and_ref_test() != 0)
		return -1;
	
	if(qr_decomposition_test() != 0)
		return -1;
	
	if(eigenvalues_test() != 0)
		return -1;
		
	if(definite_and_jordan_form() != 0)
		return -1;

	return 0;
}
