//Zapoctovy program pro predmet Programovani v C++ v roce 2013/2014
//Bedrich Pisl

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

	if(test_aritmeticke_operace() != 0)
		return -1;
	
	if(template_iterator() != 0)
		return -1;
	
	if(inversion_test() != 0)
		return -1;
	
	if(typedef_test() != 0)
		return -1;
	
	if(matice_matic_test() != 0)
		return -1;
	
	//if(checked_test() != 0)
	//	return -1;
	
	if(determinant_and_ref_test() != 0)
		return -1;
	
	if(qr_rozklad_test() != 0)
		return -1;
	
	if(vlastni_cisla() != 0)
		return -1;
		
	if(definitnost_a_jordanova_forma() != 0)
		return -1;

	return 0;
}
