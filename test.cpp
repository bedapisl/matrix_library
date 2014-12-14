//Zapoctovy program pro predmet Programovani v C++ v roce 2013/2014
//Bedrich Pisl

#include "test.h"

#define FAIL(text) {cout << __FILE__ << "::" << __LINE__ << " Test failed in function: " <<  __func__ << "(): " << text << endl; return -1;} 

using namespace std;

bool almost_equal(const matrix<double> & first, const matrix<double> & second)
{
	if((first.row_size() != second.row_size()) || (first.column_size() != second.column_size()))
		return false;

	for(int i=0; i<first.col_size(); ++i)
	{
		for(int j=0; j<first.row_size(); ++j)
		{
			if(abs(first[i][j] - second[i][j]) > 0.0000001)
				return false;
		}
	}
	return true;
}

int test1()
{
	cout.precision(1000);
	matrix<double> a();
	matrix<int> m(3,5);
	matrix<matrix<int>> crazy_matrix(0, 10, m);

	int cislo = 0;
	for(matrix<int>::row_iterator r_it = m.row_begin(); r_it != m.row_end(); ++r_it)
	{
		for(matrix<int>::row::iterator it = r_it->begin(); it != r_it->end(); ++it)
		{
			*it = cislo;
			cislo++;
		}
	}

	int nevim = m[1].at(2);
	if(nevim != 7)
		FAIL("");
	
	if(m[1][2] != 7)
		FAIL("")
	
	m = m.transposition();

	if(m[2][1] != 7)
		FAIL("");

	
	if(m.column_size() != 3)
		FAIL("transpozice nefunguje");
	
	return 0;
}

int test2()
{
	unique_ptr<matrix<double>> m (new matrix<double>(2, 5, 1.5));

	matrix<double>::row_iterator r_it = m->row_begin();
	
	r_it++;
	r_it->front() = 128;

	if(m->column_front().back() != 128)	
		FAIL("")
	
	return 0;
}

void add(int i, int & sum)
{
	sum += i;
}

int test3()
{
	matrix<int> m(3, 3, 0);
	m.at(1, 2) = 3;
	
	int sum = 0;
	auto last_column = m.column_back();
	for(auto it = last_column.begin(); it != last_column.end(); it++)
		sum += *it;
	
	bool exception = false;
	if(sum != 3)
		FAIL("")
	try{
		m.at(13, 1);
	}
	catch (out_of_range ex)
	{
		exception = true;
	}

	if(!exception)
		FAIL("")
	
	return 0;
}

int test4()		//copy constructor a assignement operator
{
	bool vypis = false;

	matrix<int>* jednotkova = new matrix<int>(4,4);
	for(size_t i=0; i<jednotkova->row_size(); i++)
	{
		for(size_t j=0; j<jednotkova->column_size(); j++)
		{
			if(i == j)
				(*jednotkova)[i][j] = 1;
			else
				(*jednotkova)[i][j] = 0;

		}
	}
	if(vypis)
	{
		cout << "jednotkova" << endl;
		cout << *jednotkova;
	}

	matrix<int> jina_matice = *jednotkova;	//copy constructor
	matrix<int> treti_matice(3, 3, 3);

	if(vypis)
	{
		cout << "treti matice" << endl;
		cout << treti_matice;
	}
	jina_matice[1][1] = 33;
	
	treti_matice = *jednotkova;		//assignement operator
	
	if(vypis)
	{
		cout << "treti matice: " << endl;
		cout << treti_matice;
	}

	if(treti_matice[2][2] != 1)
	{
		cout << "treti_matice[2][2] = " << treti_matice[2][2] << " - should be 1" << endl;
		FAIL("")
	}
	
	if(vypis)
	{
		cout << "treti matice: " << endl;
		cout << treti_matice;
	}

	if(vypis)
	{
		cout << "jednotkove matice: " << endl;
		cout << *jednotkova;
	}


	if((*jednotkova) != treti_matice)
		FAIL("")

	delete jednotkova;
	return 0;
}

int test_aritmeticke_operace()		
{
	matrix<int> A(2, 2, 2);
	matrix<int> B(2, 2, 3);
	B.row_back() = A.row_front();

	matrix<int> zero_and_one(2, 2, 0);
	zero_and_one[1][0] = 1;
	zero_and_one[1][1] = 1;

	matrix<int> vysledek = B - A + zero_and_one;
	matrix<int> spravny_vysledek(2, 2, 1);

	if(vysledek != spravny_vysledek)
		FAIL("");
	
	matrix<int> C = A + B + A;
	matrix<int> D = C*zero_and_one;
	A *= B;
	C += A;

	matrix<complex<float>> m1 (1, 1, 1);
	matrix<complex<float>> m2 (1, 1, 2);
	matrix<complex<float>> m3 = m1 + m2;

	matrix<complex<float>> m4;
	m4 = m2 * m2;
	matrix<complex<float>> m5 (1,1,0);
	m5 += m2 += m3;

	return 0;
}

int template_iterator()
{
	matrix<int> m(2,2);

	int value = 0;
	for(matrix<int>::iterator it = m.begin(); it < m.end(); ++it)
	{
		value++;
		*it = value;
	}
	
	auto it = m.begin();
	auto it2 = it + 2;
	
	
	if(*it2 != 3)
		FAIL("");

	it++;
	if((1 + it2) != (m.end() - 1))
		FAIL("");
	
	auto cit = m.column_cend();
	cit--;
	matrix<int>::col sloupec = *cit;
	matrix<int>::col::iterator sloupcovy_iterator = sloupec.begin();
	int a = *sloupcovy_iterator;
	if(a != 2)
		FAIL("");	

	matrix<int>::iterator::pointer pointer_to_int = new int;
	
	return 0;
}

int inversion_test()
{
	
	matrix<double> matice;
	matice = matrix<double>(3,3,5);
	matice[0][0] = 1;
	matice[0][1] = 3;
	matice[0][2] = 5;
	matice[1][0] = matice[1][1] = matice[1][2] = 1;
	matice[2][0] = 1;
	matice[2][1] = 2;
	matice[2][2] = 4;
	
	matrix<double> jednotkova(3, 3, 0);
	for(int i=0; i<3; ++i)
		jednotkova[i][i] = 1;
	
	if(!almost_equal(matice * matice.inversion(), jednotkova))
		FAIL("");
	
	int i=0;
	for(auto it = matice.begin(); it != matice.end(); ++it)
	{
		++i;
		*it = i * i;
	}
	
	matrix<double> inverzni = matice.inversion();	
	matrix<double> transponovana = matice.transposition();
	
	matrix<double> m6(2,2,2);
	if(m6.inversion() != matrix<double>(0, 0))
		FAIL("");

	matrix<complex<double>> m7(4, 4);
	for(auto it = m7.begin(); it != m7.end(); ++it)
	{
		*it = complex<double>(rand() % 20, rand() % 20);
	}
	m7.inversion();

	return 0;
}

int typedef_test()
{
	matrix<int> m(3,3,1);
	auto it = m.begin();
	
	matrix<int>::iterator::pointer p = &(it[0]);
	matrix<int>::value_type hodnota = *p;
	matrix<int>::row::value_type a = 1;

	if(hodnota != a)
		FAIL("");

	return 0;
}

int matice_matic_test()
{
	matrix<int> m(2, 2, 1);
	matrix<matrix<int>> main_matrix(2, 2, m);

	main_matrix[0][0][0][0] = 22;

	return 0;
}
/*		
int checked_test()		//can be ru only if function matrix<T>::the_arithmetics() is public
{
	try
	{
		matrix<int>::the_arithmetics().checked_addition(33, 55);
	}
	catch (exception)
	{
		FAIL("");
	}
	try
	{
		matrix<int>::the_arithmetics().checked_addition(-33, -1255);
	}
	catch (exception)
	{
		FAIL("");
	}
	bool thrown = false;
	try{
		matrix<int>::the_arithmetics().checked_addition(2000000000, 2000000000);
	}
	catch (exception)
	{
		thrown = true;
	}
	if(!thrown)
		FAIL("");
	 
	thrown = false;
	try{
		matrix<int>::the_arithmetics().checked_addition(-2000000000, -2000000000);
	}
	catch (exception)
	{
		thrown = true;
	}
	if(!thrown)
		FAIL("");
	

	try
	{
		matrix<int>::the_arithmetics().checked_multiplication(33, 55);
	}
	catch (exception)
	{
		FAIL("");
	}
	thrown = false;
	try{
		matrix<int>::the_arithmetics().checked_multiplication(20000000, 2000000);
	}
	catch (exception)
	{
		thrown = true;
	}
	if(!thrown)
		FAIL("");
	
	return 0;
}
*/
int determinant_and_ref_test()
{
	matrix<double> m1(2, 2, 1);
	if(m1.determinant() != 0)
		FAIL("");

	matrix<double> m2(3, 3, 0);
	m2[0][0] = -2;
	m2[0][1] = 2;
	m2[0][2] = -3;
	m2[1][0] = -1;
	m2[1][1] = 1;
	m2[1][2] = 3;
	m2[2][0] = 2;
	m2[2][1] = 0;
	m2[2][2] = -1;

	if(m2.determinant() != 18)
		FAIL("");
	
	if(m2.determinant(false) != 18)
		FAIL("");
	
	matrix<double> r = m2.row_echelon_form(false);
	if((r[1][0] != 0) || (r[2][0] != 0) || (r[2][1] != 0))
		FAIL("");

	matrix<double> rectangular_matrix(3, 5, 1);
	rectangular_matrix[2][2] = 0;

	matrix<double> rectangular_matrix2(5, 3, 1);
	rectangular_matrix2[1][2] = 0;

	srand(11);

	matrix<long long> random(4, 7);
	for(auto it = random.begin(); it != random.end(); ++it)
		*it = (long long)(rand() % 10);

	random.row_echelon_form(true);
	
	matrix<double> random2(10, 10);
	for(auto it = random2.begin(); it != random2.end(); ++it)
		*it = (double)(rand() % 10);

	random2.inversion(false);
	random2.row_echelon_form();

	
	matrix<complex<double>> m7(4, 4);
	for(auto it = m7.begin(); it != m7.end(); ++it)
	{
		*it = complex<double>(rand() % 20, rand() % 20);
	}
	
	m7.row_echelon_form();
	m7.determinant();
	
	return 0;
}

int qr_rozklad_test()
{
	matrix<double> m(3, 3);
	m[0][0] = 0;
	m[0][1] = -20;
	m[0][2] = -14;
	m[1][0] = 3;
	m[1][1] = 27;
	m[1][2] = -4;
	m[2][0] = 4;
	m[2][1] = 11;
	m[2][2] = -2;

	auto matrices = m.qr_decomposition();

	matrix<double> q = matrices.first;
	matrix<double> r = matrices.second;
	
	for(int i=0; i<m.row_size(); ++i)
	{
		for(int j=0; j<m.column_size(); ++j)
		{
			if(((q*r)[i][j] + 0.0000001 < m[i][j]) || ((q*r)[i][j] - 0.00000001 > m[i][j]))
			{
				cout << (q*r)[i][j] << " " << m[i][j] << endl;
				cout << q*r << endl << m << endl;
				FAIL("");
			}
		}
	}
	
	matrix<double> m2(2, 3, 4);
	
	pair<matrix<double>, matrix<double>> matrices2 = m2.qr_decomposition();
	matrix<double> q2 = matrices2.first;
	matrix<double> r2 = matrices2.second;
	
	for(int i=0; i<m2.row_size(); ++i)
	{
		for(int j=0; j<m2.column_size(); ++j)
		{
			if(((q2*r2)[i][j] + 0.0000001 < m2[i][j]) || ((q2*r2)[i][j] - 0.00000001 > m2[i][j]))
			{
				cout << (q2*r2)[i][j] << " " << m2[i][j] << endl;
				cout << q2*r2 << endl << m2 << endl;
				FAIL("");
			}
		}
	}
	return 0;
}

int vlastni_cisla()
{
	
	{
		matrix<double> m(2, 2, 0);
		m[0][1] = -1;
		m[1][0] = -1;

		vector<complex<double>> vlastni_cisla =  m.eigenvalues();
		vector<complex<double>> spravna_vlastni_cisla;
		
		spravna_vlastni_cisla.push_back(complex<double>(1, 0));
		spravna_vlastni_cisla.push_back(complex<double>(-1, 0));
		
		if(spravna_vlastni_cisla != vlastni_cisla)
			FAIL("");
	
		vector<vector<complex<double>>> vlastni_vektory = m.eigenvectors(vlastni_cisla);
			//There should be two real eigenvectors (-1, 1) and (1, 1) or there multiple. Program stores them as cmoplex numbers so there would be few more zeros. 
		//print_vector_of_vectors(vlastni_vektory);
	}
	
	{
		matrix<double> m(2, 2, 0);
		m[0][0] = 3;
		m[0][1] = -1;
		m[1][0] = 2;
		m[1][1] = 0;

		vector<complex<double>> vlastni_cisla = m.eigenvalues();

		if(abs(vlastni_cisla[0] - complex<double>(2, 0)) > 0.000001)
			FAIL("");
		
		if(abs(vlastni_cisla[1] - complex<double>(1, 0)) > 0.000001)
			FAIL("");
	
		vector<vector<complex<double>>> vlastni_vektory = m.eigenvectors(vlastni_cisla);
			//vectors should be (1, 1) and (1, 2) or multiples.
	
		//print_vector_of_vectors(vlastni_vektory);
	}
	
	{
		matrix<double> m(2, 2, 0);
		m[0][0] = 3;
		m[0][1] = -2;
		m[1][0] = 4;
		m[1][1] = -1;

		vector<complex<double>> vlastni_cisla = m.eigenvalues();
	
		if(abs(vlastni_cisla[0] - complex<double>(1, 2)) > 0.000001)
			FAIL("");
		
		if(abs(vlastni_cisla[1] - complex<double>(1, -2)) > 0.000001)
			FAIL("");
		
		vector<vector<complex<double>>> vlastni_vektory = m.eigenvectors(vlastni_cisla);
			//vectors should be (1, 1 + i) and (1, 1 - i) or multiples.

		matrix<complex<double>> special_matrix = m.to_complex();
		
		if(special_matrix.eigenvalues() != vlastni_cisla)
			FAIL("");

		special_matrix.eigenvectors(vlastni_cisla);

	}
	return 0;
}

int definitnost_a_jordanova_forma()
{
	matrix<double> m(3, 3, -1);
	m[0][0] = 2;
	m[0][2] = 0;
	m[2][0] = 0;
	m[1][1] = 2;
	m[2][2] = 2;

	if(!m.positive_definite())
		FAIL("");

	if(!matrix<double>::identity_matrix(5).positive_definite())
		FAIL("");

	m[1][1] = -100;
	if(m.positive_definite())
		FAIL("");

	if(m.positive_semi_definite())
		FAIL("");
	
	if(!matrix<double>(2, 2, 0).positive_semi_definite())
		FAIL("");

	matrix<double> m2(3, 3, 0);
	m2[0][0] = 3;
	m2[0][1] = 1;
	m2[1][0] = -1;
	m2[1][1] = 1;
	m2[2][0] = 3;
	m2[2][1] = 2;
	m2[2][2] = 2;

	vector<complex<double>> vlastni_cisla;
	for(int i=0; i<3; ++i)
		vlastni_cisla.push_back(complex<double>(2, 0));

	matrix<double> spravny_jordan(3, 3, 0);
	spravny_jordan[0][0] = 2;
	spravny_jordan[1][1] = 2;
	spravny_jordan[2][2] = 2;
	spravny_jordan[0][1] = 1;
	spravny_jordan[1][2] = 1;

	if(spravny_jordan.to_complex() != m2.jordan_form(vlastni_cisla))
		FAIL("");

	matrix<double> m3(3, 3, 0);
	m3[0][0] = 2;
	m3[0][1] = 4;
	m3[0][2] = -8;
	m3[1][2] = 4;
	m3[2][1] = -1;
	m3[2][2] = 4;

	spravny_jordan[0][1] = 0;

	if(spravny_jordan.to_complex() != m3.jordan_form(vlastni_cisla))
		FAIL("");
	
	matrix<complex<double>> special = m3.to_complex();
	if(special.jordan_form(vlastni_cisla) != spravny_jordan.to_complex())
		FAIL("");
	

	return 0;
}






