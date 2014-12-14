//Zapoctovy program pro predmet Programovani v C++ v roce 2013/2014
//Bedrich Pisl

#ifndef BEDAS_HEADER_H
#define BEDAS_HEADER_H

template <typename T>
class matrix_general;

template <typename T>
class matrix;

template <typename T, bool built_in_type>
class arithmetics;

template <typename T>
void swap(const matrix_general<T>& first, const matrix_general<T>& second);

template <typename T>
bool operator== (const matrix_general<T>& first, const matrix_general<T>& second);

template <typename T>
bool operator!= (const matrix_general<T>& first, const matrix_general<T>& second);

template <typename T>
matrix<T> operator+ (const matrix_general<T>& first, const matrix_general<T>& second);

template <typename T>
matrix<T> operator- (const matrix_general<T>& first, const matrix_general<T>& second);

template <typename T>
matrix<T> operator* (const matrix_general<T>& first, const matrix_general<T>& second);

template <typename T>
matrix<T> operator* (const T& multiplicator, const matrix_general<T>& m);

template <typename T>
matrix<T> operator* (const matrix_general<T>& m, const T& multiplicator);

template <typename T>
std::ostream& operator<< (std::ostream &str, const matrix_general<T> &m);

template <typename T>
class matrix_general
{
public:
	class iterator;
	class row_iterator;
	class col_iterator;
	class const_iterator;
	class row_const_iterator;
	class col_const_iterator;
	class row;
	class col;
	class matrix_base;
	
	typedef T value_type;
	typedef T& reference;
	typedef const T& const_reference;
	typedef std::ptrdiff_t difference_type;
	typedef std::size_t size_type;
	typedef row row_value_type;
	typedef col column_value_type;
	typedef row& row_reference;
	typedef col& column_reference;
	typedef const row& row_const_reference;
	typedef const col& column_const_reference;
	typedef col_iterator column_iterator;
	typedef col_const_iterator column_const_iterator;
	typedef std::ptrdiff_t row_difference_type;
	typedef std::ptrdiff_t column_difference_type;
	typedef std::size_t row_size_type;
	typedef std::size_t column_size_type;
	
	matrix_general();
	matrix_general(std::size_t rows, std::size_t cols); 
	matrix_general(std::size_t rows, std::size_t cols, const T& initial_value); 
	matrix_general(const matrix_general<T>& other) : base(new matrix_base(other.base->deep_copy())) {}	//copy constructor
	matrix_general(matrix_general<T>&& other) : base(new matrix_base(*other.base)) {*(other.base) = matrix_base(NULL, 0, 0);}		//move constructor
	~matrix_general() {}
	const matrix_general<T>& operator= (const matrix_general<T>& other);
	const matrix_general<T>& operator= (matrix_general<T>&& other);

	friend bool operator== <> (const matrix_general<T>& first, const matrix_general<T>& second); 
	friend bool operator!= <> (const matrix_general<T>& first, const matrix_general<T>& second);
	friend matrix<T> operator+ <> (const matrix_general<T>& first, const matrix_general<T>& second);
	friend matrix<T> operator- <> (const matrix_general<T>& first, const matrix_general<T>& second);
	friend matrix<T> operator* <> (const matrix_general<T>& first, const matrix_general<T>& second);
	friend matrix<T> operator* <> (const T& multiplicator, const matrix_general<T>& m);
	friend matrix<T> operator* <> (const matrix_general<T>& m, const T& multiplicator);
	
	matrix<T>& operator+= (const matrix_general<T>& other);
	matrix<T>& operator-= (const matrix_general<T>& other);
	matrix<T>& operator*= (const matrix_general<T>& other);
	
	friend void swap <> (const matrix_general<T>& first, const matrix_general<T>& second);
	void swap(matrix_general<T>& other);

	friend std::ostream& operator<< <> (std::ostream & str, const matrix_general<T> &m);	
	row operator[] (std::size_t index) {return row(base, index*base->row_length);}
	const row operator[] (std::size_t index) const {return row(base, index * base->row_length);}
	T& at(std::size_t row, std::size_t col);
	const T& at(std::size_t row, std::size_t col) const;
	row row_at(std::size_t index) {return row(base, index * base->row_length);}
	const row row_at(std::size_t index) const {return row(base, index * base->row_length);}
	col column_at(std::size_t index) {return col(base, index);}
	const col column_at(std::size_t index) const {return col(base, index);}
	
	iterator begin() {return iterator(base, 0);}
	iterator end() {return iterator(base, base->row_length * base->col_length);}
	const_iterator begin() const {return const_iterator(base, 0);} 
	const_iterator end() const {return const_iterator(base, base->row_length * base->col_length);}
	const_iterator cbegin() const {return const_iterator(base, 0);} 
	const_iterator cend() const {return const_iterator(base, base->row_length * base->col_length);}
	
	T& front() {return base->data[0];}
	T& back() {return base->data[(base->row_length * base->col_length) - 1];}
	const T& front() const {return front();}
	const T& back() const {return back();}

	row_iterator row_begin() {return row_iterator(base, 0);}
	row_iterator row_end() {return row_iterator(base, base->row_length * base->col_length);}	
	row_const_iterator row_begin() const {return row_const_iterator(base, 0);}
	row_const_iterator row_end() const {return row_const_iterator(base, base->row_length * base->col_length);}
	row_const_iterator row_cbegin() const {return row_const_iterator(base, 0);}
	row_const_iterator row_cend() const {return row_const_iterator(base, base->row_length * base->col_length);}
	
	row row_front() {return row(base, 0);}
	row row_back() {return row(base, base->row_length * (base->col_length - 1));}
	const row row_front() const {return row_front();}
	const row row_back() const {return row_back();}

	col_iterator column_begin() {return col_iterator(base, 0);}
	col_iterator column_end() {return col_iterator(base, base->row_length);}
	col_const_iterator column_begin() const {return col_const_iterator(base, 0);}
	col_const_iterator column_end() const {return col_const_iterator(base, base->row_length);}
	col_const_iterator column_cbegin() const {return col_const_iterator(base, 0);}
	col_const_iterator column_cend() const {return col_const_iterator(base, base->row_length);}
	
	col column_front() {return col(base, 0);}
	col column_back() {return col(base, base->row_length - 1);}
	const col column_front() const {return column_front();}
	const col column_back() const {return column_back();}

	std::size_t size() const {return base->row_length * base->col_length;}
	std::size_t row_size() const {return base->col_length;}
	std::size_t column_size() const {return base->row_length;}
	std::size_t col_size() const {return base->row_length;}
	
	bool empty() const {return size() == 0;}
	bool row_empty() const {return row_size() == 0;}
	bool column_empty() const {return column_size() == 0;}

	matrix<T> transposition() const;
	matrix<T> inversion(bool check_overflow = true) const;	//Overflow checking works only for built-in types. If matrix is not regular, function returns matrix 0 x 0.
	T determinant(bool check_overflow = true) const;	
	matrix<T> row_echelon_form(bool check_overflow = true) const;
	std::pair<matrix<T>, matrix<T>>  qr_decomposition() const;	//doesnt work for integers (because it uses numerical method), first matrix is ortogonal, second is triangular
				
	bool positive_definite() const;
	bool positive_semi_definite() const;		

	std::size_t rank() const;
	
	static matrix<T> identity_matrix(std::size_t row_size);
	template <typename Iterator>
	static matrix<T> householder_matrix(const Iterator first, const Iterator last);
	template <typename Iterator1, typename Iterator2>
	static matrix<T> outer_product(const Iterator1 first1, const Iterator1 last1, const Iterator2 first2, const Iterator2 last2);
	template <typename Iterator1, typename Iterator2>
	static T dot_product(const Iterator1 first1, const Iterator1 last1, const Iterator2 first2, const Iterator2 last2);
	matrix<std::complex<T>> to_complex() const;

protected:
	void nullify_first_element(std::size_t row_to_nullify, std::size_t row_to_add, std::size_t first_non_zero_element);
	void checked_nullify_first_element(std::size_t row_to_nullify, std::size_t row_to_add, std::size_t first_non_zero_element);
	void add_row(matrix_general<T>::row row_to_add, std::size_t row_index);
	void checked_add_row(matrix_general<T>::row row_to_add, std::size_t row_index);
	void multiply_row(T multiplicator, std::size_t row_index);
	void checked_multiply_row(T multiplicator, std::size_t row_index);
	void add_row_multiplication(T multiplicator, std::size_t row_to_add, std::size_t destination_row);	
	void swap_rows_elements(std::size_t first_row_index, std::size_t second_row_index);

	static arithmetics<T, std::numeric_limits<T>::is_specialized> & the_arithmetics();	//defined in other.h

	//following 5 functions will be called only if T is complex<U> for some type U
	std::vector<T> eigenvalues_impl(unsigned int iterations) const;	
	std::vector<std::vector<T>> eigenvectors_impl(std::vector<T> eigevalues) const;	
	std::vector<T> find_eigenvector(const matrix_general<T> & m, std::size_t nonzero_variable, std::size_t first_substitution) const;
	matrix<T> jordan_form_impl(std::vector<T> eigenvalues) const; 
	std::vector<std::size_t> numbers_of_jordans_cells(T eigenvalue, std::size_t multiplicity) const;

	matrix<T>* to_matrix_ptr() { return static_cast<matrix<T>*>(this);}
	const matrix<T>* to_matrix_ptr() const { return static_cast<const matrix<T>*>(this);}

	std::shared_ptr<matrix_base> base;
};

template <typename T>
class matrix : public matrix_general<T>		//for all types except complex<X> 
{
public:
	matrix() : matrix_general<T>() {}
	matrix(std::size_t rows, std::size_t cols) : matrix_general<T>(rows, cols) {}
	matrix(std::size_t rows, std::size_t cols, const T& initial_value) : matrix_general<T>(rows, cols, initial_value) {}
	matrix(const matrix<T>& other) : matrix_general<T>(other) {}
	matrix(const matrix_general<T>& other) : matrix_general<T>(other) {}
	matrix(matrix_general<T>&& other) : matrix_general<T>(other) {}		//move constructor
	matrix(matrix<T>&& other) : matrix_general<T>(other) {}
	
	const matrix<T>& operator= (matrix<T>&& other) { matrix_general<T>::operator= (other); return *this;}
	const matrix<T>& operator= (const matrix<T>& other) { matrix_general<T>::operator= (other); return *this;}

	using matrix_general<T>::the_arithmetics;
	using matrix_general<T>::row_size;
	using matrix_general<T>::col_size;

	std::vector<std::complex<T>> eigenvalues(unsigned int iterations = 100) const {return this->to_complex().eigenvalues(iterations);} //calls partial specialization for complex numbers	

	std::vector<std::vector<std::complex<T>>> eigenvectors(std::vector<std::complex<T>> eigenvalues) const {return this->to_complex().eigenvectors(eigenvalues);}
	
	matrix<std::complex<T>> jordan_form(std::vector<std::complex<T>> eigenvalues) const {return this->to_complex().jordan_form(eigenvalues);}

};

template <typename T>
class matrix<std::complex<T>> : public matrix_general<std::complex<T>> //partial specialization for complex numbers
{
public:
	matrix() : matrix_general<std::complex<T>>() {}
	matrix(std::size_t rows, std::size_t cols) : matrix_general<std::complex<T>>(rows, cols) {}
	matrix(std::size_t rows, std::size_t cols, const T& initial_value) : matrix_general<std::complex<T>>(rows, cols, initial_value) {}
	matrix(const matrix<std::complex<T>>& other) : matrix_general<std::complex<T>>(other) {}
	matrix(const matrix_general<std::complex<T>>& other) : matrix_general<std::complex<T>>(other) {}
	matrix(matrix_general<std::complex<T>>&& other) : matrix_general<std::complex<T>>(other) {}		//move constructor
	matrix(matrix<std::complex<T>>&& other) : matrix_general<std::complex<T>>(other) {}

	const matrix<std::complex<T>>& operator= (matrix<std::complex<T>>&& other) { matrix_general<std::complex<T>>::operator= (other); return *this;}
	const matrix<std::complex<T>>& operator= (const matrix<std::complex<T>>& other) { matrix_general<std::complex<T>>::operator= (other); return *this;}

	using matrix_general<std::complex<T>>::the_arithmetics;
	using matrix_general<std::complex<T>>::row_size;
	using matrix_general<std::complex<T>>::col_size;

	std::vector<std::complex<T>> eigenvalues(unsigned int iterations = 100) const {return this->eigenvalues_impl(iterations);}	

	std::vector<std::vector<std::complex<T>>> eigenvectors(std::vector<std::complex<T>> eigenvalues) const {return this->eigenvectors_impl(eigenvalues);}
	
	matrix<std::complex<T>> jordan_form(std::vector<std::complex<T>> eigenvalues) const {return this->jordan_form_impl(eigenvalues);}
};
	

#endif
