//Bedrich Pisl - Programming in C++, MFF, 2013/2014

#ifndef BEDAS_IMPL_H
#define BEDAS_IMPL_H

template <typename T>
matrix_general<T>::matrix_general() : base(new matrix_base(NULL, 0, 0))
{
}

template <typename T>
matrix_general<T>::matrix_general(std::size_t rows, std::size_t cols) : base(new matrix_base(new T [rows* cols], cols, rows))
{
}

template <typename T>
matrix_general<T>::matrix_general(std::size_t rows, std::size_t cols, const T& initial_value) : base(new matrix_base(new T [rows* cols], cols, rows))
{	
	for(std::size_t i=0; i<rows * cols; ++i)
		base->data[i] = initial_value;
}

template <typename T>
T& matrix_general<T>::at(std::size_t row, std::size_t col)
{
	if((row >= base->col_length) || (col >= base->row_length))
	{
		throw std::out_of_range("matrix::at - range check");
	}
	return *(base->data + base->row_length*row + col);
}

template <typename T>
const T& matrix_general<T>::at(std::size_t row, std::size_t col) const
{
	if((row >= base->col_length) || (col >= base->row_length))
	{
		throw std::out_of_range("matrix::at - range check");
	}
	return *(base->data + base->row_length*row + col);
}

template <typename T>
const matrix_general<T>& matrix_general<T>::operator= (const matrix_general<T>& other) 
{
	if(base->row_length * base->col_length != other.base->row_length * other.base->col_length)
	{
		base = std::shared_ptr<matrix_base>(new matrix_base(other.base->deep_copy()));
	}
	else
	{
		base->row_length = other.base->row_length;
		base->col_length = other.base->col_length;
		for(std::size_t i=0; i<base->row_length * base->col_length; ++i)
			base->data[i] = other.base->data[i];
	}
	return *this;
}

template <typename T>
const matrix_general<T>& matrix_general<T>::operator= (matrix_general<T>&& other)
{
	if(this == &other)
		return *this;
	
	auto tmp = base;	//swap
	base = other.base;
	other.base = tmp;

	*other.base = matrix_base(NULL, 0, 0);

	return *this;
}

template <typename T>
bool operator== (const matrix_general<T>& first, const matrix_general<T>& second) 
{
	return *(first.base) == *(second.base);
}

template <typename T>
bool operator!= (const matrix_general<T>& first, const matrix_general<T>& second) 
{
	return !(*(first.base) == *(second.base));
}

template <typename T>
void matrix_general<T>::swap(matrix_general<T>& other) 
{
	matrix_base tmp = base;
	base = other.base;
	other.base = tmp;
}

template <typename T>
void swap(const matrix_general<T>& first, const matrix_general<T>& second)
{
	typename matrix<T>::matrix_base tmp = first.base;
	first.base = second.base;
	second.base = tmp;
}

template <typename T>
matrix<T> operator+ (const matrix_general<T>& first, const matrix_general<T>& second)
{
	if((first.base->row_length != second.base->row_length) || (first.base->col_length != second.base->col_length))
		throw matrix_exception("Cannot add matrices of different sizes.");

	matrix_general<T> result(first.base->row_length, first.base->col_length);
	for(std::size_t i=0; i<first.base->row_length * first.base->col_length; ++i)
		result.base->data[i] = first.base->data[i] + second.base->data[i];
	
	return result;
}

template <typename T>
matrix<T> operator- (const matrix_general<T>& first, const matrix_general<T>& second)
{
	if((first.base->row_length != second.base->row_length) || (first.base->col_length != second.base->col_length))
		throw matrix_exception("Cannot subtract matrices of different sizes.");

	matrix_general<T> result(first.base->row_length, first.base->col_length);
	for(std::size_t i=0; i<first.base->row_length * first.base->col_length; ++i)
		result.base->data[i] = first.base->data[i] - second.base->data[i];
	
	return result;
}

template <typename T>
matrix<T> operator* (const matrix_general<T>& first, const matrix_general<T>& second)
{
	if(first.base->row_length != second.base->col_length)
		throw matrix_exception("Invalid size of matrices for multiplication");

	matrix_general<T> result(first.base->col_length, second.base->row_length);
	for(std::size_t i=0; i<result.base->col_length; ++i)
	{	
		for(std::size_t j=0; j<second.base->row_length; ++j)
		{
			T sum = first.the_arithmetics().zero;
			for(std::size_t k=0; k<first.base->col_length; k++)
				sum += (first.base->data[i * first.base->row_length + k]) * (second.base->data[k * second.base->row_length + j]);
			
			result.base->data[i*result.base->row_length + j] = sum;
		}
	}
	return result;
}

template <typename T>
matrix<T> operator* (const T& multiplicator, const matrix_general<T>& m)
{
	matrix_general<T> return_matrix = m;
	for(auto it = return_matrix.begin(); it != return_matrix.end(); ++it)
	{
		(*it) = (*it) * multiplicator;
	}
	return return_matrix;
}

template <typename T>
matrix<T> operator* (const matrix_general<T>& m, const T& multiplicator)
{
	return multiplicator * m;
}

template <typename T>
matrix<T>& matrix_general<T>::operator+= (const matrix_general<T>& other)
{
	if((base->row_length != other.base->row_length) || (base->col_length != other.base->col_length))
		throw matrix_exception("Cannot add matrices of different sizes.");

	for(std::size_t i=0; i<base->row_length * base->col_length; ++i)
		base->data[i] += other.base->data[i];
	
	return *this->to_matrix_ptr();
}

template <typename T>
matrix<T>& matrix_general<T>::operator-= (const matrix_general<T>& other)
{
	if((base->row_length != other.base->row_length) || (base->col_length != other.base->col_length))
		throw matrix_exception("Cannot subtract matrices of different sizes.");

	for(std::size_t i=0; i<base->row_length * base->col_length; ++i)
		base->data[i] -= other.base->data[i];
	
	return *this->to_matrix_ptr();
}

template <typename T>			
matrix<T>& matrix_general<T>::operator*= (const matrix_general<T>& other)
{
	*this = (*this) * (other);
	return *this->to_matrix_ptr();
}

template <typename T>
std::ostream& operator<< (std::ostream & str, const matrix_general<T> & m)
{
	for(std::size_t i=0; i<m.base->col_length; ++i)
	{
		for(std::size_t j=0; j<m.base->row_length; ++j)
			str << m.base->data[i * m.base->row_length + j] << "\t";
		
		str << std::endl;
	}
	return str;
}

template <typename T>
matrix<T> matrix_general<T>::transposition() const
{
	const matrix<T> &m = *this;
	matrix_general<T> transposed(m.column_size(), m.row_size());	
	
	for(std::size_t i=0; i<m.row_size(); ++i)
	{
		for(std::size_t j=0; j<m.column_size(); ++j)
			transposed[j][i] = m[i][j];
	}
	return transposed;
}

//if *this is regular matrix returns inversion, otherwise returns matrix 0x0, 
//Using Gauss - Jordan elimination to convert matrix m to identity matrix and performing same row operations
//to identity matrix to convert it to inverse matrix
template <typename T>
matrix<T> matrix_general<T>::inversion(bool overflow_check) const
{
	matrix_general<T> m = *this;
	T one = the_arithmetics().one;
	T zero = the_arithmetics().zero;
	matrix_general<T> inverted = identity_matrix(m.row_size());
	
	bool check = overflow_check && std::numeric_limits<T>::is_specialized;

	if(m.row_size() != m.column_size())
		return matrix_general<T>(0, 0);

	for(std::size_t i=0; i<inverted.row_size(); ++i)	//Inverted becomes identity matrix	
		inverted[i][i] = one;
	
	for(std::size_t i=0; i<m.row_size(); ++i)
	{
		
		bool only_zeros = true;
		std::size_t max_index = the_arithmetics().find_nonzero_row(m, i, i, only_zeros);
		
		if(only_zeros)
			return matrix_general<T>(0, 0);		//Singular matrix
	
		m.swap_rows_elements(i, max_index);
		inverted.swap_rows_elements(i, max_index);

		T multiplicator = one / m[i][i];
		if(check)
		{	
			m.checked_multiply_row(multiplicator, i);
			inverted.checked_multiply_row(multiplicator, i);
		}
		else
		{
			m.multiply_row(multiplicator, i);
			inverted.multiply_row(multiplicator, i);
		}
	
		for(std::size_t j=0; j<m.row_size(); ++j)
		{
			if((m[j][i] != zero) && (i != j))
			{
				T divisor = - m[j][i];
				if(check)
				{
					m.checked_multiply_row(one / divisor, j);
					inverted.checked_multiply_row(one / divisor, j);
					m.checked_add_row(m[i], j);
					inverted.checked_add_row(inverted[i], j);
					m.checked_multiply_row(divisor, j);
					inverted.checked_multiply_row(divisor, j);	
				}
				else
				{
					m.multiply_row(one / divisor, j);
					inverted.multiply_row(one / divisor, j);
					m.add_row(m[i], j);
					inverted.add_row(inverted[i], j);
					m.multiply_row(divisor, j);
					inverted.multiply_row(divisor, j);	
				}
			}
		}
	}
	return inverted;
}

/*Matrix is converted to identity matrix by Gauss - Jordan elimination. Determinant of identity matrix is one
and because of determinant and elementary row operation characteristics we can compute determinant.
If matrix cannot be converted to identity matrix, then it is singular and determinant is zero.*/

template <typename T>
T matrix_general<T>::determinant(bool check_overflow) const 
{
	matrix_general<T> m = *this;
	T one = the_arithmetics().one;
	T zero = the_arithmetics().zero;

	bool check = check_overflow && std::numeric_limits<T>::is_specialized;

	if(m.row_size() != m.column_size())
		return zero;

	T determinant_divisor = one;
	T determinant = one;

	for(std::size_t i=0; i<m.row_size(); ++i)
	{
		bool only_zeros = true;
		std::size_t max_index = the_arithmetics().find_nonzero_row(m, i, i, only_zeros);
	
		if(only_zeros)		//singular matrix has determinant equal zero
			return zero;
		
		if(i != max_index) 			//change of two rows
			determinant_divisor = - determinant_divisor;

		m.swap_rows_elements(i, max_index);
		
		for(std::size_t j = i + 1; j<m.row_size(); ++j)
		{
			if(m[j][i] != zero)		
			{
				T multiplicator = m[j][i];
				
				if(check)
				{
					m.checked_multiply_row(m[i][i], j);
					for(std::size_t k=i; k<m.row_size(); ++k)
					{
						T subtrahend = the_arithmetics().checked_multiplication(m[i][k], multiplicator);
						m[j][k] = the_arithmetics().checked_subtraction(m[j][k], subtrahend);
					}
					determinant_divisor *= m[i][i];
				}
				else
				{
					m.multiply_row(m[i][i], j);
					for(std::size_t k=i; k<m.row_size(); ++k)
						m[j][k] -= m[i][k] * multiplicator;

					determinant_divisor *= m[i][i];
				}
			}
		}
	}
	for(std::size_t i=0; i<m.row_size(); ++i)
	{
		if(check)
			determinant = the_arithmetics().checked_multiplication(determinant, m[i][i]);

		else
			determinant *= m[i][i];
	}

	determinant = determinant / determinant_divisor;

	return determinant;
}

template <typename T>
matrix<T> matrix_general<T>::row_echelon_form(bool check_overflow) const
{
	matrix_general<T> ref = *this;
	T zero = the_arithmetics().zero;
	
	bool check = check_overflow && std::numeric_limits<T>::is_specialized;

	int starting_row = 0;
	for(std::size_t col_index = 0; col_index < ref.column_size(); ++col_index)	
	{
		for(std::size_t row_index = starting_row; row_index<ref.row_size(); ++row_index)
		{
			if(ref[row_index][col_index] != zero)
			{	
				ref.swap_rows_elements(starting_row, row_index);

				for(std::size_t i=starting_row + 1; i<ref.row_size(); ++i)
				{
					if(ref[i][col_index] != zero)
					{
						if(check)
							ref.checked_nullify_first_element(i, row_index, col_index);
						
						else
							ref.nullify_first_element(i, row_index, col_index);
					}
				}
				starting_row++;
				break;
			}
		}
	}
	return ref;
}

//using householder transformation
//will not work with integers, because square root is needed during computation
//first return value is orthogonal matrix, second is triangular
template <typename T>
std::pair<matrix<T>, matrix<T>> matrix_general<T>::qr_decomposition() const
{
	matrix_general<T> triangular = *this;
	matrix_general<T> orthogonal = identity_matrix(triangular.row_size());

	std::size_t steps = (triangular.row_size() > triangular.column_size()) ? triangular.row_size() : triangular.column_size();

	T two = 2;

	for(std::size_t i = 0; i<steps; ++i)
	{
		typename col::iterator vector_first = (*(triangular.column_begin() + i)).begin();
		typename col::iterator vector_middle = (*(triangular.column_begin() + i)).begin() + i;
		typename col::iterator vector_last = (*(triangular.column_begin() + i)).end();

		T vector_norm = the_arithmetics().square_root(dot_product<typename col::iterator, typename col::iterator>(vector_middle, vector_last, vector_middle, vector_last));
		
		std::vector<T> v(vector_first, vector_last);
		
		for(int i=0; i<vector_middle - vector_first; ++i)
			v[i] = the_arithmetics().zero;

		v[i] -= vector_norm;

		matrix_general<T> identity = identity_matrix(triangular.row_size());
		
		matrix_general<T> outer_product_matrix = outer_product<typename std::vector<T>::iterator, typename std::vector<T>::iterator>(v.begin(), v.end(), v.begin(), v.end());
	
		T norm = dot_product<typename std::vector<T>::iterator, typename std::vector<T>::iterator>(v.begin(), v.end(), v.begin(), v.end());
		
		if(the_arithmetics().equal_zero(norm))
			continue;
		
		for(matrix_general<T>::iterator it = outer_product_matrix.begin(); it != outer_product_matrix.end(); ++it)
			*it = (*it) * two / norm;
	
		matrix_general<T> householder = identity - outer_product_matrix;
	
		triangular = householder * triangular;
		orthogonal = orthogonal * householder;
	}

	return std::pair<matrix_general<T>, matrix_general<T>>(orthogonal, triangular);
}

//computes eigenvalues with qr algorithm
//it is numerical method so the return values are only approximations
template <typename T>		//T is complex<U> for some U
std::vector<T> matrix_general<T>::eigenvalues_impl(unsigned int iterations) const
{
	if(row_size() != column_size())	//only square matrix has eigenvalues
		return std::vector<T>(0);
	
	matrix_general<T> m = *this;
	matrix_general<T> triangular;
	matrix_general<T> orthogonal;

	//T zero = matrix_general<T>::the_arithmetics().zero;
	//std::complex<T> complex_zero(zero, zero);

	T four = 4;
	
	for(unsigned int i = 0; i<iterations; ++i)
	{
		std::pair<matrix_general<T>, matrix_general<T>> decomposed = m.qr_decomposition();
		orthogonal = decomposed.first;
		triangular = decomposed.second;
		m = triangular * orthogonal;
	}
	
	std::vector<T> eigenvalues;

	for(std::size_t i=0; i<matrix_general<T>::row_size(); ++i)
	{
		if((i == matrix_general<T>::row_size() - 1) || (m[i][i+1] == the_arithmetics().zero))
			eigenvalues.push_back(m[i][i]);			//everything at i column except m[i][i] should be aproximately 0, so m[i][i] is eigenvalue
		
		else
		{
					//eigenvalues will be computed as roots of characteristic polynom of 2x2 submatrix
			
			T determinant = pow(m[i][i] - m[i+1][i+1], 2) + four * (m[i][i+1]) * (m[i+1][i]);
			T eigenvalue1 = (m[i][i] + m[i+1][i+1] + sqrt(determinant)) / 2.0;		//formula for quadratic equation
			T eigenvalue2 = (m[i][i] + m[i+1][i+1] - sqrt(determinant)) / 2.0;

			eigenvalues.push_back(eigenvalue1);
			eigenvalues.push_back(eigenvalue2);
			
			i++;
		}
	}

	return eigenvalues;
}
//computes eigenvectors by finding all nonzero solution of (A - eI)x = 0 for all eigenvalues e 
template <typename T>		//T is complex<U>
std::vector<std::vector<T>> matrix_general<T>::eigenvectors_impl(std::vector<T> eigenvalues) const
{
	if(row_size() != column_size())	//only square matrices has eigenvalues
		return std::vector<std::vector<T>>(0);
	
	std::sort(eigenvalues.begin(), eigenvalues.end(), [](T first, T second){ 		//lambda to compare complex numbers
		if(first.real() == second.real())
			return first.imag() < second.imag();
		return first.real() < second.real();
		});

	std::vector<std::vector<T>> eigenvectors;
	
	matrix_general<T> complex_matrix = *this;//to_complex();
	matrix_general<T> complex_identity = identity_matrix(row_size());

	for(std::size_t i=0; i<eigenvalues.size(); ++i)
	{
		if((i != 0) && (eigenvalues[i-1] == eigenvalues[i]))		//same eigenvalues has same eigenvectors
			continue;

		matrix_general<T> reduced = (complex_matrix - eigenvalues[i] * complex_identity).row_echelon_form();
	
		std::size_t last_nonzero_row = reduced.row_size() - 1;

		while(the_arithmetics().equal_zero(reduced[last_nonzero_row][reduced.column_size() - 1]))
			last_nonzero_row--;
		last_nonzero_row++;

		std::size_t col_index = reduced.col_size();
		std::size_t row_index = last_nonzero_row;
		
		while(row_index != 0)			//solve equation reduced * x = 0 with substitution method
		{
			--row_index; --col_index;
			reduced.multiply_row( - (the_arithmetics().one / reduced[row_index][col_index]), row_index);	//convert x_n where n = "col_index" in equation on row "row_index" on the other side of equation
			reduced[row_index][col_index] = 0;
			
			for(std::size_t j=0; j < row_index; j--)		
			{							//substitution for x_n where n = col_index 
				reduced.add_row_multiplication(reduced[j][col_index], row_index, j);
				reduced[j][col_index] = 0;
			}
		}
		
		std::size_t j = col_index;
		while(j != 0)
		{
			--j;
			if(!the_arithmetics().equal_zero(reduced[0][j]))
			{
				eigenvectors.push_back( this->find_eigenvector(reduced, j, col_index) );
			}
		}
	}
	return eigenvectors;
}

//returns if matrix is positive definite
//needs operator<= on T
template <typename T>
bool matrix_general<T>::positive_definite() const
{
	if(row_size() != col_size())
		return false;
	
	if(size() == 0)	
		return true;

	matrix_general<T> symetric = (*this + this->transposition());

	for(std::size_t i = 0; i<row_size(); ++i)
	{
		if(symetric[i][i] <= the_arithmetics().zero)
			return false;

		for(std::size_t row_index = i + 1; row_index < symetric.row_size(); ++row_index)
		{
			for(std::size_t col_index = i + 1; col_index < symetric.col_size(); ++col_index)
			{
				symetric[row_index][col_index] -= (symetric[row_index][i] * symetric[i][col_index]) / symetric[i][i];
			}
		}
	}

	return true;
}

template <typename T>
bool matrix_general<T>::positive_semi_definite() const
{	
	std::vector<std::complex<T>> eigenvalue = to_matrix_ptr()->eigenvalues();	
	for(auto it = eigenvalue.begin(); it != eigenvalue.end(); ++it)
	{
		if((it->real() < the_arithmetics().zero) || (!the_arithmetics().equal_zero(it->imag())))
			return false;
	}

	return true;
}

template <typename T>
class functor_less_for_complex 
{
public:
	bool operator() (const T & first, const T & second)		//T should be std::complex<?> 
	{
		if(first.real() == second.real())
			return first.imag() < second.imag();

		return first.real() < second.real();
	}
};

//returns jordan normal form of matrix
//eigenvalues should be accurate
template <typename T>
matrix<T> matrix_general<T>::jordan_form_impl(std::vector<T> eigenvalues) const
{
	if(row_size() != col_size())
		return matrix<T>(0, 0);

	std::map<T, std::size_t, functor_less_for_complex<T>> eigenvalues_multiplicity;

	matrix_general<T> jordan_form(row_size(), col_size(), 0);

	for(std::size_t i = 0; i<eigenvalues.size(); ++i)
		eigenvalues_multiplicity[eigenvalues[i]] += 1;		//value is number of occurences of key in eigenvalues

	std::size_t position = 0;
	for(auto it = eigenvalues_multiplicity.begin(); it != eigenvalues_multiplicity.end(); ++it)	//for all eigenvalues
	{
		std::vector<std::size_t> number_of_cells = numbers_of_jordans_cells(it->first, it->second); 
		
		for(std::size_t j=0; j<number_of_cells.size(); ++j)		//for all cells sizes with eigenvalue "it->first"
		{
			std::size_t cell_size = j + 1;
			for(std::size_t cells = 0; cells < number_of_cells[j]; ++cells)	//for all cell with eigenvalue "it->first" and size "cell_size"
			{
				jordan_form[position][position] = it->first;		//write one cell
				++position;

				for(std::size_t k = 1; k<cell_size; ++k)	
				{
					jordan_form[position][position] = it->first;
					jordan_form[position - 1][position] = the_arithmetics().one;
					++position;
				}
			}
		}
	}

	return jordan_form;
}

//computes number of jordan cells of all sizes with given eigenvalue. 
//Number of cells with size 1 is in first position (i.e. vector[0]) in vector, the size 2 in second, and so on.
template <typename T>
std::vector<std::size_t> matrix_general<T>::numbers_of_jordans_cells(T eigenvalue, std::size_t multiplicity) const
{
	std::size_t known_cells = 0;
	
	matrix_general<T> m = *this - eigenvalue * identity_matrix(row_size());
	matrix_general<T> m_power = m;
	
	std::vector<std::size_t> number_of_cells;

	std::vector<std::size_t> ranks;
	ranks.push_back(row_size());
	ranks.push_back(m.rank());
	
	std::size_t cell_size = 1;
	while(known_cells < multiplicity)
	{
		m_power *= m;
		ranks.push_back(m_power.rank());
		
		std::size_t last = ranks.size();
		number_of_cells.push_back( ranks[last - 3] - 2*ranks[last - 2] + ranks[last - 1]);

		known_cells += number_of_cells.back() * cell_size;
		++cell_size;
	}

	if(known_cells != multiplicity)
		throw matrix_exception("unknown error in jordan form computation");
	
	return number_of_cells;
}

template <typename T>
std::size_t matrix_general<T>::rank() const
{
	matrix_general<T> m = this->row_echelon_form();
							//last nonzero row in row_echelon_form determines rank			
	for(long long int i = m.row_size() - 1; i >= 0; --i)
	{	
		for(std::size_t j = 0; j<m.col_size(); ++j)
		{	
			if(!the_arithmetics().equal_zero(m[i][j]))
			{
				return i + 1;
			}
		}
	}

	return 0;
}

//first_substitution is n such that in first equation (first row) x_n was on the left side 
template <typename T>
std::vector<T> matrix_general<T>::find_eigenvector(const matrix_general<T> & m, std::size_t nonzero_variable, std::size_t first_substitution) const
{	
	std::vector<T> eigenvector(m.col_size(), 0);

	eigenvector[nonzero_variable] = the_arithmetics().one;
	eigenvector[first_substitution] = eigenvector[nonzero_variable] * m[0][nonzero_variable];

	for(std::size_t i = first_substitution + 1; i<m.col_size(); ++i)
	{	
		eigenvector[i] = eigenvector[i-1] * m[i - first_substitution][i-1];
	}
	
	return eigenvector;
}

//returns identity matrix with "row_size" rows and also "row_size" columns
template <typename T>
matrix<T> matrix_general<T>::identity_matrix(std::size_t row_size)
{
	matrix_general<T> m(row_size, row_size, the_arithmetics().zero);
	for(std::size_t i = 0; i<row_size; ++i)
		m[i][i] = the_arithmetics().one;

	return m;
}

//creates householder matrix from vector given by iterator to begining and end
template <typename T>
template <typename Iterator>
matrix<T> matrix_general<T>::householder_matrix(Iterator first, Iterator last) 
{
	matrix_general<T> m = outer_product<Iterator, Iterator>(first, last, first, last);
	T two = the_arithmetics().one + the_arithmetics().one;
	T t =  two / dot_product<Iterator, Iterator>(first, last, first, last);

	for(auto it = m.begin(); it != m.end(); ++it)
		(*it) = t * (*it);

	m = m - identity_matrix(m.row_size());

	return m;
}

//returns outer product of two vectors
template <typename T>
template <typename Iterator1, typename Iterator2>
matrix<T> matrix_general<T>::outer_product(Iterator1 first1, Iterator1 last1, Iterator2 first2, Iterator2 last2)
{
	int matrix_size = last1 - first1;
	
	matrix_general<T> m(matrix_size, matrix_size);
	
	iterator it = m.begin();
	for(Iterator1 it1 = first1; it1 != last1; ++it1)
	{
		for(Iterator2 it2 = first2; it2 != last2; ++it2)
		{
			*it = (*it1) * (*it2);
			++it;
		}
	}
	return m;
}

template <typename T>
template <typename Iterator1, typename Iterator2>
T matrix_general<T>::dot_product(Iterator1 first1, Iterator1 last1, Iterator2 first2, Iterator2 last2)
{
	T result = the_arithmetics().zero;
	
	while(first1 != last1)
	{
		result += (*first1) * (*first2);
		++first1;
		++first2;
	}
	
	return result;
}

//converts matrix_general<T> to matrix_general<complex<T>> with all imaginary numbers equal zero and same real numbers
template <typename T>
matrix<std::complex<T>> matrix_general<T>::to_complex() const
{
	matrix_general<std::complex<T>> complex_matrix(row_size(), col_size());

	matrix_general<T>::const_iterator real_it = cbegin();
	typename matrix_general<std::complex<T>>::iterator complex_it = complex_matrix.begin();

	for(; real_it != cend(); ++real_it, ++complex_it)
	{
		*complex_it = std::complex<T>(*real_it, the_arithmetics().zero);
	}
	return complex_matrix;
}

//adds multiplication of "row_to_add" to multiplication of "row_to_nullify" such that element on "first_non_zero_element" position in row_to_nullify will be zero.
//Rows are given by their indeces in matrix.
template <typename T>
void matrix_general<T>::nullify_first_element(std::size_t row_to_nullify, std::size_t row_to_add, std::size_t first_non_zero_element)
{
	matrix_general<T>& m = *this;
	T add_row_multiplicator = m[row_to_nullify][first_non_zero_element];		//this shouldnt be zero
	T nullify_row_multiplicator = m[row_to_add][first_non_zero_element];		//shouldnt be zero

	for(std::size_t i = first_non_zero_element; i<m.column_size(); ++i)	
		m[row_to_nullify][i] = m[row_to_nullify][i] * nullify_row_multiplicator - m[row_to_add][i] * add_row_multiplicator;
	
}

//The same as above but throws exception when arithmetics overflow occurs.
template <typename T>
void matrix_general<T>::checked_nullify_first_element(std::size_t row_to_nullify, std::size_t row_to_add, std::size_t first_non_zero_element)
{
	matrix_general<T>& m = *this;
	T add_row_multiplicator = m[row_to_nullify][first_non_zero_element];		//this shouldnt be zero
	T nullify_row_multiplicator = m[row_to_add][first_non_zero_element];		//shouldnt be zero

	for(std::size_t i = first_non_zero_element; i<m.column_size(); ++i)	
	{
		T first = the_arithmetics().checked_multiplication(m[row_to_nullify][i], nullify_row_multiplicator);
		T second = the_arithmetics().checked_multiplication(m[row_to_add][i], add_row_multiplicator);
		m[row_to_nullify][i] = the_arithmetics().checked_subtraction(first, second);
	}
}

template <typename T>
void matrix_general<T>::checked_add_row(matrix_general<T>::row row_to_add, std::size_t row_index)
{
	matrix_general<T>& m = *this;
	if(row_to_add.size() != m.row_size())
		throw matrix_exception("Cannot add row to matrix of different row size.");

	for(std::size_t i=0; i<m.column_size(); ++i)
		m[row_index][i] = the_arithmetics().checked_addition(m[row_index][i], row_to_add[i]);
}

template <typename T>
void matrix_general<T>::add_row(matrix_general<T>::row row_to_add, std::size_t row_index)
{
	matrix_general<T>& m = *this;
	if(row_to_add.size() != m.row_size())
		throw matrix_exception("Cannot add row to matrix of different row size.");

	for(std::size_t i=0; i<m.column_size(); ++i)
		m[row_index][i] += row_to_add[i];
}

template <typename T>
void matrix_general<T>::checked_multiply_row(T multiplicator, std::size_t row_index)
{
	matrix_general<T>& m = *this;
	for(std::size_t i=0; i<m.column_size(); ++i)
		m[row_index][i] = the_arithmetics().checked_multiplication(m[row_index][i],  multiplicator);
}

template <typename T>
void matrix_general<T>::multiply_row(T multiplicator, std::size_t row_index)
{
	matrix_general<T>& m = *this;
	for(std::size_t i=0; i<m.column_size(); ++i)
		m[row_index][i] *= multiplicator;
}

template <typename T>
void matrix_general<T>::add_row_multiplication(T multiplicator, std::size_t row_to_add, std::size_t destination_row)
{
	matrix_general<T>& m = *this;
	for(std::size_t i=0; i<m.column_size(); ++i)
		m[destination_row][i] += m[row_to_add][i] * multiplicator;
}

template <typename T>
void matrix_general<T>::swap_rows_elements(std::size_t first_row_index, std::size_t second_row_index)
{
	matrix_general<T>& m = *this;
	for(std::size_t i=0; i<m.row_size(); ++i)
	{
		T tmp = m[first_row_index][i];
		m[first_row_index][i] = m[second_row_index][i];
		m[second_row_index][i] = tmp;
	}
}

#endif 
