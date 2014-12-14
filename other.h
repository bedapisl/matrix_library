//Zapoctovy program pro predmet Programovani v C++ v roce 2013/2014
//Bedrich Pisl

#ifndef BEDAS_OTHER_H
#define BEDAS_OTHER_H

class matrix_exception : public std::runtime_error
{
public:
	matrix_exception(std::string s) : std::runtime_error(s) {}
};

template <typename T>
class matrix_general<T>::row
{
public:	
	class iterator;
	class const_iterator;

	typedef T value_type;
	typedef T& reference;
	typedef const T& const_reference;
	typedef std::ptrdiff_t difference_type;
	typedef std::size_t size_type;

	row() : base(new matrix_base(NULL, 0, 0)) {}
	row(const row& other) : base(other.base) {}
	row(std::shared_ptr<matrix_base> base, std::size_t position) : base(base), position(position) {}
	row& operator= (row other);

	bool operator== (typename matrix_general<T>::row & other) const;
	bool operator!= (typename matrix_general<T>::row & other) const {return *this != other;}

	T& operator[] (std::size_t index) {return *(base->data + position + index);}
	const T& operator[] (std::size_t index) const {return *(base->data + position + index);}
	T& at(std::size_t index);
	const T& at(std::size_t index) const {return at(index);}

	iterator begin(){return iterator(base, position);}
	iterator end(){return iterator(base, position + base->row_length);}
	const_iterator begin() const {return const_iterator(base, position);}
	const_iterator end() const {return const_iterator(base, position + base->row_length);}
	const_iterator cbegin() const {return const_iterator(base, position);}
	const_iterator cend() const {return const_iterator(base, position + base->row_length);}

	T& front() {return *(base->data + position);}
	T& back() {return *(base->data + position + base->row_length - 1);}
	const T& front() const {return front();}
	const T& back() const {return back();}

	std::size_t size() const {return base->row_length;}
	bool empty() const {return size() == 0;}

private:
	std::shared_ptr<matrix_base> base;
	std::size_t position;
};

template <typename T>
typename matrix_general<T>::row& matrix_general<T>::row::operator= (row other)
{
	if(base->row_length != other.base->row_length)
		throw matrix_exception("wrong row length");

	for(std::size_t i=0; i<base->row_length; ++i)
		base->data[position + i] = other.base->data[other.position + i];

	return *this;
}

template <typename T>
bool matrix_general<T>::row::operator== (typename matrix_general<T>::row & other) const 
{
	if(base->row_length != other.base->row_length)
		return false;
	
	for(std::size_t i=0; i<base->row_length; ++i)
		if(base->data[i] != other.base->data[i])
			return false;
	return true;
}

template <typename T>
T& matrix_general<T>::row::at(std::size_t index)
{	
	if(index >= base->row_length)
		throw std::out_of_range("matrix_general<T>::row::at - range check failed");
	
	return (*this)[index];
}

template <typename T>
class matrix_general<T>::col
{
public:	
	class iterator;
	class const_iterator;

	typedef T value_type;
	typedef T& reference;
	typedef const T& const_reference;
	typedef std::ptrdiff_t difference_type;
	typedef std::size_t size_type;

	col() : base(new matrix_base(NULL, 0, 0)) {}
	col(const col& other) : base(other.base) {}
	col(std::shared_ptr<matrix_base>  base, std::size_t position) : base(base), position(position) {}
	col& operator= (col& other);
	
	bool operator== (typename matrix_general<T>::col & other) const;
	bool operator!= (typename matrix_general<T>::col & other) const {return *this != other;}

	T& operator[] (std::size_t index) {return *(base->data + position + index*base->row_length);}
	const T& operator[] (std::size_t index) const {return *(base->data + position + index*base->row_length);}
	T& at(std::size_t index);
	const T& at(std::size_t index) const {return at(index);}

	iterator begin() {return iterator(base, position);}
	iterator end() {return iterator(base, position + (base->col_length * base->row_length));}	
	const_iterator begin() const {return const_iterator(base, position);}
	const_iterator end() const {return const_iterator(base, position + (base->col_length * base->row_length));}
	const_iterator cbegin() const {return const_iterator(base, position);}
	const_iterator cend() const {return const_iterator(base, position + (base->col_length * base->row_length));}

	T& front() {return *(base->data + position);}
	T& back() {return *(base->data + position + (base->col_length - 1) * base->row_length);}
	const T& front() const {return front();}
	const T& back() const {return back();}

	std::size_t size() {return base->col_length;}
	bool empty() {return size() == 0;}

private:
	std::shared_ptr<matrix_base> base;
	std::size_t position;
};

template <typename T>
typename matrix_general<T>::col& matrix_general<T>::col::operator= (col& other)
{
	if(base->col_length != other.base->col_length)
		throw matrix_exception("wrong column length");

	for(std::size_t i=0; i<base->col_length; ++i)
		base->data[position + i * base->row_length] = other.base->data[position + i * base->row_length];

	return *this;
}

template <typename T>
bool matrix_general<T>::col::operator== (typename matrix_general<T>::col & other) const 
{
	if(base->col_length != other.base->col_length)
		return false;
	
	for(std::size_t i=0; i<base->col_length; ++i)
		if(base->data[i] != other.base->data[i])
			return false;
	return true;
}

template <typename T>
T& matrix_general<T>::col::at(std::size_t index) 
{
	if(index >= base->col_length)
		throw std::out_of_range("matrix::col::at - range checked failed");
	
	return (*this)[index];
}

template <typename T>
class matrix_general<T>::matrix_base
{
public:
	T* data;
	std::size_t row_length;
	std::size_t col_length;

	friend matrix_general<T>;

	matrix_base(T* data, std::size_t row_length, std::size_t col_length) : data(data), row_length(row_length), col_length(col_length) {}
	matrix_base(const matrix_base& other);
	~matrix_base();
	matrix_general<T>::matrix_base deep_copy() const;
	bool operator== (const matrix_general<T>::matrix_base& other);
};

template <typename T>
bool matrix_general<T>::matrix_base::operator== (const matrix_general<T>::matrix_base& other)
{
	if((row_length != other.row_length) || (col_length != other.col_length))
		return false;
	
	for(std::size_t i=0; i<row_length * other.col_length; ++i)
	{
		if(data[i] != other.data[i])
			return false;
	}
	return true;
}

template <typename T>
matrix_general<T>::matrix_base::matrix_base(const matrix_general<T>::matrix_base& other) 
			: data(new T[other.row_length * other.col_length]), row_length(other.row_length), col_length(other.col_length)
{
	for(std::size_t i=0; i<row_length * col_length; ++i)
		data[i] = other.data[i];
}

template <typename T>
matrix_general<T>::matrix_base::~matrix_base()
{
	if(data != NULL)
		delete[] data;
}

template <typename T>
typename matrix_general<T>::matrix_base matrix_general<T>::matrix_base::deep_copy() const 
{
	T* new_data = new T [row_length * col_length];
	for(std::size_t i=0; i<row_length * col_length; i++)
		new_data[i] = data[i];
	
	return matrix_base(new_data, row_length, col_length);
}

template <typename T>
arithmetics<T, std::numeric_limits<T>::is_specialized> & matrix_general<T>::the_arithmetics()
{
	static arithmetics<T, std::numeric_limits<T>::is_specialized> instance;
	return instance;
}

//"singleton" class for computations with given type.
template <typename T, bool built_in_type>		//for user-defined types. Specialization for built-in types is below.
class arithmetics
{
public:
	arithmetics()
	{
		one = 1;
	}
	T checked_addition(const T& first, const T& second) {return first + second;}
	T checked_subtraction(const T& first, const T& second) {return first - second;}
	T checked_multiplication(const T& first, const T& second) {return first * second;}

	T square_root(const T& number)	
	{
		if(number == zero)
			return number;

		T estimate = one;
		T old_estimate;
		T two = one + one;

		while(old_estimate != estimate)
		{
			old_estimate = estimate;
			estimate = (old_estimate + number/old_estimate) / two;
		}
		
		return estimate;
	}
	
	std::size_t find_nonzero_row(const matrix_general<T>& m, std::size_t column_index, std::size_t first_row, bool & only_zeros)
	{
		only_zeros = true;
		for(std::size_t i = first_row; i < m.col_size(); ++i)
		{
			if(m[i][column_index] != zero)
			{
				only_zeros = false;
				return i;
			}
		}
		return 0;
	}
	
	bool equal_zero(const std::complex<T> & number)
	{
		return (number == std::complex<T>(zero, zero));
	}

	T zero;
	T one;
};

template <typename T>			//partial specialization for built-in types
class arithmetics<T, true> 
{
public:
	arithmetics()
	{
		one = 1;

		max = std::numeric_limits<T>::max();
		if(std::numeric_limits<T>::is_integer)
			min = std::numeric_limits<T>::min();

		else
			min = -max;
	}
	
	T checked_addition(T first, T second)
	{
		if(((second >= zero) && (first > max - second)) || ((second < zero) && (first < min - second)))
			throw std::range_error("addition caused arithmetic overflow");
		return first + second;
	}
	T checked_subtraction(T first, T second)
	{
		if(((second >= zero) && (first < min + second)) || ((second < zero) && (first > max + second)))
			throw std::range_error("subtraction caused arithmetic overflow");
		return first - second;
	}
	T checked_multiplication(T first, T second)
	{	
		if(second == zero)
			return zero;
		if(std::abs(first) > max/std::abs(second))
			throw std::range_error("multiplication caused arithmetic overflow");
		
		return first*second;
	}
	
	T square_root(T number)
	{
		if(number == zero)
			return number;

		T estimate = one;
		T old_estimate;
		T two = one + one;

		while(old_estimate != estimate)
		{
			old_estimate = estimate;
			estimate = (old_estimate + number/old_estimate) / two;
		}
		
		return estimate;
	}
	
	std::size_t find_nonzero_row(const matrix_general<T>& m, std::size_t column_index, std::size_t first_row, bool & only_zeros)
	{	
		only_zeros = true;
		std::size_t max_index = 0;
		T max_abs_value = 0;
		for(std::size_t i = first_row; i < m.col_size(); ++i)
		{
			if(m[i][column_index] != zero)
			{
				if(std::abs(m[i][column_index]) > max_abs_value)
				{
					max_abs_value = std::abs(m[i][column_index]);
					max_index = i;
					only_zeros = false;
				}
			}
		}

		return max_index;
	}

	bool equal_zero(const std::complex<T> & number)
	{
		return number == std::complex<T>(zero, zero);
	}
	
	T one;
	T zero;
	T max;
	T min;
};

template <>
inline bool arithmetics<double, true>::equal_zero(const std::complex<double> & number)
{
	return std::abs(number) < 0.00001;
}

template <>
inline bool arithmetics<long double, true>::equal_zero(const std::complex<long double> & number)
{
	return std::abs(number) < 0.00001;
}

template <>
inline bool arithmetics<float, true>::equal_zero(const std::complex<float> & number)
{
	return std::abs(number) < 0.00001;
}

template <>
inline double arithmetics<double, true>::square_root(double number)
{
	return std::sqrt(number);
}

template <>
inline float arithmetics<float, true>::square_root(float number)
{
	return std::sqrt(number);
}

template <>
inline long double arithmetics<long double, true>::square_root(long double number)
{
	return std::sqrt(number);
}


#endif
