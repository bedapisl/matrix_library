//Zapoctovy program pro predmet Programovani v C++ v roce 2013/2014
//Bedrich Pisl

#ifndef BEDAS_ITERATORS_H
#define BEDAS_ITERATORS_H

template <typename T, typename real_type>		//this class should be abstract base class and "real_type" should always be type of derived iterator
class base_iterator : public std::iterator<std::random_access_iterator_tag, T>
{
public:		
	real_type& operator++ ()		//prefix
	{
		position += difference;
		return *( (real_type*) this);
	}
	real_type operator++ (int)		//postfix
	{
		real_type result(base, position);
		position += difference;
		return result;
	}
	friend real_type operator+ (real_type it, int distance)
	{
		real_type result = it;
		result.position += result.difference * distance;
		return result;
	}
	friend real_type operator+ (int distance, real_type it)
	{
		return it + distance;
	}
	real_type& operator+= (int distance)
	{
		position += distance * difference;
		return *((real_type*) this);
	}
	real_type& operator-- ()
	{
		position -= difference;
		return *this;
	}
	real_type operator-- (int)
	{
		real_type result(base, position);
		position -= difference;
		return result;
	}
	real_type& operator-= (int parameter)
	{
		return *this += (-parameter);
	}
	friend real_type operator- (real_type it, int distance)
	{
		it.position -= distance * it.difference;
		return it;
	}
	friend int operator- (real_type first, real_type second)
	{
		return (first.position - second.position) / first.difference;
	}	
	bool operator==	(real_type other)
	{
		return position == other.position;
	}
	bool operator!= (real_type other)
	{
		return position != other.position;
	}
	bool operator< (real_type other)
	{
		return position < other.position;
	}
	bool operator> (real_type other)
	{
		return position > other.position;
	}
	bool operator>= (real_type other)
	{
		return position >= other.position;
	}
	bool operator<= (real_type other)
	{
		return position <= other.position;
	}

protected:
	base_iterator(std::shared_ptr<typename matrix_general<T>::matrix_base> base, int position, int difference) : base(base), position(position), difference(difference)
	{}
	std::shared_ptr<typename matrix_general<T>::matrix_base> base;
	int position;
	int difference;
};

//all elements iterator
template <typename T>
class matrix_general<T>::iterator : public base_iterator<T, typename matrix_general<T>::iterator>
{
public:
	iterator(std::shared_ptr<matrix_base> base, int position) : base_iterator<T, typename matrix_general<T>::iterator>(base, position, 1) {}
	T& operator* () {return *(this->base->data + this->position);}
	T* operator-> () {return (this->base->data + this->position);}		
	T& operator[] (int index) {return *(this->base->data + this->position);}
};

//const version
template <typename T>
class matrix_general<T>::const_iterator : public base_iterator<T, typename matrix_general<T>::const_iterator>
{
public:
	const_iterator(std::shared_ptr<matrix_base> base, int position) : base_iterator<T, typename matrix_general<T>::const_iterator>(base, position, 1) {}
	const T& operator* () {return *(this->base->data + this->position);}
	const T* operator-> () {return (this->base->data + this->position);}
	const T& operator[] (int index) {return *(this->base->data + this->position);}
};

//elements in row iterator
template <typename T>
class matrix_general<T>::row::iterator : public base_iterator<T, typename matrix_general<T>::row::iterator>
{
public:
	iterator(std::shared_ptr<matrix_base> base, int position) : base_iterator<T, typename matrix_general<T>::row::iterator>(base, position, 1) {}
	T& operator* () {return *(this->base->data + this->position);}
	T* operator-> () {return (this->base->data + this->position);}		
	T& operator[] (int index) {return *(this->base->data + this->position);}
};

//const version
template <typename T>
class matrix_general<T>::row::const_iterator : public base_iterator<T, typename matrix_general<T>::row::const_iterator>
{
public:
	const_iterator(std::shared_ptr<matrix_base> base, int position) : base_iterator<T, typename matrix_general<T>::row::const_iterator>(base, position, 1) {}
	const T& operator* () {return *(this->base->data + this->position);}
	const T* operator-> () {return (this->base->data + this->position);}		
	const T& operator[] (int index) {return *(this->base->data + this->position);}
};

//elements in column iterator
template <typename T>
class matrix_general<T>::col::iterator : public base_iterator<T, typename matrix_general<T>::col::iterator>
{
public:
	iterator(std::shared_ptr<matrix_base> base, int position) : base_iterator<T, typename matrix_general<T>::col::iterator>(base, position, base->row_length) {}
	T& operator* () {return *(this->base->data + this->position);}
	T* operator-> () {return (this->base->data + this->position);}		
	T& operator[] (int index) {return *(this->base->data + this->position);}
};

//const version
template <typename T>
class matrix_general<T>::col::const_iterator : public base_iterator<T, typename matrix_general<T>::col::const_iterator>
{
public:
	const_iterator(std::shared_ptr<matrix_base> base, int position) : base_iterator<T, typename matrix_general<T>::col::const_iterator>(base, position, base->row_length) {}
	const T& operator* () {return *(this->base->data + this->position);}
	const T* operator-> () {return (this->base->data + this->position);}		
	const T& operator[] (int index) {return *(this->base->data + this->position);}
};

//rows iterator
template <typename T>
class matrix_general<T>::row_iterator : public base_iterator<T, typename matrix_general<T>::row_iterator>
{
public:
	row_iterator(std::shared_ptr<matrix_base> base, int position) : base_iterator<T, typename matrix_general<T>::row_iterator>(base, position, base->row_length) {}
	row operator* () {return row(this->base, this->position);}
	row* operator-> () 
	{
		static void* aux_row = operator new(sizeof(row));	//allocates memory
		aux_row = new (aux_row) row(this->base, this->position);	//call constructor
		return (matrix_general<T>::row*)aux_row;
	}		
	row operator[] (int index) {return *(this->base->data + this->position);}
};

//const version
template <typename T>
class matrix_general<T>::row_const_iterator : public base_iterator<T, typename matrix_general<T>::row_const_iterator>
{
public:
	row_const_iterator(std::shared_ptr<matrix_base> base, int position) : base_iterator<T, typename matrix_general<T>::row_const_iterator>(base, position, base->row_length) {}
	const row operator* () {return row(this->base, this->position);}
	const row* operator-> () 
	{
		static void* aux_row = operator new(sizeof(row));	//allocates memory
		aux_row = new (aux_row) row(this->base, this->position);	//call constructor
		return (matrix_general<T>::row*)aux_row;
	}		
	const row operator[] (int index) {return *(this->base->data + this->position);}
};

//columns iterator
template <typename T>
class matrix_general<T>::col_iterator : public base_iterator<T, typename matrix_general<T>::col_iterator>
{
public:
	col_iterator(std::shared_ptr<matrix_base> base, int position) : base_iterator<T, typename matrix_general<T>::col_iterator>(base, position, 1)
	{}
	col operator* () {return col(this->base, this->position);}
	col* operator-> () 
	{
		static void* aux_col = operator new(sizeof(col));	//allocates memory
		aux_col = new (aux_col) col(this->base, this->position);	//call constructor
		return (matrix_general<T>::col*)aux_col;
	}		
	col operator[] (int index) {return *(this->base->data + this->position);}
};

//const verison
template <typename T>
class matrix_general<T>::col_const_iterator : public base_iterator<T, typename matrix_general<T>::col_const_iterator>
{
public:
	col_const_iterator(std::shared_ptr<matrix_base> base, int position) : base_iterator<T, typename matrix_general<T>::col_const_iterator>(base, position, 1)
	{}
	const col operator* () {return col(this->base, this->position);}
	const col* operator-> () 
	{
		static void* aux_col = operator new(sizeof(col));	//allocates memory
		aux_col = new (aux_col) col(this->base, this->position);	//call constructor
		return (matrix_general<T>::col*)aux_col;
	}		
	const col operator[] (int index) {return *(this->base->data + this->position);}
};

#endif

