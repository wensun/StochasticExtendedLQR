/*
* matrix.h
* ExtendedLQR Library
*
* Copyright (c) 2013 University of Utah.
* All rights reserved.
*
* Permission to use, copy, modify, and distribute this software and its
* documentation for educational, research, and non-profit purposes, without
* fee, and without a written agreement is hereby granted, provided that the
* above copyright notice, this paragraph, and the following four paragraphs
* appear in all copies.
*
* Permission to incorporate this software into commercial products may be
* obtained by contacting the authors <berg@cs.utah.edu> or the Technology
* and Venture Commercialization office at the University of Utah
* <801-581-7792> <http://tvc.utah.edu>.
*
* This software program and documentation are copyrighted by the University of
* Utah. The software program and documentation are supplied "as is," without 
* any accompanying services from the University of Utah or the authors. The 
* University of Utah and the authors do not warrant that the operation of
* the program will be uninterrupted or error-free. The end-user understands
* that the program was developed for research purposes and is advised not to
* rely exclusively on the program for any reason.
*
* IN NO EVENT SHALL THE UNIVERSITY OF UTAH OR THE AUTHORS BE LIABLE TO ANY 
* PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, 
* INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS 
* DOCUMENTATION, EVEN IF THE UNIVERSITY OF UTAH OR THE AUTHORS HAVE BEEN 
* ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
* THE UNIVERSITY OF UTAH AND THE AUTHORS SPECIFICALLY DISCLAIM ANY WARRANTIES, 
* INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
* FITNESS FOR A PARTICULAR PURPOSE AND ANY STATUTORY WARRANTY OF NON-INFRINGEMENT. 
* THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF 
* UTAH AND THE AUTHORS HAVE NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, 
* UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
*
* Please send all bug reports to <berg@cs.unc.edu>.
*
* The authors may be contacted via:
*
* Jur van den Berg
* School of Computing
* 50 S. Central Campus Drive, 
* Salt Lake City, UT 84112 
* United States of America
*
* <http://arl.cs.utah.edu/research/extendedlqr/>
*/

#ifndef __MATRIX_H__
#define __MATRIX_H__

#define _USE_MATH_DEFINES

//#include <cmath>
#include <math.h>
#include <float.h>
#include <iostream>
#include <assert.h>
#include <limits>
#include <iomanip>

template <size_t _size> class SymmetricMatrix;

template <size_t _numRows, size_t _numColumns = 1> class Matrix {

private:
	double _elems[_numRows * _numColumns];

public:
	// Retrieval
	inline size_t numRows() const { 
		return _numRows; 
	}
	inline size_t numColumns() const { 
		return _numColumns; 
	}

	// Subscript operator
	__forceinline double& operator () (size_t row, size_t column) {
		assert(row < _numRows && column < _numColumns);
		return _elems[row * _numColumns + column]; 
	}
	__forceinline double  operator () (size_t row, size_t column) const {
		assert(row < _numRows && column < _numColumns);
		return _elems[row * _numColumns + column]; 
	}

	__forceinline double& operator [] (size_t elt) {
		assert(elt < _numRows * _numColumns);
		return _elems[elt]; 
	}
	__forceinline double  operator [] (size_t elt) const {
		assert(elt < _numRows * _numColumns);
		return _elems[elt]; 
	}

	// Reset to zeros
	inline void reset() { 
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			_elems[i] = double(0);
		}
	}

	// Submatrix
	template <size_t numRows, size_t numColumns>
	inline Matrix<numRows, numColumns> subMatrix(size_t row, size_t column) const {
		assert(row + numRows <= _numRows && column + numColumns <= _numColumns);
		Matrix<numRows, numColumns> m;
		for (size_t i = 0; i < numRows; ++i) {
			for (size_t j = 0; j < numColumns; ++j) {
				m(i, j) = (*this)(row + i, column + j);
			}
		}
		return m;
	}

	template <size_t numRows>
	inline Matrix<numRows> subMatrix(size_t row, size_t column) const {
		assert(row + numRows <= _numRows && column < _numColumns);
		Matrix<numRows> m;
		for (size_t i = 0; i < numRows; ++i) {
			m[i] = (*this)(row + i, column);
		}
		return m;
	}

	inline Matrix<_numRows> column(size_t columnNr) const {
		assert(columnNr < _numColumns);
		Matrix<_numRows> m;
		for (size_t i = 0; i < _numRows; ++i) {
			m[i] = (*this)(i, columnNr);
		}
		return m;
	}

	inline Matrix<1, _numColumns> row(size_t rowNr) const {
		assert(rowNr < _numRows);
		Matrix<1, _numColumns> m;
		for (size_t i = 0; i < _numColumns; ++i) {
			m[i] = (*this)(rowNr, i);
		}
		return m;
	}

	// Insert
	template <size_t numRows, size_t numColumns>
	inline void insert(size_t row, size_t column, const Matrix<numRows, numColumns>& q) {
		assert(row + numRows <= _numRows && column + numColumns <= _numColumns);
		for (size_t i = 0; i < numRows; ++i) {
			for (size_t j = 0; j < numColumns; ++j) {
				(*this)(row + i, column + j) = q(i, j);
			}
		}
	}

	template <size_t numRows>
	inline void insert(size_t row, size_t column, const SymmetricMatrix<numRows>& q) {
		assert(row + numRows <= _numRows && column + numRows <= _numColumns);
		for (size_t i = 0; i < numRows; ++i) {
			for (size_t j = 0; j < numRows; ++j) {
				(*this)(row + i, column + j) = q(i, j);
			}
		}
	}

	// Matrix addition
	inline const Matrix<_numRows, _numColumns>& operator+=(const Matrix<_numRows, _numColumns>& q) { 
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			_elems[i] += q._elems[i];
		}
		return *this; 
	}

	// Matrix subtraction
	inline const Matrix<_numRows, _numColumns>& operator-=(const Matrix<_numRows, _numColumns>& q) { 
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			_elems[i] -= q._elems[i];
		}
		return *this; 
	}

	// Scalar multiplication
	inline const Matrix<_numRows, _numColumns>& operator*=(double a) { 
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			_elems[i] *= a;
		}
		return *this;
	}

	// Scalar division
	inline const Matrix<_numRows, _numColumns>& operator/=(double a) { 
		for (size_t i = 0; i < _numRows * _numColumns; ++i) {
			_elems[i] /= a;
		}
		return *this;
	}
};

template <size_t _size> class SymmetricMatrix {
private:
	double _elems[((_size+1)*_size)/2];

public:
	// Retrieval
	inline size_t numRows() const { 
		return _size; 
	}
	inline size_t numColumns() const { 
		return _size; 
	}

	// Subscript operator
	__forceinline double& operator () (size_t row, size_t column) {
		assert(row < _size && column < _size);
		if (row >= column) {
			return _elems[_size * column + row - ((column + 1)*column) / 2];
		} else {
			return _elems[_size * row + column - ((row + 1)*row) / 2];
		}
	}
	__forceinline double operator () (size_t row, size_t column) const {
		assert(row < _size && column < _size);
		if (row >= column) {
			return _elems[_size * column + row - ((column + 1)*column) / 2];
		} else {
			return _elems[_size * row + column - ((row + 1)*row) / 2];
		}
	}

	__forceinline double& operator [] (size_t elt) {
		assert(elt < ((_size+1)*_size)/2);
		return _elems[elt]; 
	}
	__forceinline double  operator [] (size_t elt) const {
		assert(elt < ((_size+1)*_size)/2);
		return _elems[elt]; 
	}

	inline void reset() { 
		for (size_t i = 0; i < ((_size+1)*_size)/2; ++i) {
			_elems[i] = double(0);
		}
	}

	// Matrix addition
	inline const SymmetricMatrix<_size>& operator+=(const SymmetricMatrix<_size>& M) { 
		for (size_t i = 0; i < ((_size+1)*_size)/2; ++i) {
			_elems[i] += M._elems[i];
		}
		return *this;
	}

	// Matrix subtraction
	inline const SymmetricMatrix<_size>& operator-=(const SymmetricMatrix<_size>& M) { 
		for (size_t i = 0; i < ((_size+1)*_size)/2; ++i) {
			_elems[i] -= M._elems[i];
		}
		return *this;
	}

	// Scalar multiplication
	inline const SymmetricMatrix<_size>& operator*=(double a) { 
		for (size_t i = 0; i < ((_size+1)*_size)/2; ++i) {
			_elems[i] *= a;
		}
		return *this;
	}

	// Scalar division
	inline const SymmetricMatrix<_size>& operator/=(double a) { 
		for (size_t i = 0; i < ((_size+1)*_size)/2; ++i) {
			_elems[i] /= a;
		}
		return *this;
	}

	// Extract symmetric subMatrix starting from a given diagonal element 
	template <size_t _numRows>
	inline SymmetricMatrix<_numRows> subSymmetricMatrix(size_t diag) const {
		assert(diag + _numRows <= _size);
		SymmetricMatrix<_numRows> L;
		size_t index = 0;
		for (size_t i = 0; i < _numRows; ++i) {
			for(size_t j = 0; j <= i; ++j) {
				L(i,j) = (*this)(diag + i, diag + j);
			}
		}
		return L;
	}

	template <size_t _numRows>
	inline void insert(size_t diag, const SymmetricMatrix<_numRows>& M) {
		assert(diag + _numRows <= _size);
		for (size_t i = 0; i < _numRows; ++i) {
			for(size_t j = 0; j <= i; ++j) {
				(*this)(diag + i, diag + j) = M(i,j);
			}
		}
	}

	// Cast
	inline operator Matrix<_size, _size>() const {
		Matrix<_size,_size> M;
		for (size_t i = 0; i < _size; ++i) {
			for (size_t j = 0; j < _size; ++j) {
				M(i,j) = (*this)(i,j);
			}
		}
		return M;
	}
};

// Unary minus
template <size_t _numRows, size_t _numColumns>
inline Matrix<_numRows, _numColumns> operator-(const Matrix<_numRows, _numColumns>& M) {
	Matrix<_numRows, _numColumns> L;
	for (size_t i = 0; i < _numRows * _numColumns; ++i) {
		L[i] = -M[i];
	}
	return L;
}

template <size_t _size>
inline SymmetricMatrix<_size> operator-(const SymmetricMatrix<_size>& M) {
	SymmetricMatrix<_size> L;
	for (size_t i = 0; i < ((_size+1)*_size)/2; ++i) {
		L[i] = -M[i];
	}
	return L;
}

// Unary plus
template <size_t _numRows, size_t _numColumns>
inline const Matrix<_numRows, _numColumns>& operator+(const Matrix<_numRows, _numColumns>& M) { 
	return M; 
}

template <size_t _size>
inline const SymmetricMatrix<_size>& operator+(const SymmetricMatrix<_size>& M) { 
	return M; 
}

// Transpose
template <size_t _numRows, size_t _numColumns>
inline Matrix<_numColumns, _numRows> operator~(const Matrix<_numRows, _numColumns>& M) {
	Matrix<_numColumns, _numRows> L;
	for (size_t i = 0; i < _numColumns; ++i) {
		for (size_t j = 0; j < _numRows; ++j) {
			L(i, j) = M(j, i);
		}
	}
	return L;
}

template <size_t _size>
inline const SymmetricMatrix<_size>& operator~(const SymmetricMatrix<_size>& M) { 
	return M; 
}

// Matrix trace
template <size_t _size>
inline double tr(const Matrix<_size, _size>& q) { 
	double trace = double(0);
	for (size_t i = 0; i < _size; ++i){
		trace += q(i, i);
	}
	return trace;
}

template <size_t _size>
inline double tr(const SymmetricMatrix<_size>& q) { 
	double trace = double(0);
	for (size_t i = 0; i < _size; ++i){
		trace += q(i, i);
	}
	return trace;
}


// Trace of a product p*q
template <size_t _numRows, size_t _numColumns>
inline double trProd(const Matrix<_numRows, _numColumns>& p, const Matrix<_numColumns, _numRows>& q) { 
	double trace = double(0);
	for (size_t i = 0; i < _numRows; ++i) {
		for (size_t j = 0; j < _numColumns; ++j) {
			trace += p(i, j) * q(j, i);
		}
	}
	return trace;
}

template <size_t _size>
inline double trProd(const SymmetricMatrix<_size>& p, const SymmetricMatrix<_size>& q) { 
	double trace = double(0);
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			trace += p(i, j) * q(j, i);
		}
	}
	return trace;
}

template <size_t _size>
inline double trProd(const Matrix<_size, _size>& p, const SymmetricMatrix<_size>& q) { 
	double trace = double(0);
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			trace += p(i, j) * q(j, i);
		}
	}
	return trace;
}

template <size_t _size>
inline double trProd(const SymmetricMatrix<_size>& p, const Matrix<_size, _size>& q) { 
	double trace = double(0);
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			trace += p(i, j) * q(j, i);
		}
	}
	return trace;
}


// Matrix 1-norm
template <size_t _numRows, size_t _numColumns>
inline double norm(const Matrix<_numRows, _numColumns>& q) {
	double norm1 = double(0);
	for (size_t j = 0; j < _numColumns; ++j) {
		double colabssum = double(0);
		for (size_t i = 0; i < _numRows; ++i) {
			colabssum += abs(q(i,j));
		}
		if (colabssum > norm1) {
			norm1 = colabssum;
		}
	}
	return norm1;
}

// Identity matrix
template <size_t _size>
inline SymmetricMatrix<_size> identity() {
	SymmetricMatrix<_size> m;
	for (size_t j = 0; j < _size; ++j) {
		for (size_t i = j; i < _size; ++i) {
			m(i,j) = (i == j ? 1.0 : 0.0);
		}
	}
	return m;
}

// Zero matrix
template <size_t _numRows, size_t _numColumns>
inline Matrix<_numRows, _numColumns> zeros() {
	Matrix<_numRows, _numColumns> m;
	m.reset();
	return m;
}

template <size_t _size>
inline SymmetricMatrix<_size> zeros() {
	SymmetricMatrix<_size> L;
	L.reset();
	return L;
}

template <size_t _numRows>
inline Matrix<_numRows> zero() {
	Matrix<_numRows> m;
	m.reset();
	return m;
}

// Matrix determinant
template <size_t _size>
inline double det(const Matrix<_size, _size>& q) { 
	Matrix<_size, _size> m(q);
	double D = double(1);

	size_t row_p[_size];
	size_t col_p[_size];
	for (size_t i = 0; i < _size; ++i) {
		row_p[i] = i; col_p[i] = i;
	}

	// Gaussian elimination
	for (size_t k = 0; k < _size; ++k) {
		// find maximal pivot element
		double maximum = double(0); size_t max_row = k; size_t max_col = k;
		for (size_t i = k; i < _size; ++i) {
			for (size_t j = k; j < _size; ++j) {
				double abs_ij = std::abs(m(row_p[i], col_p[j]));
				if (abs_ij > maximum) {
					maximum = abs_ij; max_row = i; max_col = j;
				}
			}
		}

		// swap rows and columns
		if (k != max_row) {
			size_t swap = row_p[k]; row_p[k] = row_p[max_row]; row_p[max_row] = swap;
			D = -D;
		}
		if (k != max_col) {
			size_t swap = col_p[k]; col_p[k] = col_p[max_col]; col_p[max_col] = swap;
			D = -D;
		}

		D *= m(row_p[k], col_p[k]);
		if (D == double(0)) {
			return double(0);
		}

		// eliminate column
		for (size_t i = k + 1; i < _size; ++i) {
			double factor = m(row_p[i], col_p[k]) / m(row_p[k], col_p[k]);
			for (size_t j = k + 1; j < _size; ++j) {
				m(row_p[i], col_p[j]) -= factor * m(row_p[k], col_p[j]);
			}
		}  
	}

	return D;
}

template <size_t _size>
inline double det(const SymmetricMatrix<_size>& q) { 
	Matrix<_size, _size> m = q;
	double D = double(1);

	// Gaussian elimination
	for (size_t k = 0; k < _size; ++k) {
		D *= m(k, k);
		if (D == double(0)) {
			return double(0);
		}

		// eliminate column
		for (size_t i = k + 1; i < _size; ++i) {
			double factor = m(i, k) / m(k, k);
			for (size_t j = k + 1; j < _size; ++j) {
				m(i, j) -= factor * m(k, j);
			}
		}  
	}

	return D;
}


// P%Q solves PX = Q for X
template <size_t _size, size_t _numColumns>
inline Matrix<_size, _numColumns> operator%(const Matrix<_size, _size>& p, const Matrix<_size, _numColumns>& q) {
	
	Matrix<_size, _size> m(p);
	Matrix<_size, _numColumns> inv(q);

	size_t row_p[_size];
	size_t col_p[_size];
	for (size_t i = 0; i < _size; ++i) {
		row_p[i] = i; col_p[i] = i;
	}

	// Gaussian elimination
	for (size_t k = 0; k < _size; ++k) {
		// find maximal pivot element
		double maximum = double(0); size_t max_row = k; size_t max_col = k;
		for (size_t i = k; i < _size; ++i) {
			for (size_t j = k; j < _size; ++j) {
				double abs_ij = abs(m(row_p[i], col_p[j]));
				if (abs_ij > maximum) {
					maximum = abs_ij; max_row = i; max_col = j;
				}
			}
		}

		assert(maximum != double(0));

		// swap rows and columns
		std::swap(row_p[k], row_p[max_row]);
		std::swap(col_p[k], col_p[max_col]);

		// eliminate column
		for (size_t i = k + 1; i < _size; ++i) {
			double factor = m(row_p[i], col_p[k]) / m(row_p[k], col_p[k]);
			for (size_t j = k + 1; j < _size; ++j) {
				m(row_p[i], col_p[j]) -= factor * m(row_p[k], col_p[j]);
			}
			for (size_t j = 0; j < _numColumns; ++j) {
				inv(row_p[i], j) -= factor * inv(row_p[k], j);
			}
		} 
	}

	// Backward substitution
	for (size_t k = _size - 1; k != -1; --k) {
		double quotient = m(row_p[k], col_p[k]);

		for (size_t j = 0; j < _numColumns; ++j) {
			inv(row_p[k], j) /= quotient;
		}

		for (size_t i = 0; i < k; ++i) {
			double factor = m(row_p[i], col_p[k]);
			for (size_t j = 0; j < _numColumns; ++j) {
				inv(row_p[i], j) -= factor * inv(row_p[k], j);
			}
		}
	}

	// reshuffle result
	size_t invrow_p[_size];
	for (size_t i = 0; i < _size; ++i) {
		invrow_p[row_p[i]] = i;
	}
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _numColumns; ++j) {
			std::swap(inv(col_p[i], j), inv(row_p[i], j));
		}
		row_p[invrow_p[col_p[i]]] = row_p[i];
		invrow_p[row_p[i]] = invrow_p[col_p[i]];
	}

	return inv;
}

template <size_t _size, size_t _numColumns>
inline Matrix<_size, _numColumns> operator%(const SymmetricMatrix<_size>& p, const Matrix<_size, _numColumns>& q) {
	// Cholesky factorization p = L*~L
	SymmetricMatrix<_size> L; // abuse SymmetricMatrix for triangular matrix
	L.reset();
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = i; j < _size; ++j) {
			double sum = p(j,i);
			for (size_t k = 0; k < i; ++k) {
				sum -= L(j,k)*L(i,k);
			}
			if (i == j) {
				assert(sum > 0.0);
				L(i,i) = sqrt(sum);
			} else {
				L(j,i) = sum / L(i,i);
			}
		}
	}

	// Backward and forward substitution
	Matrix<_size, _numColumns> M;
	for (size_t i = 0; i < _size; ++i) {
		for (size_t k = 0; k < _numColumns; ++k) {
			double sum = q(i,k);
			for (size_t j = 0; j < i; ++j) {
				sum -= L(i,j)*M(j,k);
			}
			assert(L(i,i) != 0);
			M(i,k) = sum / L(i,i);
		}
	}
	for (size_t i = _size - 1; i != -1; --i) {
		for (size_t k = 0; k < _numColumns; ++k) {
			double sum = M(i,k);
			for (size_t j = i + 1; j < _size; ++j) {
				sum -= L(j,i)*M(j,k);
			}
	
			assert(L(i,i) != 0);
			M(i,k) = sum / L(i,i);
		}
	}
	return M;
}


template <size_t _size>
inline Matrix<_size, _size> operator%(const SymmetricMatrix<_size>& p, const SymmetricMatrix<_size>& q) {
	
	// Cholesky factorization p = L*~L
	SymmetricMatrix<_size> L; // abuse SymmetricMatrix for triangular matrix
	L.reset();
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = i; j < _size; ++j) {
			double sum = p(j,i);
			for (size_t k = 0; k < i; ++k) {
				sum -= L(j,k)*L(i,k);
			}
			if (i == j) {
				assert(sum > 0.0);
				L(i,i) = sqrt(sum);
			} else {
				L(j,i) = sum / L(i,i);
			}
		}
	}

	// Backward and forward substitution
	Matrix<_size, _size> M;
	for (size_t i = 0; i < _size; ++i) {
		for (size_t k = 0; k < _size; ++k) {
			double sum = q(i,k);
			for (size_t j = 0; j < i; ++j) {
				sum -= L(i,j)*M(j,k);
			}
			M(i,k) = sum / L(i,i);
		}
	}
	for (size_t i = _size - 1; i != -1; --i) {
		for (size_t k = 0; k < _size; ++k) {
			double sum = M(i,k);
			for (size_t j = i + 1; j < _size; ++j) {
				sum -= L(j,i)*M(j,k);
			}
			M(i,k) = sum / L(i,i);
		}
	}
	return M;
}


template <size_t _size>
inline SymmetricMatrix<_size> operator!(const SymmetricMatrix<_size>& p) { 
	return 0.5*SymSum(p % (Matrix<_size, _size>) identity<_size>() );
}

template <size_t _size, size_t _numRows>
inline Matrix<_numRows, _size> operator/(const Matrix<_numRows, _size>& p, const Matrix<_size, _size>& q) {
	return ~(~q%~p);
}

template <size_t _size, size_t _numRows>
inline Matrix<_numRows, _size> operator/(const Matrix<_numRows, _size>& p, const SymmetricMatrix<_size>& q) {
	return ~(q%~p);
}

// Matrix pseudo-inverse
template <size_t _numRows, size_t _numColumns>
inline Matrix<_numColumns, _numRows> pseudoInverse(const Matrix<_numRows, _numColumns>& q) { 
	if (_numColumns <= _numRows) {
		Matrix<_numColumns, _numColumns> Vec;
		SymmetricMatrix<_numColumns> Val;
		jacobi(SymProd(~q,q), Vec, Val);

		for (size_t i = 0; i < _numColumns; ++i) {
			if (abs(Val(i,i)) <= sqrt(DBL_EPSILON)) {
				Val(i,i) = 0.0;
			} else {
				Val(i,i) = 1.0 / Val(i,i);
			}
		}
		return SymProd(Vec,Val*~Vec)*~q;
	} else {
		Matrix<_numRows, _numRows> Vec;
		SymmetricMatrix<_numRows> Val;
		jacobi(SymProd(q,~q), Vec, Val);

		for (size_t i = 0; i < _numRows; ++i) {
			if (abs(Val(i,i)) <= sqrt(DBL_EPSILON)) {
				Val(i,i) = 0.0;
			} else {
				Val(i,i) = 1.0 / Val(i,i);
			}
		}
		return ~q*SymProd(Vec,Val*~Vec);
	}
}

template <size_t _size>
inline SymmetricMatrix<_size> pseudoInverse(const SymmetricMatrix<_size>& q) { 
	Matrix<_size, _size> Vec;
	SymmetricMatrix<_size> Val;
	jacobi(q, Vec, Val);

	for (size_t i = 0; i < _size; ++i) {
		if (abs(Val(i,i)) <= sqrt(DBL_EPSILON)) {
			Val(i,i) = 0.0; //1000000.0;
		} else {
			Val(i,i) = 1.0 / Val(i,i);
		}
	}
	return SymProd(Vec,Val*~Vec);
}

// Matrix inverse
template <size_t _size>
inline Matrix<_size, _size> operator!(const Matrix<_size, _size>& q) { 
	Matrix<_size, _size> m(q);
	Matrix<_size, _size> inv = identity<_size>();

	size_t row_p[_size];
	size_t col_p[_size];
	for (size_t i = 0; i < _size; ++i) {
		row_p[i] = i; col_p[i] = i;
	}

	// Gaussian elimination
	for (size_t k = 0; k < _size; ++k) {
		// find maximal pivot element
		double maximum = double(0); size_t max_row = k; size_t max_col = k;
		for (size_t i = k; i < _size; ++i) {
			for (size_t j = k; j < _size; ++j) {
				double abs_ij = abs(m(row_p[i], col_p[j]));
				if (abs_ij > maximum) {
					maximum = abs_ij; max_row = i; max_col = j;
				}
			}
		}

		// swap rows and columns
		size_t swap = row_p[k]; row_p[k] = row_p[max_row]; row_p[max_row] = swap;
		swap = col_p[k]; col_p[k] = col_p[max_col]; col_p[max_col] = swap;

		// eliminate column
		assert(maximum != double(0));
		for (size_t i = k + 1; i < _size; ++i) {
			double factor = m(row_p[i], col_p[k]) / m(row_p[k], col_p[k]);

			for (size_t j = k + 1; j < _size; ++j) {
				m(row_p[i], col_p[j]) -= factor * m(row_p[k], col_p[j]);
			}

			for (size_t j = 0; j < k; ++j) {
				inv(row_p[i], row_p[j]) -= factor * inv(row_p[k], row_p[j]);
			}
			inv(row_p[i], row_p[k]) = -factor;
		} 
	}

	// Backward substitution
	for (size_t k = _size - 1; k != -1; --k) {
		double quotient = m(row_p[k], col_p[k]);

		for (size_t j = 0; j < _size; ++j) {
			inv(row_p[k], j) /= quotient;
		}

		for (size_t i = 0; i < k; ++i) {
			double factor = m(row_p[i], col_p[k]);
			for (size_t j = 0; j < _size; ++j) {
				inv(row_p[i], j) -= factor * inv(row_p[k], j);
			}
		}
	}

	// reshuffle result
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			m(col_p[i], j) = inv(row_p[i], j);
		}
	}

	return m; 
}

// Matrix exponentiation
#define _MATRIX_B0 1729728e1
#define _MATRIX_B1 864864e1
#define _MATRIX_B2 199584e1
#define _MATRIX_B3 2772e2
#define _MATRIX_B4 252e2
#define _MATRIX_B5 1512e0
#define _MATRIX_B6 56e0
#define _MATRIX_B7 1e0
#define _NORMLIM 9.504178996162932e-1

template <size_t _size>
inline Matrix<_size, _size> exp(const Matrix<_size, _size>& q) { 
	Matrix<_size, _size> A(q);
	int s = (int) std::max(double(0), ceil(log(norm(A)/_NORMLIM)*M_LOG2E)); 

	A /= pow(2.0,s);
	Matrix<_size, _size> A2(A*A);
	Matrix<_size, _size> A4(A2*A2);
	Matrix<_size, _size> A6(A2*A4);
	Matrix<_size, _size> U( A*(A6*_MATRIX_B7 + A4*_MATRIX_B5 + A2*_MATRIX_B3 + identity<_size>()*_MATRIX_B1) );
	Matrix<_size, _size> V( A6*_MATRIX_B6 + A4*_MATRIX_B4 + A2*_MATRIX_B2 + identity<_size>()*_MATRIX_B0 ); 
	Matrix<_size, _size> R7 = (V - U)%(V + U);

	for (int i = 0; i < s; ++i) {
		R7 = R7*R7;
	}
	return R7;
}

// Input stream
template <size_t _numRows, size_t _numColumns>
inline std::istream& operator>>(std::istream& is, Matrix<_numRows, _numColumns>& q) {
	for (size_t i = 0; i < _numRows; ++i) {
		for (size_t j = 0; j < _numColumns; ++j) {
			is >> q(i,j);
		}
	}
	return is;
}

// Output stream
template <size_t _numRows, size_t _numColumns>
inline std::ostream& operator<<(std::ostream& os, const Matrix<_numRows, _numColumns>& q) {
	for (size_t i = 0; i < _numRows; ++i) {
		for (size_t j = 0; j < _numColumns; ++j) {
			os  << std::left << std::setw(13) << q(i,j) << " ";
		}
		os << std::endl;
	}
	return os;
}

template <size_t _size>
inline std::ostream& operator<<(std::ostream& os, const SymmetricMatrix<_size>& q) {
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			os << std::left << std::setw(13) << q(i,j) << " ";
		}
		os << std::endl;
	}
	return os;
}


// Hilbert matrix
template <size_t _size>
inline Matrix<_size, _size> hilbert() {
	Matrix<_size, _size> m;
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			m(i, j) = double(1) / (double) (i + j + 1);
		}
	}
	return m;
}

// Cholesky decomposition
template <size_t _size>
inline Matrix<_size, _size> chol(const SymmetricMatrix<_size>& p)
{
	Matrix<_size, _size> L; // abuse SymmetricMatrix for triangular matrix
	L.reset();
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = i; j < _size; ++j) {
			double sum = p(j,i);
			for (size_t k = 0; k < i; ++k) {
				sum -= L(j,k)*L(i,k);
			}
			if (i == j) {
				assert(sum > 0.0);
				L(i,i) = sqrt(sum);
			} else {
				L(j,i) = sum / L(i,i);
			}
		}
	}
	return L;
}

template <size_t _size>
inline void jacobi(const SymmetricMatrix<_size>& m, Matrix<_size, _size>& V, SymmetricMatrix<_size>& D) {
	D = m;
	V = (Matrix<_size,_size>) identity<_size>();

	if (_size <= 1) {
		return;
	}

	size_t pivotRow = 0;
	size_t zeroCount = 0;

	while (true) {
		double maximum = 0;
		size_t p, q;
		for (size_t i = 0; i < pivotRow; ++i) {
			if (abs(D(i,pivotRow)) > maximum) {
				maximum = abs(D(i,pivotRow));
				p = i; q = pivotRow;
			}
		}
		for (size_t j = pivotRow + 1; j < _size; ++j) {
			if (abs(D(pivotRow,j)) > maximum) {
				maximum = abs(D(pivotRow,j));
				p = pivotRow; q = j;
			}
		}

		pivotRow = (pivotRow + 1) % _size;

		if (maximum <= DBL_EPSILON) {
			++zeroCount;
			if (zeroCount == _size) {
				break;
			} else {
				continue;
			}
		} else {
			zeroCount = 0;
		}

		double theta = 0.5*(D(q,q) - D(p,p)) / D(p,q);
		double t = 1 / (abs(theta) + _hypot(theta, 1));
		if (theta < 0) t = -t;
		double c = 1 / _hypot(t, 1);
		double s = c*t;
		double tau = s / (1 + c);

		// update D // 

		// update D(r,p) and D(r,q)
		for (size_t r = 0; r < p; ++r) {
			double Drp = D(r,p); double Drq = D(r,q);
			D(r,p) -= s*(Drq + tau*Drp);
			D(r,q) += s*(Drp - tau*Drq);
		}
		for (size_t r = p + 1; r < q; ++r) {
			double Drp = D(p,r); double Drq = D(r,q);
			D(p,r) -= s*(Drq + tau*Drp);
			D(r,q) += s*(Drp - tau*Drq);
		}
		for (size_t r = q + 1; r < _size; ++r) {
			double Drp = D(p,r); double Drq = D(q,r);
			D(p,r) -= s*(Drq + tau*Drp);
			D(q,r) += s*(Drp - tau*Drq);
		}

		// update D(p,p), D(q,q), and D(p,q);
		D(p,p) -= t*D(p,q);
		D(q,q) += t*D(p,q);
		D(p,q) = 0;

		// Update V
		for (size_t r = 0; r < _size; ++r) {
			double Vrp = V(r,p); double Vrq = V(r,q);
			V(r,p) -= s*(Vrq + tau*Vrp);
			V(r,q) += s*(Vrp - tau*Vrq);
		}
	}

	// clean up D
	for (size_t i = 0; i < _size - 1; ++i) {
		for (size_t j = 0; j < i; ++j) {
			D(i,j) = 0;
		}
	}
}




// Matrix addition
template <size_t _numRows, size_t _numColumns>
inline Matrix<_numRows, _numColumns> operator+(const Matrix<_numRows, _numColumns>& M, const Matrix<_numRows, _numColumns>& N) { 
	Matrix<_numRows, _numColumns> L;
	for (size_t i = 0; i < _numRows * _numColumns; ++i) {
		L[i] = M[i] + N[i];
	}
	return L;
}

template <size_t _size>
inline Matrix<_size, _size> operator+(const SymmetricMatrix<_size>& M, const Matrix<_size, _size>& N) { 
	Matrix<_size, _size> L;
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			L(i,j) = M(i,j) + N(i,j);
		}
	}
	return L;
}

template <size_t _size>
inline Matrix<_size, _size> operator+(const Matrix<_size, _size>& M, const SymmetricMatrix<_size>& N) { 
	Matrix<_size, _size> L;
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			L(i,j) = M(i,j) + N(i,j);
		}
	}
	return L;
}

template <size_t _size>
inline SymmetricMatrix<_size> operator+(const SymmetricMatrix<_size>& M, const SymmetricMatrix<_size>& N) { 
	SymmetricMatrix<_size> L;
	for (size_t i = 0; i < ((_size+1)*_size)/2; ++i) {
		L[i] = M[i] + N[i];
	}
	return L;
}

// Matrix subtraction
template <size_t _numRows, size_t _numColumns>
inline Matrix<_numRows, _numColumns> operator-(const Matrix<_numRows, _numColumns>& M, const Matrix<_numRows, _numColumns>& N) { 
	Matrix<_numRows, _numColumns> L;
	for (size_t i = 0; i < _numRows * _numColumns; ++i) {
		L[i] = M[i] - N[i];
	}
	return L;
}

template <size_t _size>
inline Matrix<_size, _size> operator-(const SymmetricMatrix<_size>& M, const Matrix<_size, _size>& N) { 
	Matrix<_size, _size> L;
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			L(i,j) = M(i,j) - N(i,j);
		}
	}
	return L;
}

template <size_t _size>
inline Matrix<_size, _size> operator-(const Matrix<_size, _size>& M, const SymmetricMatrix<_size>& N) { 
	Matrix<_size, _size> L;
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			L(i,j) = M(i,j) - N(i,j);
		}
	}
	return L;
}

template <size_t _size>
inline SymmetricMatrix<_size> operator-(const SymmetricMatrix<_size>& M, const SymmetricMatrix<_size>& N) { 
	SymmetricMatrix<_size> L;
	for (size_t i = 0; i < ((_size+1)*_size)/2; ++i) {
		L[i] = M[i] - N[i];
	}
	return L;
}

// Scalar multiplication
template <size_t _numRows, size_t _numColumns>
inline Matrix<_numRows, _numColumns> operator*(const Matrix<_numRows, _numColumns>& M, double a) { 
	Matrix<_numRows, _numColumns> L;
	for (size_t i = 0; i < _numRows * _numColumns; ++i) {
		L[i] = M[i] * a;
	}
	return L;
}

template <size_t _numRows, size_t _numColumns>
inline Matrix<_numRows, _numColumns> operator*(double a, const Matrix<_numRows, _numColumns>& M) { 
	Matrix<_numRows, _numColumns> L;
	for (size_t i = 0; i < _numRows * _numColumns; ++i) {
		L[i] = a * M[i];
	}
	return L;
}

template <size_t _size>
inline SymmetricMatrix<_size> operator*(const SymmetricMatrix<_size>& M, double a) { 
	SymmetricMatrix<_size> L;
	for (size_t i = 0; i < ((_size+1)*_size)/2; ++i) {
		L[i] = M[i] * a;
	}
	return L;
}

template <size_t _size>
inline SymmetricMatrix<_size> operator*(double a, const SymmetricMatrix<_size>& M) { 
	SymmetricMatrix<_size> L;
	for (size_t i = 0; i < ((_size+1)*_size)/2; ++i) {
		L[i] = a * M[i];
	}
	return L;
}

// Scalar division
template <size_t _numRows, size_t _numColumns>
inline Matrix<_numRows, _numColumns> operator/(const Matrix<_numRows, _numColumns>& M, double a) { 
	Matrix<_numRows, _numColumns> L;
	for (size_t i = 0; i < _numRows * _numColumns; ++i) {
		L[i] = M[i] / a;
	}
	return L;
}

template <size_t _size>
inline SymmetricMatrix<_size> operator/(const SymmetricMatrix<_size>& M, double a) { 
	SymmetricMatrix<_size> L;
	for (size_t i = 0; i < ((_size+1)*_size)/2; ++i) {
		L[i] = M[i] / a;
	}
	return L;
}

// Matrix multiplication
template <size_t _numRows, size_t _numColumns, size_t __numColumns>
inline Matrix<_numRows, __numColumns> operator*(const Matrix<_numRows, _numColumns>& M, const Matrix<_numColumns, __numColumns>& N) {
	Matrix<_numRows, __numColumns> L;
	for (size_t i = 0; i < _numRows; ++i) {
		for (size_t j = 0; j < __numColumns; ++j) {
			double temp = double(0);
			for (size_t k = 0; k < _numColumns; ++k) {
				temp += M(i, k) * N(k, j);
			}
			L(i, j) = temp;
		}
	}
	return L;
}

template <size_t _size, size_t _numColumns>
inline Matrix<_size, _numColumns> operator*(const SymmetricMatrix<_size>& M, const Matrix<_size, _numColumns>& N) {
	Matrix<_size, _numColumns> L;
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _numColumns; ++j) {
			double temp = double(0);
			for (size_t k = 0; k < _size; ++k) {
				temp += M(i, k) * N(k, j);
			}
			L(i, j) = temp;
		}
	}
	return L;
}

template <size_t _size, size_t _numRows>
inline Matrix<_numRows, _size> operator*(const Matrix<_numRows, _size>& M, const SymmetricMatrix<_size>& N) {
	Matrix<_numRows, _size> L;
	for (size_t i = 0; i < _numRows; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			double temp = double(0);
			for (size_t k = 0; k < _size; ++k) {
				temp += M(i, k) * N(k, j);
			}
			L(i, j) = temp;
		}
	}
	return L;
}

template <size_t _size>
inline Matrix<_size, _size> operator*(const SymmetricMatrix<_size>& M, const SymmetricMatrix<_size>& N) {
	Matrix<_size, _size> L;
	for (size_t i = 0; i < _size; ++i) {
		for (size_t j = 0; j < _size; ++j) {
			double temp = double(0);
			for (size_t k = 0; k < _size; ++k) {
				temp += M(i, k) * N(k, j);
			}
			L(i, j) = temp;
		}
	}
	return L;
}

// Compute product M*N*~M where N is symmetric
template<size_t _size, size_t _numRows>
inline SymmetricMatrix<_numRows> sqrForm(const Matrix<_numRows, _size>& M, const SymmetricMatrix<_size>& N) {
	return SymProd(M*N,~M);
}

// Compute product M*~M
template<size_t _numRows, size_t _numColumns>
inline SymmetricMatrix<_numRows> sqrForm(const Matrix<_numRows, _numColumns>& M) {
	return SymProd(M,~M);
}

// Compute product M*N of which one knows that the result is symmetric (and save half the computation)
template <size_t _size, size_t _numRows>
inline SymmetricMatrix<_size> SymProd(const Matrix<_size, _numRows>& M, const Matrix<_numRows, _size>& N) {
	SymmetricMatrix<_size> S;
	for (size_t j = 0; j < _size; ++j) {
		for (size_t i = j; i < _size; ++i) {
			double temp = double(0);
			for (size_t k = 0; k < _numRows; ++k) {
				temp += M(i, k) * N(k, j);
			}
			S(i, j) = temp;
		}
	}
	return S;
}

template <size_t _size>
inline SymmetricMatrix<_size> SymProd(const SymmetricMatrix<_size>& M, const SymmetricMatrix<_size>& N) {
	SymmetricMatrix<_size> S;
	for (size_t j = 0; j < _size; ++j) {
		for (size_t i = j; i < _size; ++i) {
			double temp = double(0);
			for (size_t k = 0; k < _size; ++k) {
				temp += M(i, k) * N(k, j);
			}
			S(i, j) = temp;
		}
	}
	return S;
}

template <size_t _size>
inline SymmetricMatrix<_size> SymProd(const SymmetricMatrix<_size>& M, const Matrix<_size, _size>& N) {
	SymmetricMatrix<_size> S;
	for (size_t j = 0; j < _size; ++j) {
		for (size_t i = j; i < _size; ++i) {
			double temp = double(0);
			for (size_t k = 0; k < _size; ++k) {
				temp += M(i, k) * N(k, j);
			}
			S(i, j) = temp;
		}
	}
	return S;
}

template <size_t _size>
inline SymmetricMatrix<_size> SymProd(const Matrix<_size, _size>& M, const SymmetricMatrix<_size>& N) {
	SymmetricMatrix<_size> S;
	for (size_t j = 0; j < _size; ++j) {
		for (size_t i = j; i < _size; ++i) {
			double temp = double(0);
			for (size_t k = 0; k < _size; ++k) {
				temp += M(i, k) * N(k, j);
			}
			S(i, j) = temp;
		}
	}
	return S;
}


// Compute sum M+M^T
template <size_t _size>
inline SymmetricMatrix<_size> SymSum(const Matrix<_size, _size>& M) {
	SymmetricMatrix<_size> S;
	for (size_t j = 0; j < _size; ++j) {
		for (size_t i = j; i < _size; ++i) {
			S(i,j) = M(i,j) + M(j,i);
		}
	}
	return S;
}

// Principal square root
template <size_t _size>
inline SymmetricMatrix<_size> sqrt(const SymmetricMatrix<_size>& M) {
	Matrix<_size,_size> V;
	SymmetricMatrix<_size> D;
	jacobi(M, V, D);
	for (size_t i = 0; i < _size; ++i) {
		if (D(i,i) > 0) {
			D(i,i) = sqrt(D(i,i));
		} else {
			D(i,i) = 0;
		}
	}
	return SymProd(V, D*~V);
}

inline double scalar(const Matrix<1,1>& M) {
	return M[0];
}


#endif