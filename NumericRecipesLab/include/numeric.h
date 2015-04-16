/*!
* \file  numeric.h
* \brief Class template of NumericRecipesLab which is designed for basic Numerical Analysis
* \author Talia
* \version 1.1
* \date 2015-04-10 
*/

#ifndef NUMERIC_H
#define NUMERIC_H


#include <vector>
#include <complex>
#include "matrix.h"

using namespace std;
using namespace matrixlab;

namespace NumericRecipesLab
{

	template <typename Type>
	class Numeric
	{
	public:
		/// constructor and deconstructor
		Numeric();
		~Numeric();

		/// General Gaussian elimination
		bool GenGaussElimation( const Matrix<Type> &A, const vector<Type> &b );
		/// Gaussian elimination with full pivoting
		bool FullPivotGaussElimation( const Matrix<Type> &A, const vector<Type> &b );
		/// Gaussian elimination with partital column pivoting
		bool PartialPivotGaussElimation( const Matrix<Type> &A, const vector<Type> &b );
		/// Gaussian-Jordan elimination 
		bool GaussJordanElimation( const Matrix<Type> &A, const vector<Type> &b );
		/// Compute symbolic matrix inverse with Gaussian-Jordan Elimintaion
		bool InvMwithGaussJordan( const Matrix<Type> &A);
		/// L-U decomposition of Matrix of L
		bool MatLUdec( const Matrix<Type> &A);
		/// L-U decomposition of Matrix of U
		bool MatLUPdec( const Matrix<Type> &A);
		/// get the the result of equations of m by m matrix
		vector<Type> getvX() const ;
		/// get the the result of a matrix's inv
		Matrix<Type> getinvM() const; 
		/// get the the result of a matrix's L
		Matrix<Type> getMatL() const; 
		/// get the the result of a matrix's U
		Matrix<Type> getMatU() const; 
		/// get the the result of a matrix's P
		Matrix<Type> getMatP() const; 
		/// 		
		Type operator=( const Type &x );

	private:
		/// the result of equations of m by m matrix
		vector<Type> vX;
		/// the result of matrix's inv
		Matrix<Type> invM;
		/// the result of matrix's L
		Matrix<Type> LM;
		/// the result of matrix's U
		Matrix<Type> UM;
		/// the result of matrix's P
		Matrix<Type> PM;

	};
	#include "numeric_impl.h"

}

#endif