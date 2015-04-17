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
		Numeric( const Matrix<Type> &M );
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
		/// L-U decomposition of Matrix 
		bool MatLUdec( const Matrix<Type> &A);
		/// L-U-P decomposition of Matrix
		bool MatLUPdec( const Matrix<Type> &A);
		/// solve the multiply funciton with L-U-P 
		bool LUPsolveFun( const Matrix<Type> &A, const vector<Type> &b);
		/// solve the multiply funciton with LeastSquares( Overdetermined  Functions ) 
		bool LeasetSquaresSolveFun( const Matrix<Type> &A, const vector<Type> &b);
		/// inv of L
		bool InvL( Matrix<Type> &L);
		/// inv of U;
		bool InvU( Matrix<Type> &U);
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
		/// get the the result of a matrix's L's inv
		Matrix<Type> getMatinvL() const; 
		/// get the the result of a matrix's U's inv
		Matrix<Type> getMatinvU() const; 
		///template<typename Type>
		Matrix<Type> PerMatrix( int size, int *pcol ); 		
		Type operator=( const Type &x );

	private:
		/// 
		Matrix<Type> M;
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
		/// the result of matrix's L's inv
		Matrix<Type> invLM;
		/// the result of matrix's U's inv
		Matrix<Type> invUM;

	};
	/// produce a Permutation matrix with 1's location
	

	#include "numeric_impl.h"

}

#endif