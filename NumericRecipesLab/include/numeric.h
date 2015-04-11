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
		/// get the the result of equations of m by m matrix
		vector<Type> getvX() const ;

		/// 		
		Type operator=( const Type &x );

	private:
		/// the result of equations of m by m matrix
		vector<Type> vX;

	};
	#include "numeric_impl.h"

}

#endif