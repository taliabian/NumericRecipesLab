/*!
* \file  matrix.h
* \brief Implementation for NumericRecipesLab class
* \author Talia
* \version 1.1
* \date 2015-04-10 
*/

template <typename Type>
Numeric<Type>::Numeric()
{
	//
}

template <typename Type>
Numeric<Type>::~Numeric()
{
	//
}

/// General Gaussian elimination
template <typename Type>
bool Numeric<Type>:: GenGaussElimation( const Matrix<Type> &A, const vector<Type> &b)
{
	int num = A.rows();
	Matrix<Type> tmpA = A;
	Matrix<Type> tmpB = A;
	vector<Type> tmpVa = b;
	vector<Type> tmpVb = b;
	//vector<Type> tmpX(num);
	vX = vector<Type>( num );
	int i,j;
	for ( int m = 0; m<(num-1); m++)
	{
		for ( i=(m+1); i<num; i++)
		{
			for( j=m; j<num; j++)
			{
				if( j == m )
					tmpB[i][j] = 0;
				else
				{
					if ( tmpA[m][m] == 0)
						return false;
					else
						tmpB[i][j] = tmpB[i][j] - tmpA[m][j]*tmpA[i][m]/tmpA[m][m];
				}
			}
			tmpVb[i] = tmpVb[i] - tmpVa[m]*tmpA[i][m] /tmpA[m][m];
		}
		tmpA = tmpB;
		tmpVa = tmpVb;
	}
	for ( i=(num-1); i>=0; i--)
	{
		if ( i==(num-1) )
		{
			vX[num-1] = tmpVb[num-1];
		}
		else
		{
			vX[i] = tmpVb[i];
			for ( j=(num-1); j>i; j--)
			{
				vX[i] =  vX[i] - tmpB[i][j]*vX[j]; 
			}
		}
		vX[i] = vX[i]/tmpB[i][i];
	}
	return true;
}

/// Gaussian elimination with full pivoting
template <typename Type>
bool FullPivotGaussElimation( const Matrix<Type> &A, const vector<Type> &b )
{

}

///  get the the result of equations of m by m matrix
template <typename Type>
vector<Type> Numeric<Type>::getvX() const
{
	return vX;
}

/// \brief 重载"=",系数赋值
/// \param x 右值
/// \return 赋值左参数
template <typename Type>
Type Numeric<Type>::operator=( const Type &x )
{
	return tmp;
}