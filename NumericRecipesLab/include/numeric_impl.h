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
	vX = vector<Type>( num );
	int i,j;
	for ( int m = 0; m<(num-1); m++)
	{
		for ( i=(m+1); i<num; i++)
		{
			for( j=m; j<num; j++)
			{
				if( j == m )
				{
					if( tmpB[i][j] == 0)
						break;
					else
						tmpB[i][j] = 0;
				}
				else
				{
					if ( tmpA[m][m] == 0)
						return false;
					else
						tmpB[i][j] = tmpB[i][j] - tmpA[m][j]*tmpA[i][m]/tmpA[m][m];
				}
			}
			if(  tmpA[i][m] == 0 )
				continue;
			else
				tmpVb[i] = tmpVb[i] - tmpVa[m]*tmpA[i][m] /tmpA[m][m];// i's row
		}
		tmpA = tmpB;
		tmpVa = tmpVb;
	}
	// compute the x[i]
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
bool Numeric<Type>::FullPivotGaussElimation( const Matrix<Type> &A, const vector<Type> &b )
{
	
	int num = A.rows();
	Matrix<Type> tmpA = A;
	Matrix<Type> tmpB;
	Matrix<Type> tmpC;
	vector<Type> tmpVa = b;
	vector<Type> tmpVb = b;
	vX = vector<Type>( num );
	int *maxrow = new int[num-1];
	int *maxcol = new int[num-1];
	Data_Pos<Type> maxandpos;
	int i,j;
	Type tmp;
	for ( int m = 0; m<(num-1); m++) // the m'th cycle
	{
		/// find the max pivoting element of matrix
		tmpC = CopyFromMatrix(tmpA, m, m, num-1, num-1);
		maxandpos = FindMaxandPos(abs(tmpC));
		maxrow[m] = maxandpos.row + m;
		maxcol[m] = maxandpos.col + m;
		/// change the rows and then cols 
		tmpA = ExchangeRows(tmpA,m,maxrow[m]);
		tmpA = ExchangeCols(tmpA,m,maxcol[m]);
		/// change the right constant
		Type tmp = tmpVa.at(m);
		tmpVa.at(m) = tmpVa.at(maxrow[m]);
		tmpVa.at(maxrow[m]) = tmp;
		tmpB = tmpA;
		tmpVb = tmpVa;
		for ( i=(m+1); i<num; i++)
		{
			for( j=m; j<num; j++)
			{
				if( j == m )
				{
					if( tmpB[i][j] == 0)
						break;
					else
						tmpB[i][j] = 0;
				}
				else
				{
					if ( tmpA[m][m] == 0)
						return false;
					else
						tmpB[i][j] = tmpB[i][j] - tmpA[m][j]*tmpA[i][m]/tmpA[m][m];
				}
			}
			if(  tmpA[i][m] == 0 )
				continue;
			else
				tmpVb[i] = tmpVb[i] - tmpVa[m]*tmpA[i][m] /tmpA[m][m];// i's row
		}
		tmpA = tmpB;
		tmpVa = tmpVb;
	}
	// compute the x[i]
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
	// back the x[i] origin squencese
	for( i=num-2; i>=0; i--)
	{
		tmp = vX.at(i);
		vX.at(i) = vX.at(maxcol[i]);
		vX.at(maxcol[i]) = tmp;
	}
	delete maxrow;
	delete maxcol;
	return true;
}

/// Gaussian elimination with partital column pivoting
template <typename Type>
bool Numeric<Type>::PartialPivotGaussElimation( const Matrix<Type> &A, const vector<Type> &b )
{
	int num = A.rows();
	Matrix<Type> tmpA = A;
	Matrix<Type> tmpB;
	Matrix<Type> tmpC;
	vector<Type> tmpVa = b;
	vector<Type> tmpVb = b;

	vX = vector<Type>( num );
	int *maxrow = new int[num-1];
	Data_Pos<Type> maxandpos;
	int i,j;
	for ( int m = 0; m<(num-1); m++) // the m'th cycle
	{
		/// find the max pivoting element of matrix's m column
		tmpC = CopyFromMatrix(tmpA, m, m, num-1, m);
		maxandpos = FindMaxandPos(abs(tmpC));
		maxrow[m] = maxandpos.row + m;
		/// change the rows 
		tmpA = ExchangeRows(tmpA,m,maxrow[m]);
		/// change the right constant
		Type tmp = tmpVa.at(m);
		tmpVa.at(m) = tmpVa.at(maxrow[m]);
		tmpVa.at(maxrow[m]) = tmp;
		tmpB = tmpA;
		tmpVb = tmpVa;
		for ( i=(m+1); i<num; i++)
		{
			for( j=m; j<num; j++)
			{
				if( j == m )
				{
					if( tmpB[i][j] == 0)
						break;
					else
						tmpB[i][j] = 0;
				}
				else
				{
					if ( tmpA[m][m] == 0)
						return false;
					else
						tmpB[i][j] = tmpB[i][j] - tmpA[m][j]*tmpA[i][m]/tmpA[m][m];
				}
			}
			if(  tmpA[i][m] == 0 )
				continue;
			else
				tmpVb[i] = tmpVb[i] - tmpVa[m]*tmpA[i][m] /tmpA[m][m];// i's row
		}
		tmpA = tmpB;
		tmpVa = tmpVb;
	}
	/// compute the x[i]
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
	
	delete maxrow;
	return true;
}
/// Gaussian-Jordan elimination 
template <typename Type>
bool  Numeric<Type>::GaussJordanElimation( const Matrix<Type> &A, const vector<Type> &b )
{
	int num = A.rows();
	Matrix<Type> tmpA = A;
	Matrix<Type> tmpB;
	Matrix<Type> tmpC;
	Matrix<Type> tmpD = eye(num,Type(1));
	Matrix<Type> tmpE = eye(num, Type(1));
	vector<Type> tmpVa = b;
	vector<Type> tmpVb = b;

	vX = vector<Type>( num );
	int *maxrow = new int[num-1];
	Data_Pos<Type> maxandpos;
	int i,j,k;
	for ( int m = 0; m<(num-1); m++) // the m'th cycle
	{
		/// find the max pivoting element of matrix's m column
		tmpC = CopyFromMatrix(tmpA, m, m, num-1, m);
		maxandpos = FindMaxandPos(abs(tmpC));
		maxrow[m] = maxandpos.row + m;
		/// change the rows 
		tmpA = ExchangeRows(tmpA,m,maxrow[m]);
		/// change the right constant
		Type tmp = tmpVa.at(m);
		tmpVa.at(m) = tmpVa.at(maxrow[m]);
		tmpVa.at(maxrow[m]) = tmp;
		tmpB = tmpA;
		tmpVb = tmpVa;
		/// change the tmpD( for invA ) squence
		tmpD = ExchangeRows(tmpD,m,maxrow[m]);
		tmpE = tmpD;
		for ( i=0; i<num; i++)// row
		{
			if( i != m )
			{
				for( j=m; j<num; j++)// column
				{
					if( j == m )
					{
						if( tmpB[i][j] == 0)
							break;
						else
							tmpB[i][j] = 0;
					}
					else
					{
						if ( tmpA[m][m] == 0)
							return false;
						else
							tmpB[i][j] = tmpB[i][j] - tmpA[m][j]*tmpA[i][m]/tmpA[m][m];
					}
				}
				if(  tmpA[i][m] == 0 )
					continue;
				else
				{
					tmpVb[i] = tmpVb[i] - tmpVa[m]*tmpA[i][m] /tmpA[m][m];// i's row
					for( k=0; k<num; k++)
					{
						tmpD[i][k] = tmpD[i][k] - tmpE[m][k]*tmpA[i][m] /tmpA[m][m];
					}
				}
			}
		}
		tmp = tmpB[m][m];
		/// units
		for( j=m; j<num; j++)
		{
			if( j==m )
				tmpB[m][j] = 1;
			else
				tmpB[m][j] = tmpB[m][j] / tmp;
		}
		tmpVb[m] = tmpVb[m] / tmp;
		for( k=0; k<num; k++)
		{
			tmpD[m][k] = tmpD[m][k] / tmp; //k's row
		}
		tmpA = tmpB;
		tmpVa = tmpVb;
	}
	/// compute xi( xi == vb) and inv(M)
	for( i=0; i<num-1; i++)
	{
		tmpVb[i] = tmpVb[i] - tmpVb[num-1]*tmpB[i][num-1]/tmpB[num-1][num-1];//i's row
		for( k=0; k<num; k++)
		{
			tmpD[i][k] = tmpD[i][k] - tmpD[num-1][k]*tmpB[i][num-1]/tmpB[num-1][num-1]; 
		}
		tmpB[i][num-1] = 0;
	}
	tmpVb[num-1] = tmpVb[num-1]/tmpB[num-1][num-1];
	for( k=0; k<num; k++)
	{
		tmpD[num-1][k] = tmpD[num-1][k]/tmpB[num-1][num-1];
	}
	tmpB[num-1][num-1] = 1;
	vX = tmpVb;
	invM = tmpD;
	return true;
}
/// Compute symbolic matrix inverse with Gaussian-Jordan Elimintaion
template <typename Type>
bool Numeric<Type>::InvMwithGaussJordan( const Matrix<Type> &A)
{
	int num = A.rows();
	Matrix<Type> tmpA = A;
	Matrix<Type> tmpB;
	Matrix<Type> tmpC;
	Matrix<Type> tmpD = eye(num,Type(1));
	Matrix<Type> tmpE = eye(num, Type(1));
	int *maxrow = new int[num-1];
	Data_Pos<Type> maxandpos;
	int i,j,k;
	Type tmp;
	for ( int m = 0; m<(num-1); m++) // the m'th cycle
	{
		/// find the max pivoting element of matrix's m column
		tmpC = CopyFromMatrix(tmpA, m, m, num-1, m);
		maxandpos = FindMaxandPos(abs(tmpC));
		maxrow[m] = maxandpos.row + m;
		/// change the rows 
		tmpA = ExchangeRows(tmpA,m,maxrow[m]);
		tmpB = tmpA;
		/// change the tmpD( for invA ) squence
		tmpD = ExchangeRows(tmpD,m,maxrow[m]);
		tmpE = tmpD;
		for ( i=0; i<num; i++)// row
		{
			if( i != m )
			{
				for( j=m; j<num; j++)// column
				{
					if( j == m )
					{
						if( tmpB[i][j] == 0)
							break;
						else
							tmpB[i][j] = 0;
					}
					else
					{
						if ( tmpA[m][m] == 0)
							return false;
						else
							tmpB[i][j] = tmpB[i][j] - tmpA[m][j]*tmpA[i][m]/tmpA[m][m];
					}
				}
				if(  tmpA[i][m] == 0 )
					continue;
				else
				{
					for( k=0; k<num; k++)
					{
						tmpD[i][k] = tmpD[i][k] - tmpE[m][k]*tmpA[i][m] /tmpA[m][m];
					}
				}
			}
		}
		tmp = tmpB[m][m];
		/// units
		for( j=m; j<num; j++)
		{
			if( j==m )
				tmpB[m][j] = 1;
			else
				tmpB[m][j] = tmpB[m][j] / tmp;
		}
		for( k=0; k<num; k++)
		{
			tmpD[m][k] = tmpD[m][k] / tmp; //k's row
		}
		tmpA = tmpB;
	}
	/// compute xi( xi == vb) and inv(M)
	for( i=0; i<num-1; i++)
	{
		for( k=0; k<num; k++)
		{
			tmpD[i][k] = tmpD[i][k] - tmpD[num-1][k]*tmpB[i][num-1]/tmpB[num-1][num-1]; 
		}
		tmpB[i][num-1] = 0;
	}
	for( k=0; k<num; k++)
	{
		tmpD[num-1][k] = tmpD[num-1][k]/tmpB[num-1][num-1];
	}
	tmpB[num-1][num-1] = 1;
	invM = tmpD;
	return true;
}

/// L-U decomposition of Matrix 
template <typename Type>
bool Numeric<Type>::MatLUdec( const Matrix<Type> &A)
{
	int num = A.rows();
	Matrix<Type> tmpA = A;
	Matrix<Type> tmpB ;
	Matrix<Type> tmpL1 = eye(num, Type(1));
	Matrix<Type> tmpL2 = eye(num, Type(1));
	int *maxrow = new int[num-2];
	int *maxcol = new int[num-2];
	Data_Pos<Type> maxandpos;
	LM = tmpL1;
	int i,j;
	for( j=0; j<num; j++) // colum's
	{
		/// find the max pivoting element of matrix's m column
		tmpB = CopyFromMatrix(tmpA, j, j, num-1, j);
		maxandpos = FindMaxandPos(abs(tmpB));
		maxrow[j] = maxandpos.row + j;
		/// change the rows 
		tmpA = ExchangeRows(tmpA,j,maxrow[j]);
		UM = tmpA;
		for( i=j+1; i<num; i++) // row's
		{
			tmpL1[i][j] = tmpA[i][j]/tmpA[j][j];
			tmpL2[i][j] = -tmpL1[i][j];
		}
		tmpA = tmpL2*tmpA;
		LM = LM*tmpL1;
		UM = tmpL2*UM;
		for( i=j+1; i<num; i++) // row's
		{
			tmpL1[i][j] = 0;
			tmpL2[i][j] = 0;
		}
	} 
	// back the L origin squencese
	for( i=num-2; i>=0; i--)
	{
		if( i != maxrow[i])
		{
			tmpA = CopyFromMatrix(LM, i, i, i, num-1);
			tmpB = CopyFromMatrix(LM, maxrow[i], i, maxrow[i], num-1);
			LM.ReplaceByMatrix( tmpA, maxrow[i], i);
			LM.ReplaceByMatrix( tmpB, i, i);
		}
	}
	return true;
}
/// L-U-P decomposition of Matrix
template <typename Type>
bool Numeric<Type>::MatLUPdec( const Matrix<Type> &A)
{
	int num = A.rows();
	Matrix<Type> tmpA = A;
	Matrix<Type> tmpB ;
	Matrix<Type> tmpL1 = eye(num, Type(1));
	Matrix<Type> tmpL2 = eye(num, Type(1));
	int *maxrow = new int[num-2];
	int *maxcol = new int[num-2];
	Data_Pos<Type> maxandpos;
	LM = tmpL1;
	PM = tmpL1;
	int i,j;
	for( j=0; j<num; j++) // colum's
	{
		/// find the max pivoting element of matrix's m column
		tmpB = CopyFromMatrix(tmpA, j, j, num-1, j);
		maxandpos = FindMaxandPos(abs(tmpB));
		maxrow[j] = maxandpos.row + j;
		/// change the rows 
		tmpA = ExchangeRows(tmpA,j,maxrow[j]);
		UM = tmpA;
		for( i=j+1; i<num; i++) // row's
		{
			tmpL1[i][j] = tmpA[i][j]/tmpA[j][j];
			tmpL2[i][j] = -tmpL1[i][j];
		}
		tmpA = tmpL2*tmpA;
		LM = LM*tmpL1;
		UM = tmpL2*UM;
		for( i=j+1; i<num; i++) // row's
		{
			tmpL1[i][j] = 0;
			tmpL2[i][j] = 0;
		}
	} 
	// back the L origin squencese
	for( i=num-2; i>=0; i--)
	{
		if( i != maxrow[i])
		{
			PM = ExchangeRows(PM,i,maxrow[i]);
			cout << "PM:" << PM;
		}
	}
	return true;
}

///  get the the result of equations of m by m matrix
template <typename Type>
vector<Type> Numeric<Type>::getvX() const
{
	return vX;
}
/// get the the result of a matrix's inv
template <typename Type>
Matrix<Type>  Numeric<Type>::getinvM() const
{
	return invM;
}
/// get the the result of a matrix's inv
template <typename Type>
Matrix<Type>  Numeric<Type>::getMatL() const
{
	return LM;
}
/// get the the result of a matrix's inv
template <typename Type>
Matrix<Type>  Numeric<Type>::getMatU() const
{
	return UM;
}
/// get the the result of a matrix's inv
template <typename Type>
Matrix<Type>  Numeric<Type>::getMatP() const
{
	return PM;
}
/// \brief ����"=",ϵ����ֵ
/// \param x ��ֵ
/// \return ��ֵ�����
template <typename Type>
Type Numeric<Type>::operator=( const Type &x )
{
	return tmp;
}