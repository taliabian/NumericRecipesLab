/*!
* \file NumericTest.cpp
* \brief NumericRecipesLab class testing 
* \author Talia
* \version 1.1
* \date 2015-04-10 
*/

#include <iostream>
#include <iomanip>
#include ".\include\matrix.h"
#include ".\include\numeric.h"

using namespace std;
using namespace matrixlab;
using namespace NumericRecipesLab;

typedef double  Type;
const   int     M = 3;
const   int     N = 3;

int main()
{
	/// 1
	cout << "============ General Gauss Elimnation to solve the Nonsingular matrix equation=========="<<endl;
	cout << "1st: Enter the size of matrix A :" ;
	int tmp = 0;
	cin >> tmp;
	Matrix<Type> gaussA(tmp,tmp);
	cout << "2st: Enter the element of Nonsingular matrix A : \n";
	for (int i=0; i<tmp; i++)
	{
		for (int j=0; j<tmp; j++)
		{
			cin>>gaussA[i][j];
		}
	}
	vector<Type> gaussb(tmp);
	cout << "3st: Enter the element of vector b : \n";
	for (int i=0; i<tmp; i++)
	{
		cin>>gaussb[i];
	}
	Numeric<Type> m1( gaussA ); 
	vector<Type> r1(tmp);
	Matrix<Type> inv;
	Matrix<Type> L, U, P;
	bool yn;
	yn = m1.GenGaussElimation( gaussA, gaussb);
	if( yn == true)
	{
		r1 = m1.getvX();
		cout<< "4st: The result of General Gauss Elimnation is :"<< endl;
		cout<<"vector x:"<<endl;
		cout<<r1;
	}else{
		cout<< "sorry, cannot solve the equation"<<endl;
	}
	/// 2
	cout << "============ Gaussian elimination with full pivoting =========="<<endl;
	yn = m1.FullPivotGaussElimation( gaussA, gaussb);
	if( yn == true)
	{
		r1 = m1.getvX();
		cout<< "The result of Gauss Elimnation with Full Pivoting is :"<< endl;
		cout<<"vector x:"<<endl;
		cout<<r1;
	}else{
		cout<< "sorry, cannot solve the equation"<<endl;
	}
	/// 3
	cout << "============ Gaussian elimination with partital column pivoting =========="<<endl;
	yn = m1.PartialPivotGaussElimation( gaussA, gaussb);
	if( yn == true)
	{
		r1 = m1.getvX();
		cout<< "The result of Gauss Elimnation with Partial Pivoting is :"<< endl;
		cout<<"vector x:"<<endl;
		cout<<r1;
	}else{
		cout<< "sorry, cannot solve the equation"<<endl;
	}
	/// 4
	cout << "============ Gaussian-Jordan elimination =========="<<endl;
	yn = m1.GaussJordanElimation( gaussA, gaussb);
	if( yn == true)
	{
		r1 = m1.getvX();
		cout<< "The result of Gaussian-Jordan elimination  is :"<< endl;
		cout<<"vector x:"<<endl;
		cout<<r1;
		inv = m1.getinvM();
		cout<< "The Inv of Equation Matrix paraments is :"<< endl;
		cout<<"Inv(A):"<<endl;
		cout<<inv;
	}else{
		cout<< "sorry, cannot solve the equation"<<endl;
	}
	/// 5
	cout << "============ Compute symbolic matrix inverse with Gaussian-Jordan Elimintaion =========="<<endl;
	yn = m1.InvMwithGaussJordan( gaussA);
	if( yn == true)
	{
		inv = m1.getinvM();
		cout<< "The Inv of Equation Matrix paraments is :"<< endl;
		cout<<"Inv(A):"<<endl;
		cout<<inv;
	}else{
		cout<< "sorry, cannot solve the equation"<<endl;
	}
	/// 6
	cout << "============= L-U decomposition of Matrix  =========="<<endl;
	yn = m1.MatLUdec( gaussA);
	if( yn == true)
	{
		U = m1.getMatU();
		L = m1.getMatL();
		cout<< "The L of Matrix :"<< endl;
		cout<< L;
		cout<< "The U of Matrix :"<< endl;
		cout<< U;
	}else{
		cout<< "sorry, cannot decomposition the matrix"<<endl;
	}
	cout<< "L*U: " << L*U;
	/// 7
	cout << "=============  L-U-P decomposition of Matrix  =========="<<endl;
	yn = m1.MatLUPdec( gaussA);
	if( yn == true)
	{
		U = m1.getMatU();
		L = m1.getMatL();
		P = m1.getMatP();
		cout<< "The L of Matrix :"<< endl;
		cout<< L;
		cout<< "The U of Matrix :"<< endl;
		cout<< U;
		cout<< "The P of Matrix :"<< endl;
		cout<< P;
	}else{
		cout<< "sorry, cannot decomposition the matrix"<<endl;
	}
	cout<< "L*U: " << L*U;
	cout << "P*A" << P*gaussA;
	
	/// 8
	cout << "=============  solve the equations with L-U-P decomposition of Matrix  =========="<<endl;
	yn = m1.LUPsolveFun( gaussA, gaussb);
	if( yn == true)
	{
		r1 = m1.getvX();
		cout<< "The result is :"<< endl;
		cout<<"vector x:"<<endl;
		cout<<r1;
	}else{
		cout<< "sorry, cannot decomposition the matrix"<<endl;
	}

	cout << "============ Solve the equations with Least Squares use L-U-P Decompositions and Pseudo inverse =========="<<endl;
	cout << "1st: Enter the rows and cols of matrix A :" ;
	int rows, cols;
	cout << "rows: "; 
	cin >> rows; 
	cout << "cols: "; 
	cin >> cols;
	Matrix<Type> gaussA1(rows, cols);
	cout << "2st: Enter the element of Nonsingular matrix A : \n";
	for (int i=0; i<rows; i++)
	{
		for (int j=0; j<cols; j++)
		{
			cin>>gaussA1[i][j];
		}
	}
	vector<Type> gaussb1(rows);
	cout << "3st: Enter the element of vector b : \n";
	for (int i=0; i<tmp; i++)
	{
		cin>>gaussb1[i];
	}
	Numeric<Type> m2( gaussA1 ); 
	Matrix<Type> m3( rows, rows);
	m3  = m2.LeasetSquaresSolveFun( gaussA1, gaussb1 );

	return 0;
}  