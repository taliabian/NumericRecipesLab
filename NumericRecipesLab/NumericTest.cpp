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
	Numeric<Type> m1;
	bool yn;
	yn = m1.GenGaussElimation( gaussA, gaussb);
	if( yn == true)
	{
		vector<Type> r1 = m1.getvX();
		cout<< "4st: The result of General Gauss Elimnation is :"<< endl;
		cout<<"vector x:"<<endl;
		cout<<r1;
	}else{
		cout<< "sorry, cannot solve the equation"<<endl;
	}


	return 0;
}