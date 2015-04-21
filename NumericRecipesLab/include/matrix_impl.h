/*!
* \file matrix_impl.h
* \brief Implementation for MatrixLab class
* \author Talia
* \version 1.1
* \date 2015-04-09 
*/

/// \brief �����ʼ��(pv0, prow0, prow1)
/// \param rows ��������
/// \param columns ��������
template <typename Type>
void Matrix<Type>::init( int rows, int columns )
{
	nRow = rows;
	nColumn = columns;
	nTotal = nRow * nColumn;

	pv0 = new Type[nTotal];
	prow0 = new Type*[nRow];
	prow1 = new Type*[nRow];
	Type *p = pv0;
	for( int i=0; i<nRow; ++i )
	{
		prow0[i] = p;
		prow1[i] = p-1;
		p += nColumn;
	}

	prow1--;
}

/// \brief �����������ݳ�ʼ����������pv0
/// \param v ����ָ��
template <typename Type>
inline void Matrix<Type>::copyFromArray( const Type *v )
{
	for( long i=0; i<nTotal; ++i )
		pv0[i] = v[i];
}
/// \brief ���þ�������Ϊͬһ��x
/// \param x ����ϵ��
template <typename Type>
inline void Matrix<Type>::setByScalar( const Type &x )
{
	for( long i=0; i<nTotal; ++i )
		pv0[i] = x;
}
/// \brief ��������(pv0, prow0, prow1)
template <typename Type>
void Matrix<Type>::destroy()
{
	if( pv0 == NULL )
		return ;
	else
		delete []pv0;

	if( prow0 != NULL )
		delete []prow0;

	prow1++;
	if( prow1 != NULL )
		delete []prow1;
}

/// \brief Ĭ�Ϲ��캯��
template <typename Type>
Matrix<Type>::Matrix()
: pv0(0), prow0(0), prow1(0), nRow(0), nColumn(0), nTotal(0)
{
}
/// \brief ���캯��,������֪�����ʼ��
/// \param A ��֪����
/// \see init(),copyFromArray()
template <typename Type>
Matrix<Type>::Matrix( const Matrix<Type> &A )
{
	init( A.nRow, A.nColumn );
	copyFromArray( A.pv0 );
}

/// \brief ���캯��,��֪����������������ϵ��,��ʼ������
/// \param rows ������
/// \param columns ������
/// \param x ����ϵ��,Ĭ��Ϊ0
/// \see init(),setByScalar()
template <typename Type>
Matrix<Type>::Matrix( int rows, int columns, const Type &x = Type(0) )
{
	init( rows,columns );
	setByScalar(x);
}
/// \brief ���캯��,��֪��������������������,��ʼ������
/// \param rows ������
/// \param columns ������
/// \param array ��������ָ��
/// \see init(),copyFromArray()
template <typename Type>
Matrix<Type>::Matrix( int rows, int columns, const Type *arrays )
{
	init( rows,columns );
	copyFromArray( arrays );
}
/// \brief ����������
template <typename Type>
Matrix<Type>::~Matrix()
{
	destroy();
}
/// \brief ����"=",����ֵ,such as: matrixB = matrixA
/// \param A ��ֵ�Ҳ���
/// \return ��ֵ�����
/// \remarks ���maxtrixB �� maxtrix A��ά������,��ԭ��matrixB����,��ֵΪmatrixA
template <typename Type>
Matrix<Type>& Matrix<Type>::operator=( const Matrix<Type> &A )
{
	if( pv0 == A.pv0 )
		return *this;

	if( nRow == A.nRow && nColumn == A.nColumn )
		copyFromArray( A.pv0 );
	else
	{
		destroy();
		init( A.nRow, A.nColumn );
		copyFromArray( A.pv0 );
	}

	return *this;
}
/// \brief ����"=",ϵ����ֵ,such as: matrixA = matrix(x)
/// \param x ����ϵ��
/// \return ��ֵ�����
template <typename Type>
inline Matrix<Type>& Matrix<Type>::operator=( const Type &x )
{
	setByScalar( x );
	return *this;
}

/// \brief ����"[]",����'0'�������ݷ���. 
/// \param x ����
/// \return �����i�е��׵�ַ
template <typename Type>
inline Type* Matrix<Type>::operator[]( int i )
{
	return prow0[i];
}
template <typename Type>
inline const Type* Matrix<Type>::operator[]( int i ) const
{
	return prow0[i];
}
/// \brief ����"()",����'1'�����ݷ���. 
/// \param row ����
/// \param column ����
/// \return �����i�е�j������
template <typename Type>
inline Type& Matrix<Type>::operator()( int row, int column )
{
	return  prow1[row][column];
}

template <typename Type>
inline const Type& Matrix<Type>::operator()( int row, int column ) const
{
	return  prow1[row][column];
}
/// \brief ����"*",���������׵�ַ����. 
/// \return ���������׵�ַpv0
template <typename Type>
inline Matrix<Type>::operator Type*()
{
	return pv0;
}

template <typename Type>
inline Matrix<Type>::operator const Type*() const
{
	return pv0;
}
/// \brief ���ؾ����С
template <typename Type>
inline long Matrix<Type>::size() const
{
	return nTotal;
}
/// \brief ���ؾ���ά��
/// \param dimension '1'��������,'2'��������
template <typename Type>
int Matrix<Type>::dim( int dimension ) const
{

	if( dimension == 1 )
		return nRow;
	else if( dimension == 2 )
		return nColumn;
	else
		return 0;
}
/// \brief ���ؾ�������
template <typename Type>
inline int Matrix<Type>::rows() const
{
    return nRow;
}
/// \brief ���ؾ�������
template <typename Type>
inline int Matrix<Type>::cols() const
{
    return nColumn;
}
/// \brief �����趨�����С
/// \param rows ����
/// \param columns ����
template <typename Type>
Matrix<Type>& Matrix<Type>::resize( int rows, int columns )
{
	if(  rows == nRow && columns == nColumn )
		return *this;

	destroy();
	init( rows, columns );

	return *this;
}
/// \brief ��þ���������
/// \param row ����(0��)
/// \return ����������
template <typename Type>
vector<Type> Matrix<Type>::getRow( int row ) const
{

	vector<Type> tmp( nColumn );
	for( int j=0; j<nColumn; ++j )
		tmp[j] = prow0[row][j];

	return tmp;
}
/// \brief ��þ���������
/// \param column ����(0��)
/// \return ����������
template <typename Type>
vector<Type> Matrix<Type>::getColumn( int column ) const
{
	vector<Type> tmp( nRow );
	for( int i=0; i<nRow; ++i )
		tmp[i] = prow0[i][column];

	return tmp;
}
/// \brief �趨����������
/// \param row ����(0��)
/// \param v ������
template <typename Type>
void Matrix<Type>::setRow( const vector<Type> &v, int row )
{
	for( int j=0; j<nColumn; ++j )
		prow0[row][j] = v[j];
}
/// \brief �趨����������
/// \param column ����(0��)
/// \param v ������
template <typename Type>
void Matrix<Type>::setColumn( const vector<Type> &v, int column )
{
	for( int i=0; i<nRow; ++i )
		prow0[i][column] = v[i];
}
/// \brief ����'+=', �����Ԫ�ؼ�ϵ��
/// \param x ϵ��
template <typename Type>
Matrix<Type>& Matrix<Type>::operator+=( const Type &x )
{
    Type **rowPtr = prow0;
    Type *colPtr = 0;

    for( int i=0; i<nRow; ++i )
    {
        colPtr = *rowPtr++;
        for( int j=0; j<nColumn; ++j )
            *colPtr++ += x;
    }
	return *this;
}
/// \brief ����'+=', ����A��Ԫ�ؼӶ�Ӧ����B��Ԫ��ϵ��
/// \param rhs ����B
template <typename Type>
Matrix<Type>& Matrix<Type>::operator+=( const Matrix<Type> &rhs )
{

    Type **rowPtrL = prow0;
    Type *colPtrL = 0;
    Type **rowPtrR = rhs.prow0;
    const Type *colPtrR = 0;

    for( int i=0; i<nRow; ++i )
    {
        colPtrL = *rowPtrL++;
        colPtrR = *rowPtrR++;
        for( int j=0; j<nColumn; ++j )
            *colPtrL++ += *colPtrR++;
    }

	return *this;
}
/// \brief ����'-=', �����Ԫ�ؼ�ϵ��
/// \param x ϵ��
template <typename Type>
Matrix<Type>& Matrix<Type>::operator-=( const Type &x )
{
    Type **rowPtr = prow0;
    Type *colPtr = 0;

    for( int i=0; i<nRow; ++i )
    {
        colPtr = *rowPtr++;
        for( int j=0; j<nColumn; ++j )
            *colPtr++ -= x;
    }

	return *this;
}
/// \brief ����'-=', ����A��Ԫ�ؼ���Ӧ����B��Ԫ��ϵ��
/// \param rhs ����B
template <typename Type>
Matrix<Type>& Matrix<Type>::operator-=( const Matrix<Type> &rhs )
{
    Type **rowPtrL = prow0;
    Type *colPtrL = 0;
    Type **rowPtrR = rhs.prow0;
    const Type *colPtrR = 0;

    for( int i=0; i<nRow; ++i )
    {
        colPtrL = *rowPtrL++;
        colPtrR = *rowPtrR++;
        for( int j=0; j<nColumn; ++j )
            *colPtrL++ -= *colPtrR++;
    }

	return *this;
}
/// \brief ����'*=', �����Ԫ�س�ϵ��
/// \param x ϵ��
template <typename Type>
Matrix<Type>& Matrix<Type>::operator*=( const Type &x )
{
    Type **rowPtr = prow0;
    Type *colPtr = 0;

    for( int i=0; i<nRow; ++i )
    {
        colPtr = *rowPtr++;
        for( int j=0; j<nColumn; ++j )
            *colPtr++ *= x;
    }

	return *this;
}
/// \brief ����'*=', ����A��Ԫ�س˶�Ӧ����B��Ԫ��ϵ��
/// \param rhs ����B
template <typename Type>
Matrix<Type>& Matrix<Type>::operator*=( const Matrix<Type> &rhs )
{
    Type **rowPtrL = prow0;
    Type *colPtrL = 0;
    Type **rowPtrR = rhs.prow0;
    const Type *colPtrR = 0;

    for( int i=0; i<nRow; ++i )
    {
        colPtrL = *rowPtrL++;
        colPtrR = *rowPtrR++;
        for( int j=0; j<nColumn; ++j )
            *colPtrL++ *= *colPtrR++;
    }

	return *this;
}
/// \brief ����'/=', �����Ԫ�س�ϵ��
/// \param x ϵ��
template <typename Type>
Matrix<Type>& Matrix<Type>::operator/=( const Type &x )
{
    Type **rowPtr = prow0;
    Type *colPtr = 0;

    for( int i=0; i<nRow; ++i )
    {
        colPtr = *rowPtr++;
        for( int j=0; j<nColumn; ++j )
            *colPtr++ /= x;
    }

	return *this;
}
/// \brief ����'/=', ����A��Ԫ�س���Ӧ����B��Ԫ��ϵ��
/// \param rhs ����B
template <typename Type>
Matrix<Type>& Matrix<Type>::operator/=( const Matrix<Type> &rhs )
{

    Type **rowPtrL = prow0;
    Type *colPtrL = 0;
    Type **rowPtrR = rhs.prow0;
    const Type *colPtrR = 0;

    for( int i=0; i<nRow; ++i )
    {
        colPtrL = *rowPtrL++;
        colPtrR = *rowPtrR++;
        for( int j=0; j<nColumn; ++j )
            *colPtrL++ /= *colPtrR++;
    }

	return *this;
}
/// \brief ����'<<', �������A��Ԫ��			 
/// \param A ����A
template <typename Type>
ostream& operator<<( ostream &out, const Matrix<Type> &A )
{
	int rows = A.rows();
	int columns = A.cols();

	out << "size: " << rows << " by " << columns << "\n";
	for( int i=0; i<rows; ++i )
	{
		for( int j=0; j<columns; ++j )
			out << A[i][j] << "\t";
		out << "\n";
	}

	return out;
}
/// \brief ����'<<', �������v��Ԫ��		 
/// \param v ����v
template<typename Type>
ostream& operator<<( ostream &out, const vector<Type> &v )
{
	int rows = v.size();

	out << "size: " << rows << " by " << 1 << "\n";
	for( int i=0; i<rows; ++i )
	{
			out << v.at(i) << "\n";
	}
	return out;
}
/// \brief ����'<<', ���Data_Pos��Ԫ��		 
/// \param v ����v
template<typename Type>
ostream& operator<<( ostream &out, const Data_Pos<Type> &d )
{
	out << "Data: " << d.Data << ", Pos: (" << d.row << ", " << d.col <<")";
	return out;
}
/// \brief ����'>>', �������A��Ԫ��	 
/// \param A �������
template <typename Type>
istream& operator>>( istream &in, Matrix<Type> &A )
{
	int rows, columns;
	in >> rows >> columns;

	if( !( rows == A.rows() && columns == A.cols() ) )
		A.resize( rows, columns );

	for( int i=0; i<rows; ++i )
	{
		for( int j=0; j<columns; ++j )
			in >> A[i][j];
		cout<<endl;
	}
	return in;
}
/// \brief ����'-', ����A��Ԫ��ȡ��	
/// \param A ����A
template<typename Type>
Matrix<Type> operator-( const Matrix<Type> &A )
{
	int rows = A.rows();
	int columns = A.cols();

	Matrix<Type> tmp( rows, columns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			tmp[i][j] = -A[i][j];

	return tmp;
}
/// \brief ����'+', ����A��Ԫ�ؼ�ϵ��
/// \param A ����A
/// \param x ϵ��
template<typename Type>
inline Matrix<Type> operator+( const Matrix<Type> &A, const Type &x )
{
	Matrix<Type> tmp( A );
	return tmp += x;
}
/// \brief ����'+', ϵ���Ӿ���A��Ԫ��
/// \param A ����A
/// \param x ϵ��
template<typename Type>
inline Matrix<Type> operator+( const Type &x, const Matrix<Type> &A )
{
	return A + x;
}
/// \brief ����'+', ����A1��Ԫ�ؼӶ�Ӧ����A2��Ԫ��ϵ��
/// \param A1 ����A1
/// \param A2 ����A2
template<typename Type>
inline Matrix<Type> operator+( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	Matrix<Type> tmp( A1 );
	return tmp += A2;
}
/// \brief ����'-', ����A��Ԫ�ؼ�ϵ��
/// \param A ����A
/// \param x ϵ��
template<typename Type>
inline Matrix<Type> operator-( const Matrix<Type> &A, const Type &x )
{
	Matrix<Type> tmp( A );
	return tmp -= x;
}
/// \brief ����'-', ϵ��������A��Ԫ��
/// \param A ����A
/// \param x ϵ��
template<typename Type>
inline Matrix<Type> operator-( const Type &x, const Matrix<Type> &A )
{
	Matrix<Type> tmp( A );
	return -tmp += x;
}
/// \brief ����'-', ����A1��Ԫ�ؼ���Ӧ����A2��Ԫ��ϵ��
/// \param A1 ����A1
/// \param A2 ����A2
template<typename Type>
inline Matrix<Type> operator-( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	Matrix<Type> tmp( A1 );
	return tmp -= A2;
}
/// \brief ����'*', �����Ԫ�س�ϵ��
/// \param A ����
/// \param x ϵ��
template <typename Type>
inline Matrix<Type> operator*( const Matrix<Type> &A, const Type &x )
{
	Matrix<Type> tmp( A );
	return tmp *= x;
}
/// \brief ����'*', ϵ���˾����Ԫ��
/// \param A ����
/// \param x ϵ��
template <typename Type>
inline Matrix<Type> operator*( const Type &x, const Matrix<Type> &A )
{
	return A * x;
}
/// \brief ����'*', ����A1��Ԫ�س˾���A2
/// \param A1 ����A1
/// \param A1 ����A2
template <typename Type>
Matrix<Type> operator*( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	int rows = A1.rows();
	int columns = A2.cols();

	Matrix<Type> tmp( rows, columns );
    optMult( A1, A2, tmp );

	return tmp;
}
/// \brief ����'*', ����A��Ԫ�س�����v
/// \param A ���� 
/// \param v ���� 
template <typename Type>
vector<Type> operator*( const Matrix<Type> &A, const vector<Type> &b )
{

	int rows = A.rows();

	vector<Type> tmp(rows);

    optMult( A, b, tmp );

	return tmp;
}
/// \brief ����'*', ����v��Ԫ�س�ϵ��
/// \param x ϵ�� 
/// \param v ���� 
template<typename Type>
vector<Type> operator*( const Type &x, const vector<Type> &v )
{
	vector<Type> tmpv(v.size());
	for( int i=0; i<v.size(); i++)	
	{
		tmpv[i] = x*v[i];
	}
	return tmpv;
}

template<typename Type>
vector<Type> operator-( const vector<Type> &v1, const vector<Type> &v2 )
{
	vector<Type> tmpv(v1.size());
	for( int i=0; i<v1.size(); i++)	
	{
		tmpv.at(i) = v1.at(i) - v2.at(i);
	}
	return tmpv;
}

/// \brief ����'/', �����Ԫ�س�ϵ��
/// \param A ����
/// \param x ϵ��
template <typename Type>
inline Matrix<Type> operator/( const Matrix<Type> &A, const Type &x )
{
	Matrix<Type> tmp( A );
	return tmp /= x;
}
/// \brief ����'/', ������Ԫ�س�ϵ��
/// \param v ����
/// \param x ϵ��
template<typename Type>
vector<Type> operator/( const vector<Type> &v, const Type &x )
{
	int msize = v.size();
	vector<Type> tmp(msize);
	for (int i = 0; i < msize; i++)
	{
		tmp.at(i) = v.at(i)/x;
	}
	return tmp;
}
/// \brief ����'/', ϵ���������Ԫ��
/// \param A ����
/// \param x ϵ��
template <typename Type>
Matrix<Type> operator/( const Type &x, const Matrix<Type> &A )
{
	int rows = A.rows();
	int clumns = A.cols();

	Matrix<Type> tmp( rows,clumns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<clumns; ++j )
			tmp[i][j] = x / A[i][j];

	return tmp;
}
/// \brief ����'*', ����A��Ԫ�س˾���B
/// \param A ����
/// \param B ����
/// \param C �������
template <typename Type>
Matrix<Type>& optMult( const Matrix<Type> &A, const Matrix<Type> &B,
                    Matrix<Type> &C )
{
    int M = A.rows();
    int N = B.cols();
    int K = A.cols();

    C.resize( M, N );
    Type        sum;
    const Type  *pRow,
                *pCol;

    for( int i=0; i<M; i++ )
        for( int j=0; j<N; ++j )
        {
            pRow  = &A[i][0];
            pCol  = &B[0][j];
            sum = 0;

            for( int k=0; k<K; ++k )
            {
                sum += (*pRow) * (*pCol);
                pRow++;
                pCol += N;
            }
            C[i][j] = sum;
        }
    return C;
}
/// \brief ����'*', ����A��Ԫ�س�����b
/// \param A ���� 
/// \param b ���� 
/// \param c 
template <typename Type>
vector<Type>& optMult( const Matrix<Type> &A, const vector<Type> &b,
                    vector<Type> &c )
{
    int M = A.rows();
    int N = A.cols();

    c.resize( M );
    Type        sum;
    const Type  *pRow,
                *pCol;

    for( int i=0; i<M; i++ )
    {
        pRow  = &A[i][0];
        pCol  = &b[0];
        sum = 0;

        for( int j=0; j<N; ++j )
        {
            sum += (*pRow) * (*pCol);
            pRow++;
            pCol++;
        }
        c[i] = sum;
    }
    return c;
}
/// \brief ����A1�˾���A2
/// \param A1 ����
/// \param A2 ����
/// \return  �������
template<typename Type>
inline Matrix<Type> elemMult( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	Matrix<Type> tmp( A1 );
	return tmp *= A2;
}
/// \brief ����A1�˾���A2
/// \param A1 ����
/// \param A2 ����
/// \return  A1�������
template <typename Type>
inline Matrix<Type>& elemMultEq( Matrix<Type> &A1, const Matrix<Type> &A2 )
{
    return A1 *= A2;
}
/// \brief ����A1��Ԫ�س���Ӧ����A2��Ԫ��ϵ��
/// \param A1 ����
/// \param A2 ����
/// \return A1/A2
template <typename Type>
inline Matrix<Type> elemDivd( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	Matrix<Type> tmp( A1 );
	return tmp /= A2;
}
/// \brief ����A1��Ԫ�س���Ӧ����A2��Ԫ��ϵ��
/// \param A1 ����
/// \param A2 ����
/// \return A1
template <typename Type>
inline Matrix<Type>& elemDivdEq( Matrix<Type> &A1, const Matrix<Type> &A2 )
{
    return A1 /= A2;
}
/// \brief ʵ������Aת��
/// \param A ʵ����
/// \return trT(A)
template <typename Type>
Matrix<Type> trT( const Matrix<Type> &A )
{
	int rows = A.cols();
	int clumns = A.rows();

	Matrix<Type> tmp( rows, clumns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<clumns; ++j )
			tmp[i][j] = A[j][i];

	return tmp;
}
/// \brief ��������Aת��
/// \param A ��������
/// \return trH(A)
template <typename Type>
Matrix<Type> trH( const Matrix<Type> &A )
{
	int rows = A.cols();
	int clumns = A.rows();

	Matrix<Type> tmp( rows, clumns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<clumns; ++j )
			tmp[i][j] = conj(A[j][i]);

	return tmp;
}
/// \brief A1��ת�õ��A2:(A1^T).*A2
/// \param A1 ����
/// \param A2 ����
/// \return (A1^T).*A2
template <typename Type>
Matrix<Type> trMult( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{
	int rows = A1.cols();
	int columns = A2.cols();
	int K = A1.rows();

	Matrix<Type> tmp( rows, columns );

    for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			for( int k=0; k<K; ++k )
			   tmp[i][j] += A1[k][i] * A2[k][j];

	return tmp;
}
/// \brief A��ת�õ��v:(A^T).*v
/// \param A ����
/// \param v ����
/// \return (A^T).*v
template <typename Type>
vector<Type> trMult( const Matrix<Type> &A, const vector<Type> &v )
{
	int rows = A.rows();
	int columns = A.cols();

	vector<Type> tmp( columns );
    for( int i=0; i<columns; ++i )
		for( int j=0; j<rows; ++j )
			tmp[i] += A[j][i] * v[j];

	return tmp;
}
/// \brief A1���A2��ת��:A1.*(A2^T)
/// \param A1 ����
/// \param A2 ����
/// \return A1.*(A2^T)
template <typename Type>
Matrix<Type> multTr( const Matrix<Type> &A1, const Matrix<Type> &A2 )
{

	int rows = A1.rows();
	int columns = A2.rows();
	int K = A1.cols();

	Matrix<Type> tmp( rows, columns );
    for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			for( int k=0; k<K; ++k )
			   tmp[i][j] += A1[i][k] * A2[j][k];

	return tmp;
}
/// \brief ����a�������b��ת��:a.*(b^T)
/// \param a ����
/// \param b ����
/// \return a.*(b^T)
template <typename Type>
Matrix<Type> multTr( const vector<Type> &a, const vector<Type> &b )
{
	int rows = a.dim();
	int columns = b.dim();

	Matrix<Type> tmp( rows, columns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			tmp[i][j] = a[i]*b[j];

	return tmp;
}
/// \brief ����a��ת�õ������b:(a^T).*b
/// \param a ����
/// \param b ����
/// \return (a^T).*b
template<typename Type>
Type trMult( const vector<Type> &a, const vector<Type> &b )
{
	Type tmp = 0;
	int rows = a.size();
	for( int i = 0; i<rows; i++)
		tmp += a[i]*b[i];
	return tmp;
}
/// \brief ��������A1�Ĺ���ת�õ�˸�������A2:(A1^H).*A2
/// \param A1 ��������
/// \param A2 ��������
/// \return (A1^H).*A2
template <typename Type>
Matrix<complex<Type> > trMult( const Matrix<complex<Type> > &A1,
                               const Matrix<complex<Type> > &A2 )
{

	int rows = A1.cols();
	int columns = A2.cols();
	int K = A1.rows();

	Matrix<complex<Type> > tmp( rows, columns );
    for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			for( int k=0; k<K; ++k )
			   tmp[i][j] += conj(A1[k][i]) * A2[k][j];

	return tmp;
}
/// \brief ��������A�Ĺ���ת�õ�˸�������v:(A^H).*v
/// \param A ��������
/// \param v ��������
/// \return (A^H).*v
template <typename Type>
vector<complex<Type> > trMult( const Matrix<complex<Type> > &A, const vector<complex<Type> > &v )
{
	int rows = A.rows();
	int columns = A.cols();

	vector<complex<Type> > tmp( columns );

    for( int i=0; i<columns; ++i )
		for( int j=0; j<rows; ++j )
			tmp[i] += conj(A[j][i]) * v[j];

	return tmp;
}
/// \brief ��������A1��˸�������A2�Ĺ���ת��:A1.*(A2^H)
/// \param A1 ��������
/// \param A2 ��������
/// \return A1.*(A2^H)
template <typename Type>
Matrix<complex<Type> > multTr( const Matrix<complex<Type> > &A1,
                               const Matrix<complex<Type> > &A2 )
{

	int rows = A1.rows();
	int columns = A2.rows();
	int K = A1.cols();

	Matrix<complex<Type> > tmp( rows, columns );
    for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			for( int k=0; k<K; ++k )
			   tmp[i][j] += A1[i][k] * conj(A2[j][k]);

	return tmp;
}
/// \brief ��������a��˸�������b�Ĺ���ת��:a.*(b^H)
/// \param a ��������
/// \param b ��������
/// \return a.*(b^H)
template <typename Type>
Matrix<complex<Type> > multTr( const vector<complex<Type> > &a,
                               const vector<complex<Type> > &b )
{
	int rows = a.dim();
	int columns = b.dim();

	Matrix<complex<Type> > tmp( rows, columns );
	for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			tmp[i][j] = a[i]*conj(b[j]);

	return tmp;
}
/// \brief ��������a�Ĺ���ת�õ�˸�������b:(a^H).*b
/// \param a ��������
/// \param b ��������
/// \return (a^H).*b
template <typename Type>
complex<Type>  trMult( const vector<complex<Type> > &a, const vector<complex<Type> > &b )
{
	int rows = a.dim();

	complex<Type> tmp;
	for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			tmp += conj(a[i])*b[j];

	return tmp;
}
/// \brief ���ɵ�λ����
/// \param x ά��
template <typename Type>
Matrix<Type> eye( int N, const Type &x )
{
    Matrix<Type> tmp( N, N );
	for( int i=0; i<N; ++i )
		tmp[i][i] = x;

	return tmp;
}
/// \brief ���������
/// \param x ά��
template <typename Type>
Matrix<Type> zeros( int row, int col )
{
    Matrix<Type> tmp( row, col );
	return tmp;
}
/// \brief ���ɾ���A�Խ�Ԫ��
/// \param A ����
template <typename Type>
vector<Type> diag( const Matrix<Type> &A )
{
	int nColumn = A.rows();
	if( nColumn > A.cols() )
		nColumn = A.cols();

	vector<Type> tmp( nColumn );
	for( int i=0; i<nColumn; ++i )
		tmp[i] = A[i][i];

	return tmp;
}
/// \brief �����Խ�Ԫ�����ɾ���A
/// \param d �Խ�Ԫ������
template <typename Type>
Matrix<Type> diag( const vector<Type> &d )
{
	int N = d.size();

	Matrix<Type> tmp( N, N );
	for( int i=0; i<N; ++i )
		tmp[i][i] = d[i];

	return tmp;
}
/// \brief �����Frobenius���� = sqrt(sum(A(i,j)*A(i,j))) 
/// \param A ����
template <typename Type>
Type norm( const Matrix<Type> &A )
{
	int m = A.rows();
	int n = A.cols();

	Type sum = 0;
	for( int i=1; i<=m; ++i )
		for( int j=1; j<=n; ++j )
            sum += A(i,j) * A(i,j);

	return sqrt(sum);
}
/// \brief ������Frobenius������ģ��
/// \param v ����
template <typename Type>
Type norm( const vector<Type> &v )
{
	int m = v.size();

	Type sum = 0;
	for( int i=0; i<m; ++i )
            sum += v.at(i) * v.at(i);

	return sqrt(sum);
}
/// \brief ���������Frobenius���� = sqrt(sum(A(i,j)*conj(A(i,j)))) 
/// \param A ����
template <typename Type>
Type norm( const Matrix<complex<Type> > &A )
{
	int m = A.rows();
	int n = A.cols();
	complex<Type> tmp1,tmp2,tmp3;
	Type sum = 0;
	for( int i=1; i<=m; ++i )
		for( int j=1; j<=n; ++j)
			sum += norm(A(i,j));//c++ complex function

	return sqrt(sum);
}
/// \brief ��������lhs��rhs������
/// \param lhs ����
/// \param rhs ����
template <typename Type> void swap( Matrix<Type> &lhs, Matrix<Type> &rhs )
{
    int m = lhs.rows();
	int n = lhs.cols();

	for( int i=1; i<=m; ++i )
		for( int j=1; j<=n; ++j )
            swap( lhs(i,j), rhs(i,j) );
}
/// \brief �������Ԫ�����
/// \param A ����
/// \return ����
template <typename Type>
vector<Type> sum( const Matrix<Type> &A )
{
	int m = A.rows();
	int n = A.cols();
	vector<Type> sum(n);

	for( int j=1; j<=n; ++j )
		for( int i=1; i<=m; ++i )
            sum.at(j-1) += A(i,j);

	return sum;
}
/// \brief �������Ԫ����Сֵ
/// \param A ����
/// \return ����
template <typename Type>
vector<Type> min( const Matrix<Type> &A )
{
	int m = A.rows();
	int n = A.cols();
	vector<Type> Min(n);

	for( int j=1; j<=n; ++j )
	{
	    Type tmp = A(1,j);
        for( int i=2; i<=m; ++i )
            if( tmp > A(i,j) )
                tmp = A(i,j);
        Min.at(j-1) = tmp;
	}

	return Min;
}

/// \brief �������Ԫ�����ֵ
/// \param A ����
/// \return ����
template <typename Type>
vector<Type> max( const Matrix<Type> &A )
{
	int m = A.rows();
	int n = A.cols();
	vector<Type> Max(n);

	for( int j=1; j<=n; ++j )
	{
	    Type tmp = A(1,j);
        for( int i=2; i<=m; ++i )
            if( tmp < A(i,j) )
                tmp = A(i,j);
        Max.at(j-1) = tmp;
	}

	return Max;
}

/// \brief �������Ԫ��ƽ��ֵ
/// \param A ����
/// \return ����
template <typename Type>
inline vector<Type> mean( const Matrix<Type> &A )
{
	return sum(A) / Type(A.rows());
}

/// \brief ��ʵ�����츴������2==>2+0*i
/// \param rA ����
/// \return ��������
template <typename Type>
Matrix<complex<Type> > complexMatrix( const Matrix<Type> &rA )
{
	int rows = rA.rows();
	int columns = rA.cols();

    Matrix<complex<Type> > cA( rows, columns );
    for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			cA[i][j] = rA[i][j];

    return cA;
}
/// \brief ��ʵ����mR��mI���츴������mR+mI*i
/// \param rR ����
/// \param rI ����
/// \return ��������
template <typename Type>
Matrix<complex<Type> > complexMatrix( const Matrix<Type> &mR,
                                      const Matrix<Type> &mI )
{
	int rows = mR.rows();
	int columns = mR.cols();

    Matrix<complex<Type> > cA( rows, columns );
    for( int i=0; i<rows; ++i )
		for( int j=0; j<columns; ++j )
			cA[i][j] = complex<Type>( mR[i][j], mI[i][j] );

    return cA;
}
/// \brief ����ľ���ֵ����
/// \param A ����
/// \return ����
template <typename Type>
Matrix<Type> abs( const Matrix<complex<Type> > &A )
{
    int m = A.rows(),
        n = A.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = abs( A[i][j] );

    return tmp;
}

/// \brief ���������ģ����
/// \param A ����
/// \return ʵ����
template <typename Type>
Matrix<Type> abs( const Matrix<Type> &A )
{
    int m = A.rows(),
        n = A.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
		{
			if ( A[i][j] < 0)
				tmp[i][j] = -A[i][j];
			else
				tmp[i][j] = A[i][j];
		}

    return tmp;
}

/// \brief ��������ļ�����ǶȾ���
/// \param A ����
/// \return ʵ����
template <typename Type>
Matrix<Type> arg( const Matrix<complex<Type> > &A )
{
    int m = A.rows(),
        n = A.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = arg( A[i][j] );

    return tmp;
}

/// \brief ���������ʵ������
/// \param A ����
/// \return ʵ����
template <typename Type>
Matrix<Type> real( const Matrix<complex<Type> > &A )
{
    int m = A.rows(),
        n = A.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = A[i][j].real();

    return tmp;
}

/// \brief ����������鲿����
/// \param A ����
/// \return ʵ����
template <typename Type>
Matrix<Type> imag( const Matrix<complex<Type> > &A )
{
    int m = A.rows(),
        n = A.cols();
    Matrix<Type> tmp( m, n );

    for( int i=0; i<m; ++i )
        for( int j=0; j<n; ++j )
            tmp[i][j] = A[i][j].imag();

    return tmp;
}

template <typename Type>
Matrix<Type> strcatMatrix( const Matrix<Type> &A1, const Matrix<Type> &A2)
{
	int rows = A1.rows();
	int colums1 = A1.cols();
	int colums = A1.cols() + A2.cols();
	Matrix<Type> tmp( rows, colums );
	for ( int i=0; i<rows; i++)
		for (int j=0; j<colums; j++)
			if ( j<colums1 )
				 tmp[i][j] = A1[i][j];
			else
				tmp[i][j] = A2[i][j-colums1];
	return tmp;
}

template <typename Type>
Matrix<Type> strcatMatrix( const Matrix<Type> &A, const vector<Type> &v)
{
	int rows = A.rows();
	int colums = A.cols() + 1;
	Matrix<Type> tmp( rows, colums );
	for ( int i=0; i<rows; i++)
		for (int j=0; j<colums; j++)
			if ( j<colums-1 )
				tmp[i][j] = A[i][j];
			else
				tmp[i][j] = v[i];
	return tmp;
}

template <typename Type>
Matrix<Type> strcatMatrix( const vector<Type> &v, const Matrix<Type> &A)
{
	int rows = A.rows();
	int colums = A.cols() + 1;
	Matrix<Type> tmp( rows, colums );
	for ( int i=0; i<rows; i++)
		for (int j=0; j<colums; j++)
			if ( j==0 )
				tmp[i][j] = v[i];
			else
				tmp[i][j] = A[i][j-1];
	return tmp;
}

template <typename Type>
Matrix<Type> strcatMatrix( const vector<Type> &v1, const vector<Type> &v2)
{
	int rows = v1.size();
	Matrix<Type> tmp( rows, 2 );
	for ( int i=0; i<rows; i++)
		for ( int j=0; j<2; j++)
			if ( j==0 )
				tmp[i][j] = v1[i];
			else
				tmp[i][j] = v2[i];
	return tmp;
}
/// \brief ��������ֵ����λ��
/// \param A ����
/// \return Data_Pos���ݺ�λ�ýṹ��
template<typename Type> Data_Pos<Type> FindMaxandPos( const Matrix<Type> &A )
{
	int rows = A.rows();
	int cols = A.cols();
	Data_Pos<Type> tmp; 
 	tmp.Data = A[0][0];
	tmp.row = 0;
	tmp.col = 0;
	for( int i=0; i<rows; i++)
		for( int j=0; j<cols; j++)
		{
			if( A[i][j] >= tmp.Data)
			{
				tmp.Data = A[i][j];
				tmp.row = i;
				tmp.col = j;
			}
		}
	return tmp;
}
/// \brief �������Сֵ����λ��
/// \param A ����
/// \return Data_Pos���ݺ�λ�ýṹ��
template<typename Type> Data_Pos<Type> FindMinandPos( const Matrix<Type> &A )
{
	int rows = A.rows();
	int cols = A.cols();
	Data_Pos<Type> tmp; 
 	tmp.Data = A[0][0];
	tmp.row = 0;
	tmp.col = 0;
	for( int i=0; i<rows; i++)
		for( int j=0; j<cols; j++)
		{
			if( A[i][j] <= tmp.Data)
			{
				tmp.Data = A[i][j];
				tmp.row = i;
				tmp.col = j;
			}
		}
	return tmp;
}
template<typename Type>
Matrix<Type> ExchangeRows( Matrix<Type> &A, int row1, int row2 )
{
	Matrix<Type> tmp = A;
	int cols = A.cols();
	if( row1 != row2)
	{
		vector<Type> tmpv = A.getRow( row1 );
		tmp.setRow( tmp.getRow(row2), row1) ;
		tmp.setRow( tmpv, row2);
	}
	return tmp;
}
template<typename Type>
Matrix<Type> ExchangeCols( Matrix<Type> &A, int col1, int col2 )
{
	Matrix<Type> tmp = A;
	int rows = A.rows();
	if( col1 != col2)
	{
		vector<Type> tmpv = A.getColumn( col1 );
		tmp.setColumn( tmp.getColumn(col2), col1); 
		tmp.setColumn( tmpv, col2);
	}
	return tmp;
}

template<typename Type>
Matrix<Type> CopyFromMatrix( const Matrix<Type> &A, int row1, int col1, int row2, int col2 )
{
	int rows = row2-row1+1;
	int cols = col2-col1+1;
	Matrix<Type> tmp( rows, cols );
	for( int i=0; i<rows; i++)
		for(int j=0; j<cols; j++)
			tmp[i][j] = A[row1+i][col1+j];
	return tmp;
}
template<typename Type>
void Matrix<Type>::ReplaceByMatrix( const Matrix<Type> &A, int r, int c )
{
	int rn = A.rows();
	int cn = A.cols();
	int i,j;
	for( i=0; i<rn; i++)
		for( j=0; j<cn; j++)
			prow0[r+i][c+j] = A[i][j];
}
template<typename Type>
void Matrix<Type>::SetData( Type d, int rowd, int cold)
{
	prow0[rowd][cold] = d;
}

template<typename Type>
Matrix<Type> ExchangeRowData( Matrix<Type> &A , int r1, int c1, int r2, int c2, int size )
{
	Matrix<Type> tmp = A;
	int i;
	Type t;
	for( i=0; i<size; i++)
	{
		t = tmp[r1][i+c1];
		tmp[r1][i+c1] = tmp[r2][i+c2];
		tmp[r2][i+c2] = t;
	}
	return tmp;
}

template<typename Type>
void Vector2Vector( const vector<Type> &v, vector<Type> &v1)
{
	vector<Type> tmpv(v.size());
	for( int i=0; i<v.size(); i++)
		v1[i] = v[i];
}