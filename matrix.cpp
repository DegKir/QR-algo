#include"matrix.hpp"
#include<iostream>
#include<utility>
#include<vector>
#include<cmath>

#define SIZE 100
#define EPSILON 1E-10

Matrix::Matrix(){}

Matrix::Matrix(std::vector<std::vector<double>> coefficients) : _coefficients(coefficients){}

std::vector<std::vector<double>> Matrix::get_coefficients() const
{
	return _coefficients;
}

int Matrix::size() const
{
	return _coefficients.size();
}

bool Matrix::is_upper_triangle() const
{
	int size = this->size();
	for(int i = 0; i < size; ++i)
	{
		for(int j = 0; j< i ; ++j)
		{
			if(_coefficients[i][j]!=0)
				return false;
		}
	}
	return true;
}

double Matrix::get_norm() const
{
	return 0;
}

void Matrix::print() const
{
	std::cout.precision(5);
	int size = this->size();
	for(int i = 0;i<size;i++)
	{
		for(int j = 0; j<size;j++)
		{
			std::cout<<_coefficients[i][j]<<" ";
		}
		std::cout<<std::endl;
	}
}

Matrix Matrix::transp()
{
	std::vector<std::vector<double>> coef=_coefficients;
	double temp;
	int size = _coefficients.size();
	for(int i = 0; i < size ; i++)
		for(int j = i; j < size ; j++)
		{
			temp=coef[i][j];
			coef[i][j]=coef[j][i];	
			coef[j][i]=temp;
		}
	return Matrix(coef);
}

//Нужно дописать динамическюу память
Matrix operator*(const Matrix &A,const Matrix &B)
{
	int size = A.size();
	std::vector<std::vector<double>> ccoef(size);
	std::vector<double> row(size);
	std::vector<std::vector<double>> bcoef=B._coefficients;
	std::vector<std::vector<double>> acoef=A._coefficients;

//	Проверим, что размеры матриц совпадают
	if(A.size()!=B.size())
		std::cout<<"ERROR"<<std::endl;

//	Заполняем итоговую матрицу пустыми векторами
	for(int i = 0; i < A.size();i++)
		ccoef[i]=row;

//	Создаем и заполняем статические массивы
	double a[SIZE][SIZE];
	double b[SIZE][SIZE];
	double c[SIZE][SIZE];
	for (int i = 0; i < size; i++)
	{
		for(int j = 0; j< size;j++)
		{
			a[i][j]=acoef[i][j];
			b[i][j]=bcoef[i][j];
		}
	}

	register double sum=0;
	for(int i = 0; i < size; ++i)
	{
		for(int j = 0; j < size; ++j)
		{
			sum=0;
			c[i][j]=0;
			for(int k = 0; k< size; ++k)
			{
				sum+= a[i][k]*b[k][j];
			}
			if(sum*sum<EPSILON)
				sum=0;
			c[i][j]=sum;
		}
	}
//
//	register double sum=0;
//	int column_size=A.column_size();
//	int A_row_size=A.row_size();
//	int row_size=B.row_size();
//
//	for(int i = 0; i < column_size;i++)
//	{
//		row=A._coefficients[i];
//		for(int j = 0; j < row_size;j++)
//		{
//			for(int k = 0; k < A_row_size;k++)
//			{
//				sum+=row[k]*bcoef[k][j];
//			}
//			if(sum*sum<EPSILON)
//				sum=0;
//			ccoef[i][j]=sum;
//			sum=0;
//		}
//	}


	for(int i = 0; i < size;i++)
	{
		for(int j = 0; j< size;j++)
		{
			ccoef[i][j]=c[i][j];
		}
	}	
	Matrix result(ccoef);
	return result;
}

Matrix::~Matrix(){}

Matrix Givens_rotation(int pivot_position, int delete_position, double pivot_value, double delete_value,int size)
{
	double a1=pivot_value;
	double a2=delete_value;

	//1.Создаем единичную матрицу
	std::vector<double> row(size);
	std::vector<std::vector<double>> E(size);
	for(int i = 0; i < size; i++)
		E[i]=row;
	for(int i = 0; i < size; i++)
		E[i][i]=1;

	//2.Вычисляем синус и косинус
	double rot_cos=a1/sqrt(a1*a1+a2*a2);
	double rot_sin=-a2/sqrt(a1*a1+a2*a2);
	E[pivot_position][pivot_position]=rot_cos;
	E[delete_position][delete_position]=rot_cos;
	E[pivot_position][delete_position]=-rot_sin;
	E[delete_position][pivot_position]=rot_sin;
	return Matrix(E);
}


std::pair<Matrix,Matrix> Matrix::QR_rot_decomposition()
{
	Matrix Givens;
	Matrix R=Matrix(_coefficients);
	std::vector<std::vector<double>> coefficients=_coefficients;

	std::vector<Matrix> rotation_matrixes;
	rotation_matrixes.reserve(100);
	int k = 1;

	//Проходимся по нижне-треугольной матрице
	for(int i=0;i<coefficients.size();i++)
	{
		for(int j = 0;j<i;j++)
		{
//			std::cout<<"deleting elemnt at i="<<i<<"and at j="<<j<<std::endl;
			//Получаем матрицу Гивенса для даного элемента, выбрав опорным первый(пока что).
			int pivot_position=j;////!!!!!
			int delete_position=i;
			double pivot_value=coefficients[j][j];
			double delete_value=coefficients[i][j];
			int size = coefficients.size();
			Givens=Givens_rotation(pivot_position, delete_position, pivot_value, delete_value, size);	
			//Умножаем рабочую матрицу на матрицу Гивенса, саму матрицу Гивенса добавляем в вектор
			rotation_matrixes.push_back(Givens);
			R=Givens*R;
//			std::cout<<"G matrix is"<<std::endl;
//			Givens.print();
//			std::cout<<"R matrix is now"<<std::endl;
//			R.print();
//			std::cout<<std::endl;
			coefficients=R.get_coefficients();
		}
	}
	Matrix D;
	//Q=G_1^T*G_2^T*...*G_n^T*
	Matrix Q=rotation_matrixes[0].transp();
	for(int i = 1; i < rotation_matrixes.size() ; i++)
	{
		D=rotation_matrixes[i].transp();
		Q=Q*D;		
	}
	return std::pair<Matrix,Matrix>{Q,R};
}

Matrix QR_algorithm(Matrix A, int max_iterations)
{
	Matrix B=A;
	Matrix Q,R;
	std::pair<Matrix,Matrix> QR;
	int k = 0;
	while(!B.is_upper_triangle())
	{
		k++;
		if(k>max_iterations)
			break;
		std::cout<<k<<std::endl;
		QR=B.QR_rot_decomposition();
		Q=QR.first;
		R=QR.second;
		B=R*Q;
	}
	return B;
}

Matrix get_Hessenberg(Matrix A)
{
	std::vector<std::vector<double>> coefficients=A.get_coefficients();
	Matrix R=A;
	Matrix Givens;
//	for(int i = 0; i<coefficients.size();i++)
//	{
//		for(int j = 0; j<i-1;j++)
//		{
//			int pivot_position=j;
//			int delete_position=i;
//			double pivot_value=coefficients[j][j];
//			double delete_value=coefficients[i][j];
//			int size = coefficients.size();
//			Givens=Givens_rotation(pivot_position, delete_position, pivot_value, delete_value, size);	
//			std::cout<<"+++++++++++++++++++"<<std::endl;
//			(Givens*R).print();
//			std::cout<<"+++++++++++++++++++"<<std::endl;
//			R=Givens*R;
//			R=R*Givens.transp();
//			coefficients=R.get_coefficients();
//		}
//	}
	Matrix Temp;
	for(int j = 0; j<coefficients.size();j++)
	{
//		std::cout<<"J="<<j<<std::endl;
		for(int i = coefficients.size()-1;i>j+1;i--)
		{
//			std::cout<<"I="<<i<<std::endl;
//			std::cout<<"NEW ITERATION"<<std::endl;
			int pivot_position=i-1;
			int delete_position=i;
			double pivot_value=coefficients[i-1][j];
			double delete_value=coefficients[i][j];
			int size=coefficients.size();
//			std::cout<<"Pivot element is "<<pivot_value<<std::endl;
//			std::cout<<"Delete element is "<<delete_value<<std::endl;
			Givens=Givens_rotation(pivot_position, delete_position, pivot_value, delete_value,size);
//			std::cout<<"+++++++++++++++++++"<<std::endl;
//			std::cout<<"+++++++++++++++++++"<<std::endl;
//			(Givens*R).print();
//			std::cout<<"+++++++++++++++++++"<<std::endl;
			R=Givens*R;
			Temp=Givens.transp();
			R=R*Temp;
//			R.print();
			coefficients=R.get_coefficients();
		}
	}
	return Matrix(R);
}
