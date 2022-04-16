#include"matrix.hpp"
#include<iostream>
#include<utility>
#include<vector>
#include<cmath>

#define EPSILON 1E-10

void matrix_hello()
{
	std::cout<<"Matrix hello"<<std::endl;
}

Matrix::Matrix(){}

Matrix::Matrix(std::vector<std::vector<double>> coefficients) : _coefficients(coefficients){}

std::vector<std::vector<double>> Matrix::get_coefficients()
{
	return _coefficients;
}

int Matrix::column_size()
{
	return _coefficients.size();
}

int Matrix::row_size()
{
	return _coefficients[0].size();
}

double Matrix::get_norm()
{
	return 0;
}

void Matrix::print()
{
	std::cout.precision(3);
	for(int i = 0;i<(this->column_size());i++)
	{
		for(int j = 0; j<(this->row_size());j++)
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
	for(int i = 0; i < _coefficients.size() ; i++)
		for(int j = i; j < _coefficients.size() ; j++)
		{
			temp=coef[i][j];
			coef[i][j]=coef[j][i];	
			coef[j][i]=temp;
		}
	return Matrix(coef);
}

Matrix operator*(Matrix A,Matrix B)
{
	std::vector<std::vector<double>> ccoef;
	std::vector<double> row(B.row_size());
	ccoef.reserve(A.column_size());
	std::vector<std::vector<double>> bcoef=B._coefficients;

	//Проверим, что размеры матриц совпадают
	if(A.row_size()!=B.column_size())
		std::cout<<"ERROR"<<std::endl;

	for(int i = 0; i < A.column_size();i++)
		ccoef.push_back(row);

	Matrix C(ccoef);	


	double sum=0;
	for(int i = 0; i < C.column_size();i++)
	{
		row=A._coefficients[i];
		for(int j = 0; j < C.row_size();j++)
		{
			for(int k = 0; k < A.row_size();k++)
			{
				sum+=row[k]*bcoef[k][j];
			}
			ccoef[i][j]=sum;
			sum=0;
		}
	}
	
	for(int i = 0; i < C.column_size();i++)
		for(int j = 0;j<C.row_size();j++)
			if(ccoef[i][j]*ccoef[i][j]<EPSILON)
				ccoef[i][j]=0;

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
	//Проверяем, что матрица на входе квартная
	if(this->column_size()!=this->row_size())
		std::cout<<"ERROR ERROR ERROR"<<std::endl;

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
	//Q=G_1^T*G_2^T*...*G_n^T*
	Matrix Q=rotation_matrixes[0].transp();
	for(int i = 1; i < rotation_matrixes.size() ; i++)
	{
		Q=Q*rotation_matrixes[i].transp();		
	}
	return std::pair<Matrix,Matrix>{Q,R};
}
