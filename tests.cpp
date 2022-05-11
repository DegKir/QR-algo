#include<iostream>
#include"matrix.hpp"
#include<stdlib.h>
#include<time.h>
#include<ctime>
void print_test();
void mult_test(Matrix&,Matrix&);
Matrix generate_matrix(int size, int seed=0);
Matrix generate_sym_matrix(int size, int seed=0);
double time_matrix_function(Matrix&, void f(Matrix&));
double time_matrix_function(Matrix&,Matrix&, void f(Matrix&,Matrix&));

void test_givens_rotation();
void transp_test();
void rot_decomposition_test(int size);
void RQ_algorithm_test(Matrix&,int);
void Hessenberg_test(int size);
void RQ_Hessenberg_test(Matrix&,int);
void cut_test(Matrix&);
int main()
{
	int size = 12;
//	rot_decomposition_test(size);
	Matrix A = generate_sym_matrix(size,10);
//	cut_test(A);
	Matrix B = generate_sym_matrix(size,20);
//	time_matrix_function(A,B,mult_test);
//	rot_decomposition_test(10);
//	RQ_algorithm_test(A,10000);



	RQ_Hessenberg_test(A,10000);
	return 0;
}

void cut_test(Matrix& A)
{
	int size=A.size();
	for(int i = 0; i< size; i++)
	{
		A.print();
		std::cout<<"++++++++++++++"<<std::endl;
		A=A.cut();
	}
}

void RQ_Hessenberg_test(Matrix &A,int ticks)
{
	std::cout<<"Original matrix is"<<std::endl;
	A.print();
	Matrix B=get_Hessenberg(A);
	std::cout<<"Hessenberg matrix is"<<std::endl;
	B.print();
	std::cout<<"Result of Hessenberg QR is matrix:"<<std::endl;
	B=QR_algorithm(A,ticks);
	B.print();
}

void Hessenberg_test(int size)
{
	Matrix A=generate_sym_matrix(size);
	std::cout<<"Original matrix is"<<std::endl;
	A.print();
	std::cout<<"Now let's transorm it to the Hessenbger matrix"<<std::endl;
	get_Hessenberg(A).print();
}

void RQ_algorithm_test(Matrix& A,int ticks)
{
	std::cout<<"Originate matrix is"<<std::endl;
	A.print();
	std::cout<<"Result of simple QR is matrix:"<<std::endl;
	Matrix B = QR_algorithm(A,ticks);
	B.print();
}

void rot_decomposition_test(int size)
{
	std::pair<Matrix,Matrix> QR;
	Matrix A=generate_sym_matrix(size);
	Matrix Q,R,QT;
	int start_time=clock();

	std::cout<<"Original matrix is"<<std::endl;
	A.print();
	std::cout<<"Now let's decompose it"<<std::endl;
	QR=A.QR_rot_decomposition();
	Q=QR.first;
	R=QR.second;
	std::cout<<"Q matrix is"<<std::endl;
	Q.print();
	std::cout<<"R matrix is"<<std::endl;
	R.print();
	std::cout<<"Let's have some test::"<<std::endl;
	std::cout<<"1.Check Q*R:"<<std::endl;
	(Q*R).print();
	std::cout<<"2.Check Q is ortogonal"<<std::endl;
	std::cout<<"Q^T*Q"<<std::endl;
	QT=Q.transp();
	(Q*QT).print();
	int time=clock()-start_time;
	std::cout<<"time is"<<(double)time/CLOCKS_PER_SEC<<"s"<<std::endl;

//	(QR.first*QR.first.transp()).print();
}

void transp_test()
{
	Matrix A=generate_matrix(4);
	A.print();
	std::cout<<"now transporting"<<std::endl;
	(A.transp()).print();
}

void test_givens_rotation()
{
	int size = 3;
	int pos_x=1;
	int pos_y=1;
	Matrix A=generate_matrix(size);
	std::vector<std::vector<double>> a=A.get_coefficients();
	A.print();
	std::cout<<"Now let's delete at ("<<pos_x<<","<<pos_y<<") elemnt via rotation"<<std::endl;
	std::cout<<"G matrix is"<<std::endl;
	(Givens_rotation(0,1,a[0][pos_x],a[pos_y][pos_x],size)).print();
	std::cout<<"Rotated matrix is"<<std::endl;
//	((Givens_rotation(0,1,a[0][pos_x],a[pos_y][pos_x],size))*A).print();
}

double time_matrix_function(Matrix &A, void f(Matrix&))
{
	int start_time=clock();
	f(A);
	int time=clock()-start_time;
	std::cout<<"time is"<<(double)time/CLOCKS_PER_SEC<<"s"<<std::endl;
}

double time_matrix_function(Matrix &A,Matrix &B, void f(Matrix&, Matrix&))
{
	std::cout<<"Time function activated"<<std::endl;
	int start_time=clock();
	std::cout<<"function has begun"<<std::endl;
	f(A,B);
	int time=clock()-start_time;
	std::cout<<"time is"<<(double)time/CLOCKS_PER_SEC<<"s"<<std::endl;
}


Matrix generate_matrix(int size,int seed)
{
	std::vector<double> row(size);
	std::vector<std::vector<double>> mat(size);
	for(int i = 0 ; i < size; i++)
		mat[i]=row;
	srand(time(NULL)+seed);
	for(int i = 0; i < size; i++)
		for(int j = 0; j < size; j++)
			mat[i][j]=rand() % 10;
	return Matrix(mat);
}

Matrix generate_sym_matrix(int size,int seed)
{
	std::vector<double> row(size);
	std::vector<std::vector<double>> mat(size);
	for(int i = 0;i < size; i++)
		mat[i]=row;
	srand(time(NULL)+seed);
	for(int i = 0; i< size;i++)
	{
		for(int j = i; j< size;j++)
		{
			int k = rand() % 10;
			mat[i][j]=k;
			mat[j][i]=k;
		}
	}
	return Matrix(mat);
}

void mult_test(Matrix &A, Matrix &B)
{
	Matrix C;
	std::cout<<"Mult test activated"<<std::endl;
	C=A*B;
	std::cout<<"Mult test ended\\result is : "<<std::endl;
	C.print();
}
//
void print_test()
{
	std::vector<double> A{1,2,3};
	std::vector<double> B{4,5,6};
	std::vector<double> C{7,8,9};
	Matrix M(std::vector<std::vector<double>> {A,B,C});
	M.print();
}
