#include<iostream>
#include"matrix.hpp"
#include<stdlib.h>
#include<time.h>
#include<ctime>
void print_test();
void mult_test(int size);
Matrix generate_matrix(int size, int seed=0);
double time_int_function( void f(int), int int_value=10);
void test_givens_rotation();
void transp_test();
void rot_decomposition_test(int size);


int main()
{
//	time_int_function(mult_test,1000);
//	test_givens_rotation();
//	print_test();
//	unsigned int start_time = clock();
///	mult_test(500);
//	unsigned int time=clock()-start_time;
//	std::cout<<"time is "<<(float)time/CLOCKS_PER_SEC<<"s"<<std::endl;
//	generate_matrix(1000);	
//	transp_test();
//	test_givens_rotation();
//	rot_decomposition_test();
	time_int_function(rot_decomposition_test,20);
	return 0;
}

void rot_decomposition_test(int size)
{
	std::pair<Matrix,Matrix> QR;
	Matrix A=generate_matrix(size);
	std::cout<<"Original matrix is"<<std::endl;
	A.print();
	std::cout<<"Now let's decompose it"<<std::endl;
	QR=A.QR_rot_decomposition();
	std::cout<<"Q matrix is"<<std::endl;
	QR.first.print();
	std::cout<<"R matrix is"<<std::endl;
	QR.second.print();
	std::cout<<"Let's have some test::"<<std::endl;
	std::cout<<"1.Check Q*R:"<<std::endl;
	(QR.first*QR.second).print();
	std::cout<<"2.Check Q is ortogonal"<<std::endl;
	std::cout<<"Q^T*Q"<<std::endl;
	(QR.first*QR.first.transp()).print();
}

void transp_test()
{
	Matrix A=generate_matrix(10);
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
	((Givens_rotation(0,1,a[0][pos_x],a[pos_y][pos_x],size))*A).print();
}

double time_int_function(void f(int),int int_value)
{
	int start_time=clock();
	f(int_value);
	int time=clock()-start_time;
	std::cout<<"time is"<<(double)time/CLOCKS_PER_SEC<<"s"<<std::endl;
}
//А как много времени уходит на генерацию и как много времени уходит на само умножение ?
//Эти процессы надо разделить
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

void mult_test(int size)
{
	Matrix A = generate_matrix(size);
	Matrix B = generate_matrix(size,100);
//	std::cout<<"Multiplying test"<<std::endl;
//	std::cout<<"Matrix A is"<<std::endl;
//	A.print();
//	std::cout<<"Matrix B is"<<std::endl;
//	B.print();
//	std::cout<<"their multiplication is C:"<<std::endl;
	Matrix C;
//	std::cout<<"HEHE"<<std::endl;
	C=A*B;
//	C.print();	
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
