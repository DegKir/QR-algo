#include<utility>
#include<vector>
//Тестовая функция для проверки коректности подключения файлов
void matrix_hello();

//В этом классе реализуется объект вектора из линейной алгебры
class Vector
{
public:
	friend Vector operator*(Vector,Vector);					//Скалярное умножение векторов
	friend std::vector<Vector> ortogonalization(std::vector<Vector>);	//Ортогонализация Грамма-Шмидта
private:
	std::vector<double> _coefficients;					//Вектор коэффициентов
};


//В этом классе реалиузется объект матрицы
class Matrix
{
public:
	Matrix();
	Matrix(std::vector<std::vector<double>>);
	Matrix transp();
	std::vector<std::vector<double>> get_coefficients();
	int column_size();
	int row_size();
	double get_norm();				//Возвращает 2-норму матрицы
	void print();					//Печатает матрицу
	std::pair<Matrix,Matrix> QR_rot_decomposition();//QR разложение, полученное вращениями
	~Matrix();
	friend Matrix operator*(Matrix,Matrix); 	//Возвращает произведение матриц
	
private:
	
	std::vector<std::vector<double>> _coefficients;	//Вектор векторов коэффициентов	
};

//Функция возвращает Матрицу Гивенса
Matrix Givens_rotation(int pivot_position, int delete_position, double pivot_value, double delete_value,int size);


//Функция возвращает результат QR алгоритма от исходной матрицы
Matrix QR_algorithm(Matrix);
