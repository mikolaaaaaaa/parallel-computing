#include <iostream>
#include <mpi.h>

using namespace std;

int** result_matrix;

int** get_matrix(int size) {
	int **matrix = new int*[size];
	for (int i = 0; i < size; i++) {
		matrix[i] = new int[size];
	}
	return matrix;
}

void fill_matrix(int** matrix, int size) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			matrix[i][j] = 2;
		}
	}
}

void fill_vector(int* vect, int size) {
	for (int i = 0; i < size; i++) {
		vect[i] = 1;
	}
}

int* get_vector(int size) {
	int* vect = new int[size];
	return vect;
}

void print_matrix(int **matrix, int size) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
}

int multiply_vector(int *a, int *b, int n) {

	int result = 0;
	for (int i = 0; i < n; i++) {
		result += a[i] * b[i];
	}
	return result;

}

void multiply_matrix(int **res_matrix, int **arr, int from, int to, int n) {

	int c = 0;
	for (int i = from; i < to; i++) {
		for (int ii = 0; ii < n; ii++) {
			for (int j = 0; j < n; j++) {
				c += arr[i][j] * arr[j][ii];
			}
			res_matrix[i][ii] = c;
			c = 0;
		}
	}

}

void multiply_matrix_by_number(int **matrix, int number, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			matrix[i][j] *= number;
		}
	}
}



int main(int argc, char **argv)
{
	int myrank, size, i, m;
	int C = 0; // результат умножения векторов
	int* A = 0, * B = 0; // входные векторы
	int N; // размер векторов и массива
	int* a = 0, * b = 0; // части векторов для каждого процесса

	int** Matrix = 0; // входящая матрица 
	int** matrix = 0; // результирующая матрица
	
	 MPI_Init(&argc, &argv);
  	 MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	 MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	if (myrank == 0)
	{
		cout << "Vvedite N"; cin >> N;
	}

	MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (myrank == 0) // выделение памяди
	{
		A = get_vector(N);
		B = get_vector(N); 
		Matrix = get_matrix(N);
		result_matrix = get_matrix(N);
	}
	
	// определения количества элементов в каждом процессе
	m = N / size;
	a = get_vector(m);
	b = get_vector(m);
	matrix = get_matrix(N);
	
	if (myrank == 0)
		fill_vector(A, N); // заполнение исходной матрицы А
    	fill_vector(B, N);
		fill_matrix(Matrix, N);
		MPI_Scatter(A, m, MPI_INT, a, m, MPI_INT, 0,
		MPI_COMM_WORLD);
    	MPI_Scatter(B, m, MPI_INT, b, m, MPI_INT, 0,
		MPI_COMM_WORLD);
		MPI_Bcast(Matrix, 1, MPI_INT, 0, MPI_COMM_WORLD);
		// рассылка частей матрицы по процессам
	
		int c = multiply_vector(a, b, m); // основная функция обработки
		multiply_matrix(result_matrix, Matrix, myrank * m, (myrank + 1) * m, N);
		MPI_Reduce(&c, &C, 1, MPI_INT, MPI_SUM,
		0, MPI_COMM_WORLD);

	// сбор результатов в корневой процесс
		if (myrank == 0) {
			multiply_matrix_by_number(result_matrix, C, N);
			print_matrix(result_matrix, N);
		//	cout << C; // вывод результата
		}
	MPI_Finalize();
}


