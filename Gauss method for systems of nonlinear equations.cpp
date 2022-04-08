#include <math.h>
#include <conio.h>
#include <iostream>

// Количество уравнений
const int N = 2;

double f(int num, double x, double y) {
	switch (num) {
	case 0:
		return sin(x+y);
	case 1:
		return x;
	default:
		std::cout << "equation number " << num << "not found";
		exit(0);
	}
}

double derivative_x(int num, double x, double y) {
	double dx = 0.00001;
	return (f(num, x + dx, y) - f(num, x, y))/dx;
}

double derivative_y(int num, double x, double y) {
	double dy = 0.00001;
	return (f(num, x, y + dy) - f(num, x, y)) / dy;
}

double jacobi(double x, double y) {
	return derivative_x(0, x, y) * derivative_y(1, x, y) - derivative_y(0, x, y) * derivative_x(1, x, y);
}

double norm(double x, double y) {
	return sqrt(x * x + y * y);
}

double norm(double* x) {
	double sum = 0;
	for (int i = 0; i < N; i++) {
		sum += x[i]*x[i];
	}
	return sqrt(sum);
}

double phi(int num, double x, double y) {
	switch (num) {
	case 0:
		return f(num, x, y) + x;
	case 1:
		return f(num, x, y) + y;
	default:
		break;
	}
}
void GaussSolve(double** a, double* b, double* x) {
// Решение СЛУ методом Гаусса с выбором главного элемента по столбцу
// На входе: матрица коэффициентов, вектор правых частей, вектор решения
// На выходе: заполненный вектор решения
	const double eps = 0.0001;
	double tmpValue;
	int i, j, k, z;
	double dblLeadElement;

	// Прямой ход
	for (i = 0; i < N; i++) {

		dblLeadElement = a[i][i];

		// Находим строку, в которой элемент, стоящий под ведущим - наибольший
		double tmpMax = dblLeadElement;
		int tmpMaxNumber = i;

		for (z = i; z < N; z++) {
			if (a[z][i] > tmpMax) { tmpMax = a[z][i]; tmpMaxNumber = z; }
		}

		// Меняем местами i-ю строку и строку tmpMaxNumber
		for (z = i; z < N; z++) {
			tmpValue = a[i][z];
			a[i][z] = a[tmpMaxNumber][z];
			a[tmpMaxNumber][z] = tmpValue;

		}
		tmpValue = b[i];
		b[i] = b[tmpMaxNumber];
		b[tmpMaxNumber] = tmpValue;

		dblLeadElement = tmpMax;

		for (j = i; j < N; j++) {
			a[i][j] /= dblLeadElement;
		}
		b[i] /= dblLeadElement;

		for (k = i + 1; k < N; k++) {

			double dblToDivide = a[k][i] / a[i][i];
			for (z = i; z < N; z++) {
				a[k][z] -= a[i][z] * dblToDivide;
			}

			b[k] -= b[i] * dblToDivide;
		}
	}

	// Обратный ход
	x[N - 1] = b[N - 1];

	for (k = N - 2; k >= 0; k--) {
		double sum = b[k];

		for (j = k + 1; j < N; j++) {
			sum -= a[k][j] * x[j];
		}
		x[k] = sum;
	}
}

void Newton() {
	printf("%s", "Решение системы нелинейных уравнений методом Ньютона\n");

	int count = 0;

	// Точность
	double eps = 0.00001;											//Метод Ньютона
																	//Ищется ближайшее решение к начальному приближению
	// Матрица														//Сложность заключается в вычислении обратной к Якоби
	double** matrix = new double*[N];								//Поэтому часто используют вид x(k+1) = x(k) + dx(k)
	for (int i = 0; i < N; i++) {									//dx(k)- Находится из решения f(x(k)) + J(x(k))*dx(x(k)) = 0
		matrix[i] = new double[N];									//f(x(k)) - вектор-столбец функции
	}																//J(x(k)) - матрица Якоби от К-того приближения
																	//dx(x(k)) - типа приближения к вектору-ответу
	// Вектор правых частей											//((У меня в решении через dx используется метод Гаусса))
	double* b = new double[N];										//
																	//Для метода требуется:
	// Вектор решения												//Невырожденная матрица Якоби
	double* x = new double[N];										//Дифференцируемые функции в системе ур-ий
	double* xk = new double[N];										//
	double* p = new double[N];										//Критерий остановки:
																	//fabs(X(k+1) - X(k)) < eps, eps - требуемая точность
	// Выбираем начальное приближение								//
	xk[0] = 1; xk[1] = 1;											//КОЛИЧЕСТВО ИТЕРАЦИЙ ЗАВИСИТ ОТ ТОЧНОСТИ КВАДРАТИЧНО//сходится квадратично
																	//
	do                                                              //
	{
		count++;

		for (int i = 0; i < N; i++) {
			x[i] = xk[i];
		}

		// Заполняем матрицу A
		matrix[0][0] = derivative_x(0,x[0],x[1]);
		matrix[0][1] = derivative_y(0,x[0], x[1]);
		matrix[1][0] = derivative_x(1, x[0], x[1]); 
		matrix[1][1] = derivative_y(1, x[0], x[1]);

		// Заполняем вектор правых частей
		b[0] = -f(0, x[0], x[1]); 
		b[1] = -f(1, x[0], x[1]);

		GaussSolve(matrix, b, p);
		// Получаем следующее приближение,
		// складывая предыдущее и вновь полученное
		for (int i = 0; i < N; i++) {
			xk[i] = x[i] + p[i];
		}

	} while (fabs(norm(x[0], x[1]) - norm(xk[0],x[1])) >= eps);
	printf("%s", "Вектор решения: [ ");
	for (int i = 0; i < N; i++) {
		printf("%f ", xk[i]);
	}
	printf("%s", "]\n");
	printf("Решение найдено за %i итераций\n", count);
}

void simple_iterations() {
	std::cout << "Метод простых итераций" << std::endl;
	int count = 0;
	double eps = 0.00001;
	double* xk = new double[N];
	double* x = new double[N];
	for (int i = 0; i < N; i++) {
		x[i] = 1;
	}
	if (fabs(jacobi(x[0],x[1])) <= 1) {
		std::cout << "Метод сходится" << std::endl;
	}else{
		std::cout << "Метод не сходится" << std::endl;
		exit(0);
	}
	do {
		count++;
		for (int i = 0; i < N; i++) {
			xk[i] = x[i];
		}
		x[0] = phi(0, xk[0], xk[1]);
		x[1] = phi(1, x[0], xk[1]);
	} while (fabs(norm(x) - norm(xk)) > eps);

	std::cout << "Вектор - решение: [" << x[0];
	for (int i = 1; i < N; i++) {
		std::cout << ", " << x[i];
	}
	std::cout <<"]" << std::endl;

	std::cout << "count: " << count << std::endl;
}


int main() {
	setlocale(LC_ALL, "Rus");
	Newton();
	std::cout << std::endl;
	simple_iterations();
	return 0;
}

//fabs(norm(x) - norm(xk)
