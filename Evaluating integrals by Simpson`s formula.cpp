#define _CRT_SECURE_NO_WRARNINGS
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>

double eps = 0.0001;
int n;

double f(double x) { //первоначальная функция
	return x; //12.5
	//return 2 * sin(x) + pow(x, 2) / 10 + 2; //15.599342295
	//return x * x + x*sin(x) + sin(x); // 0-5 40.00576927
}

double F(double x) { //интеграл, посчитанный ручками
	return pow(x, 2) / 2;
	//return -2 * cos(x) + pow(x, 3) / 30 + 2 * x;
	//return pow(x, 3)/3 + sin(x) - x* cos(x) - cos(x);
}

double integral(double x_l, double x_r) {
	return F(x_r) - F(x_l);
}

double l_rectangle(double x_l, double x_r) {
	double teta = (double)1 / 3; //оценка рунге-кутты для дважды дифф интегралов
	int n = 1;
	double length = abs(x_r - x_l);
	double s_n = 0;
	double s_2n = 0;
	do { // в каждой итерации этого цикла приближенно считаем интеграл при разбиении отрезка на n частей (каждый раз n возрастает вдвое)
		double s = 0;
		double sublength = length / n;
		for (int it = 0; it < n; it++) {
			s += f(x_l + it * sublength) * sublength;
		}
		s_n = s_2n;
		s_2n = s;
		n *= 2;
	} while (teta * abs(s_2n - s_n) > eps); // ну и если по Рунге мы получили достаточную точность при таком n, то можно заканчивать

	return s_2n;
}

double c_rectangle(double x_l, double x_r) {
	double teta = (double)1 / 3;
	int n = 1;
	double length = abs(x_r - x_l);
	double s_n = 0;
	double s_2n = 0;
	do {
		double s = 0;
		double sublength = length / n;
		for (int it = 0; it < n; it++) s += f(x_l + it * sublength + 0.5 * sublength) * sublength;
		s_n = s_2n;
		s_2n = s;
		n *= 2;
	} while (teta * abs(s_2n - s_n) > eps);

	return s_2n;
}

double trapezoid(double x_l, double x_r) {


	double teta = (double)1 / 3;
	int n = 1;
	double length = abs(x_r - x_l);
	double s_n = 0;
	double s_2n = 0;
	do {
		double s = 0;
		double sublength = length / n;
		for (int it = 0; it < n; it++) s += (f(x_l + sublength * it) + f(x_l + sublength * it + sublength)) / 2 * sublength;
		s_n = s_2n;
		s_2n = s;
		n *= 2;
	} while (teta * abs(s_2n - s_n) > eps);

	return s_2n;
}

double simpson(double x_l, double x_r) {
	double teta = (double)1 / 15;
	n = 1;
	double length = abs(x_r - x_l);
	double s_n = 0;
	double s_2n = 0;
	do {
		double s = 0;
		double sublength = length / n;
		for (int it = 0; it < n; it++) s += (sublength / 6) * (f(x_l + sublength * it) + f(x_l + sublength * it + sublength / 2) * 4 + f(x_l + sublength * it + sublength));
		s_n = s_2n;
		s_2n = s;
		n *= 2;
	} while (teta * abs(s_2n - s_n) > eps);

	return s_2n;
}

int main() {
	double a = 0, b = 3;
	double tmp;
	double tmpint;
	tmpint = integral(a, b);
	std::cout << "Definite integral:\t" << tmpint << "\t" << 0 << "\n";
	tmp = l_rectangle(a, b);
	std::cout << "Left rectangles:\t" << tmp << "\t" << abs(tmp - tmpint) << "\n";
	tmp = c_rectangle(a, b);
	std::cout << "Central rectangles:\t" << tmp << "\t" << abs(tmp - tmpint) << "\n";
	tmp = trapezoid(a, b);
	std::cout << "Trapezoids:\t\t" << tmp << "\t" << abs(tmp - tmpint) <<  "\n";
	tmp = simpson(a, b);
	std::cout << "Simpson:\t\t" << tmp << "\t" << abs(tmp - tmpint) <<  "\n";
}
