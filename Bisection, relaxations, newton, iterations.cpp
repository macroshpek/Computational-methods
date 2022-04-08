#define _CRT_SECURE_NO_WRARNINGS
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>

double bisection(double x_l, double x_r);
double simple_iterations(double x_l, double x_r);
double newton(double x_l, double x_r);
double relax(double x_l, double x_r);
double myMax(double x_l, double x_r);
double myMin(double x_l, double x_r);
double deriv(double x);
double koef(double x_l, double x_r);

using namespace std;
double eps = 0.0001;
int counter;

double f(double x) {
	return 2 * sin(x) + pow(x, 2) / 10 - 1;
	//return pow(x, 5) + 2 * x + 10 * pow(x, 2);
}

int main() {
	double x_l = -1, x_r = 3;
	double x_max = myMax(x_l, x_r);
	double x_min = myMin(x_l, x_r);

	cout << "Bisection: " << bisection(x_l, x_r) << " in " << counter << " iterations \n";
	cout << "Simple iterations: " << simple_iterations(x_l, x_r) << " in " << counter << " iterations \n";
	cout << "Newton: " << newton(x_l, x_r) << " in " << counter << " iterations \n";
	cout << "Relax: " << relax(x_l, x_r) << " in " << counter << " iterations \n";
	cout << "x_max from [" << x_l << ", " << x_r << "] = " << x_max << ", f(x_max) = " << f(x_max) << "\n";
	cout << "x_min from [" << x_l << ", " << x_r << "] = " << x_min << ", f(x_min) = " << f(x_min) << "\n";
	system("pause");
}

double deriv(double x)
{
	const double mtDx = 1.0e-6;

	double x1 = f(x + mtDx);
	double x0 = f(x);

	return (x1 - x0) / mtDx;
}

double bisection(double x_l, double x_r) {
	double x_m;
	counter = 0;
	do {
		if (abs(f(x_l)) < eps) return x_l;
		if (abs(f(x_r)) < eps) return x_r;
		x_m = (x_l + x_r) / 2;
		if (f(x_l) * f(x_m) < 0)
			x_r = x_m;
		else
			x_l = x_m;
		counter++;
	} while (x_r - x_l > eps);
	return x_m;
}

double simple_iterations(double x_l, double x_r) {
	counter = 0;
	double k = 2.0 / (myMax(x_l, x_r) + 0.25);
tryAgain:
	double x, xNew = (x_l + x_r) / 2;
	do {
		x = xNew;
		xNew = x - k * f(x);
		counter++;
	} while (abs(x - xNew) > eps);

	if (xNew<x_l || xNew>x_r) {
		k = -k;
		goto tryAgain;
	}
	return xNew;
}

double newton(double x_l, double x_r) {
	counter = 0;
	double x, xNew = (x_l + x_r) / 2 + 0.25;
	do {
		x = xNew;
		xNew = x - f(x) / deriv(x);
		//cout << deriv(x) << "\n";
		counter++;
	} while (abs(x - xNew) > eps);
	return xNew;
}

double relax(double x_l, double x_r) {
	counter = 0;
	double k = 2 / koef(x_l, x_r);
	//cout << myMax(x_l, x_r) << "\n";
	//cout << myMin(x_l, x_r) << "\n";
	double x, xNew = (x_l + x_r) / 2 + 1000 * eps;
	do {
		x = xNew;
		xNew = x - k * f(x);
		counter++;
	} while (abs(x - xNew) > eps);
	return xNew;
}

double koef(double x_l, double x_r) {
	double* arr = new double[1000];
	double length = x_r - x_l;
	for (int i = 0; i < 1000; i++) {
		arr[i] = deriv(x_l + i * length / 1000);
	}
	double min = arr[0], max = arr[0];
	for (int i = 0; i < 1000; i++) {
		if (max < arr[i]) max = arr[i];
		if (min > arr[i]) min = arr[i];
	}
	//cout << min << "\n";
	//cout << max << "\n";
	return max;
}

double myMin(double x_l, double x_r) {
	double* sols = new double[10];
	counter = 0;
	time_t systime;
	time(&systime);
	srand((unsigned int)systime);
	double alpha = 0.9;
	double length = x_r - x_l;
	double x = (x_l + x_r) / 2;
	double L = f(x);
	for (int i = 0; i < 10; i++) {
		for (double T = 80; T > 0.00008; T *= alpha)
		{
			for (int i = 0; i < 500; i++)
			{
			tryAgain:
				double xNew = x + (length / 10) * ((rand() / (double)RAND_MAX) * 2 - 1);
				counter++;
				if (xNew < x_l || xNew > x_r) goto tryAgain;
				double LNew = f(xNew);

				if (LNew < L || (rand() / (double)RAND_MAX) <= exp(-(LNew - L) / T))
				{
					L = LNew;
					x = xNew;
				}
			}
		}
		sols[i] = x;
	}
	x = 0;
	for (int i = 0; i < 10; i++) x += sols[i];
	x = x / 10;
	sols[0] = x_l;
	sols[1] = x_r;
	for (int i = 0; i < 2; i++) {
		if (f(x) > f(sols[i])) x = sols[i];
	}
	//cout << counter << "\n";
	return x;
}

double myMax(double x_l, double x_r) {
	double* sols = new double[10];
	counter = 0;
	time_t systime;
	time(&systime);
	srand((unsigned int)systime);
	double alpha = 0.9;
	double length = x_r - x_l;
	double x = (x_l + x_r) / 2;
	double L = -f(x);
	for (int i = 0; i < 10; i++) {
		for (double T = 80; T > 0.00008; T *= alpha)
		{
			for (int i = 0; i < 500; i++)
			{
			tryAgain:
				double xNew = x + (length / 10) * ((rand() / (double)RAND_MAX) * 2 - 1);
				counter++;
				if (xNew < x_l || xNew > x_r) goto tryAgain;
				double LNew = -f(xNew);

				if (LNew < L || (rand() / (double)RAND_MAX) <= exp(-(LNew - L) / T))
				{
					L = LNew;
					x = xNew;
				}
			}
		}
		sols[i] = x;
	}
	x = 0;
	for (int i = 0; i < 10; i++) x += sols[i];
	x = x / 10;
	sols[0] = x_l;
	sols[1] = x_r;
	for (int i = 0; i < 2; i++) {
		if (f(x) < f(sols[i])) x = sols[i];
	}
	//cout << counter << "\n";
	return x;
}
