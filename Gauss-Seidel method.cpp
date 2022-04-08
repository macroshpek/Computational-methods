#include <iostream>
#include <cmath>
#include <ctime>
#include <iomanip> 
using namespace std;

void SLAU_output(double** A, double* b, int m, int n) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			cout << A[i][j] << " * x" << j;
			if (j < n - 1)
				cout << " + ";
		}
		cout << "  =  " << b[i] << endl;
	}
}

double get_norm(double* x1, double* x0, int n)
{
	double norm = 0;

	for (int i = 0; i < n; i++)
	{
		norm += pow(x1[i] - x0[i], 2);
	}

	return (sqrt(norm)); 
}

double rounding(double d, double digits)
{
	d = round(d * pow(10, digits)) / pow(10, digits);
	return d;
}

bool diagonal(double** a, int n)
{
	bool domination = true; 

	for (int i = 0; i < n; i++)
	{
		double sum = 0;
		for (int j = 0; j < n; j++) sum += abs(a[i][j]);
		sum -= abs(a[i][i]);

		if (abs(a[i][i]) > sum)
		{
			cout << a[i][i] << " > " << sum << " - OK." << endl;
		}
		else
		{
			domination = false;
			cout << a[i][i] << " < " << sum << " - NO." << endl;
		}
	}

	cout << endl;
	return (domination);
}

int main()
{
	setlocale(LC_ALL, "");
	unsigned int digits;
	int n;
	double** A;
	double* b;
	cout << "размерность:";
	cin >> n;
	cout << endl;
	cout << "точность:";
	cin >> digits;
	cout << endl;
	A = new double* [n];
		for (int i = 0; i < n; i++)
		{
			A[i] = new double[n];
		}
		b = new double[n];
		std::cout << endl << "Введите коэффициенты и свободные члены:" << endl;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				cin >> A[i][j];
			}
			cin >> b[i];
		}

	const char* out_line = "\n__________________________________________________________________________\n";
	std::cout << out_line << endl;
	if (diagonal(A, n)) {
		double* x1, * x0;
		x1 = new double[n];
		for (int i = 0; i < n; i++) {
			x1[i] = 0;
		}
		x0 = new double[n];
		int k = 0;
		double norm;
		do {
			for (int i = 0; i < n; i++)
				x0[i] = x1[i];

			for (int i = 0; i < n; i++) {
				double f = 0;
				for (int j = 0; j < i; j++)
					f += (A[i][j] * x1[j]);
				for (int j = i + 1; j < n; j++)
					f += (A[i][j] * x0[j]);
				x1[i] = (b[i] - f) / A[i][i];
			}

			k++;
			norm = get_norm(x1, x0, n);
		}
		while (!(norm < 1 / pow(10, digits)));
		cout << "Решение системы:" << endl << endl;
		for (int i = 0; i < n; i++)
		{
			cout << "x[" << i << "] = " << setprecision(digits + 1) << x1[i] << endl;
		}
		cout << endl << "Количество итераций: " << k;
	}

	cout << endl << out_line << endl << endl;
	return 1;
}