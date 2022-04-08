#include <iostream>
#include <fstream>

using namespace std;

double fxy(double x, double y) {
	//dy/dx = y - x
	double fxy = y - x;
	return fxy;
}

void runge() {
	double x0 = 0;
	double y0 = 1.5;
	double x = x0;
	double y = y0;
	double h = 0.25;
	double k1;
	double k2;
	double k3;
	double k4;

	for (double x = x0; x < 1.5; x += h) {
		cout << x << "\t" << y << "\n";
		k1 = h * fxy(x, y); 
		k2 = h * fxy(x + 0.5*h, y + 0.5*k1);
		k3 = h * fxy(x + 0.5 * h, y + 0.5 * k2);
		k4 = h * fxy(x +   h, y + k3);

		y = y + 1.0 / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
	}
}

void calculate(double* y) {
	double k1, k2, k3, k4, x, yi, h;
	double x0 = 0;
	double y0 = 1.5;
	h = 0.25;  // step
	yi = y0; x = x0;
	for (int i = 0; i < 4; i++) { //just runge kooooooooooooota
		k1 = h * fxy(x, yi);
		k2 = h * fxy(x + h / 2, yi + k1 / 2);
		k3 = h * fxy(x + h / 2, yi + k2 / 2);
		k4 = h * fxy(x + h, yi + k3);
		yi += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
		x += h;
		*(y + i + 1) = yi;
	}
}

void adams() {
	const int num_points = 10;
	double y[num_points + 1], h;
	double x0 = 0;
	double x = x0 + 4*0.25;
	double y0 = 1.5;
	y[0] = y0;
	// apply Runge-Kutta method
	calculate(y);
	h = 0.25;
	// extrapolating
	for (int i = 4; i < num_points; i++) {  //тут начинается с 4х ибо рунга уже всё вычислил, молодчина( начальные значения берутся из одношагового метода рунге, рунге вообще хороший мужик был)
		y[i] = y[i - 1] + h / 24 * (55 * fxy(x0 + (i - 1) * h, y[i - 1]) -
			59 * fxy(x0 + (i - 2) * h, y[i - 2]) +
			37 * fxy(x0 + (i - 3) * h, y[i - 3]) -												//это по формулам из википедии, полюбому рунга их написал, по-другому он не мог поступить
			9 * fxy(x0 + (i - 4) * h, y[i - 4]));
		x += h;
	}
	for (int i = 0; i < num_points; i++)
		printf("%f\t%f\t%f\n", (x0 + i * h), y[i], (2 * exp(x0 + i * h) - (x0 + i * h) - 1));
}

void euler() {
	double x0 = 0;
	double y0 = 1.5;
	double x = x0;
	double y = y0;
	double h = 0.25;
	for (double x = x0; x < 1.5; x += h) {
		cout << x << "\t" << y << "\n";

		y = y + h * fxy(x + h / 2, y + h / 2 * fxy(x, y));
	}

}

int main() {
	cout << "x" << "\t" << "y" << "\n";
	cout << "roongeee-cootaaa" << "\n";
	runge();
	cout << "eeeeeuuuuleeer" << "\n";
	euler();
	cout << "aaaaadaamts" << "\n";
	adams();
	
}

//самый точный метод - рунге кута
// эйлер накапливает много ошибок, и вообще рунге модифицировал его аппроксимирую остаточные члены
//у рунге коэфф остаточного члена меньше в 960 раз чем у адамса, и позволяет брать шаг в 4 раза больше чем у адамса