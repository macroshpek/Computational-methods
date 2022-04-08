#include <iostream>
#include <cmath>

using namespace std;

const double EPS = 1e-8;
const double pi = 3.1415926536;

int wrachenie(double** coefficients, int numberOfEquation,
    double** solution, double precision) {
    int result = 1;
    int i, j, k;
    int maxI, maxJ;
    double max, fi;
    double** matricaPoworota;
    matricaPoworota = new double* [numberOfEquation];
    for (i = 0; i < numberOfEquation; i++) {
        matricaPoworota[i] = new double[numberOfEquation];
    }
    double** temp;
    temp = new double* [numberOfEquation];
    for (i = 0; i < numberOfEquation; i++) {
        temp[i] = new double[numberOfEquation];
    }
    double fault = 0.0;
    for (i = 0; i < numberOfEquation; i++) {
        for (j = i + 1; j < numberOfEquation; j++) {
            fault = fault + coefficients[i][j] * coefficients[i][j];
        }
    }
    fault = sqrt(2 * fault);
    while (fault > precision) {
        max = 0.0;
        for (i = 0; i < numberOfEquation; i++) {
            for (j = i + 1; j < numberOfEquation; j++) {
                if (coefficients[i][j] > 0 && coefficients[i][j] > max) {
                    max = coefficients[i][j];
                    maxI = i;
                    maxJ = j;
                }
                else if (coefficients[i][j] < 0 && -coefficients[i][j] > max) {
                    max = -coefficients[i][j];
                    maxI = i;
                    maxJ = j;
                }
            }
        }
        for (i = 0; i < numberOfEquation; i++) {
            for (j = 0; j < numberOfEquation; j++) {
                matricaPoworota[i][j] = 0;
            }
            matricaPoworota[i][i] = 1;
        }
        if (coefficients[maxI][maxI] == coefficients[maxJ][maxJ]) {
            matricaPoworota[maxI][maxI] = matricaPoworota[maxJ][maxJ] =
                matricaPoworota[maxJ][maxI] = sqrt(2.0) / 2.0;
            matricaPoworota[maxI][maxJ] = -sqrt(2.0) / 2.0;
        }
        else {
            fi = 0.5 * atan((2.0 * coefficients[maxI][maxJ]) /
                (coefficients[maxI][maxI] - coefficients[maxJ][maxJ]));
            matricaPoworota[maxI][maxI] = matricaPoworota[maxJ][maxJ] = cos(fi);
            matricaPoworota[maxI][maxJ] = -sin(fi);
            matricaPoworota[maxJ][maxI] = sin(fi);
        }
        for (i = 0; i < numberOfEquation; i++) {
            for (j = 0; j < numberOfEquation; j++) {
                temp[i][j] = 0.0;
            }
        }
        for (i = 0; i < numberOfEquation; i++) {
            for (j = 0; j < numberOfEquation; j++) {
                for (k = 0; k < numberOfEquation; k++) {
                    temp[i][j] = temp[i][j] + matricaPoworota[k][i] * coefficients[k][j];
                }
            }
        }
        for (i = 0; i < numberOfEquation; i++) {
            for (j = 0; j < numberOfEquation; j++) {
                coefficients[i][j] = 0.0;
            }
        }
        for (i = 0; i < numberOfEquation; i++) {
            for (j = 0; j < numberOfEquation; j++) {
                for (k = 0; k < numberOfEquation; k++) {
                    coefficients[i][j] = coefficients[i][j] +
                        temp[i][k] * matricaPoworota[k][j];
                }
            }
        }
        fault = 0.0;
        for (i = 0; i < numberOfEquation; i++) {
            for (j = i + 1; j < numberOfEquation; j++) {
                fault = fault + coefficients[i][j] * coefficients[i][j];
            }
        }
        fault = sqrt(2 * fault);
        for (i = 0; i < numberOfEquation; i++) {
            for (j = 0; j < numberOfEquation; j++) {
                temp[i][j] = 0.0;
            }
        }
        for (i = 0; i < numberOfEquation; i++) {
            for (j = 0; j < numberOfEquation; j++) {
                for (k = 0; k < numberOfEquation; k++) {
                    temp[i][j] = temp[i][j] + solution[i][k] * matricaPoworota[k][j];
                }
            }
        }
        for (i = 0; i < numberOfEquation; i++) {
            for (j = 0; j < numberOfEquation; j++) {
                solution[i][j] = temp[i][j];
            }
        }
        result++;
    }
    return result;
}

void eigenVectors(double** coefficients, int size, double*& tempp) {
    setlocale(LC_ALL, "rus");
    int i, j;
    double** solution, precision;
    solution = new double* [size];
    for (i = 0; i < size; i++) {
        solution[i] = new double[size];
    }
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            solution[i][j] = 0;
        }
        solution[i][i] = 1;
    }
    precision = 0.001;
    double** temp = new double* [size];
    for (int i = 0; i < size; i++) {
        temp[i] = new double[size];
        for (int j = 0; j < size; j++) {
            temp[i][j] = coefficients[i][j];
        }
        std::cout << std::endl;
    }
    int steps = wrachenie(coefficients, size, solution, precision);
    /*
    for (i = 0; i < size; i++) {
        cout << "Собственный вектор lambda = " << coefficients[i][i] << ":\n";
        for (j = 0; j < size; j++) {
            cout << solution[j][i] << "\n";
        }
    }
    for (int i = 0; i < size; i++) {
        std::cout << "__________________" << std::endl;
        for (int k = 0; k < size; k++) {
            double sum = 0;
            for (int r = 0; r < size; r++) {
                sum += temp[k][r] * solution[r][i];
            }
            std::cout << sum << " == " << solution[k][i] * coefficients[i][i] << std::endl;
        }
    }
    */
    for (int i = 0; i < size; i++) {
        tempp[i] = coefficients[i][i];
    }
}

double max_above_d(Matrix& a) {
    double s = 0;
    for (int i = 0; i < a.rows_count(); ++i) {
        for (int j = i + 1; j < a.rows_count(); ++j) {
            s += pow(a[j][i], 2);
        }
    }
    return sqrt(s);
}

int main() {
    // Размерность матрицы
    int n;
    // Считываем размерность матрицы
    cout << "dimenti0n of matrix: " << endl;
    cin >> n;
    // Матрица системы и Q матрица в RQ разложении
    Matrix A(n);
    double** temp_m = new double* [n];
    for (int i = 0; i < n; i++) {
        temp_m[i] = new double[n];
    }
    double* temp = new double[n];

    cout.setf(ios::fixed);
    cout.precision(8);

    // Считываем A
    cout << "Enter matrix A: " << endl;
    cin >> A;
    cout << endl;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            temp_m[i][j] = A[i][j];
        }
    }

    int i;

    for (i = 0; max_above_d(A) > EPS && i < 20; i++) {
        //        A = qr(A);
        auto qrs = A.qr();
        A = qrs.second * qrs.first;
    }

    cout << endl << n << endl;
    eigenVectors(temp_m, n, temp);
    for (int i = 0; i < n; i++) {
        A[i][i] = temp[i];
    }
    cout << A;

    // (QR)x = b сводится к задаче Rx = ~Qx = b'
    cout << "Iteration count : " << i << endl;

    return 0;
}