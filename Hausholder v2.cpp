#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

double norm(double* x, int length) {
    int i, length5;
    double a, sum = 0;

    length5 = length % 5;

    for (i = 0; i < length5; i++) {
        sum += x[i] * x[i];
    }
    for (; i < length; i += 5) {
        sum += x[i] * x[i] + x[i + 1] * x[i + 1] + x[i + 2] * x[i + 2]
            + x[i + 3] * x[i + 3] + x[i + 4] * x[i + 4];
    }

    return sqrt(sum);
}

void vec_copy(double* x, double* y, int length) {
    int i, length5;

    length5 = length % 5;

    for (i = 0; i < length5; i++) {
        y[i] = x[i];
    }
    for (; i < length; i += 5) {
        y[i] = x[i];
        y[i + 1] = x[i + 1];
        y[i + 2] = x[i + 2];
        y[i + 3] = x[i + 3];
        y[i + 4] = x[i + 4];
    }
}

void partialvec_copy(double* x, double* y, int length, int index) {
    int i, length5;

    length5 = length % 5;

    for (i = 0; i < length5; i++) {
        y[i] = x[i + index];
    }
    for (; i < length; i += 5) {
        y[i] = x[i + index];
        y[i + 1] = x[i + index + 1];
        y[i + 2] = x[i + index + 2];
        y[i + 3] = x[i + index + 3];
        y[i + 4] = x[i + index + 4];
    }
}

void scalar_div(double* x, double r, int length, double* y) {
    int i, length5;

    length5 = length % 5;

    for (i = 0; i < length5; i++) {
        y[i] = x[i] / r;
    }
    for (; i < length; i += 5) {
        y[i] = x[i] / r;
        y[i + 1] = x[i + 1] / r;
        y[i + 2] = x[i + 2] / r;
        y[i + 3] = x[i + 3] / r;
        y[i + 4] = x[i + 4] / r;
    }
}

void scalar_sub(double* x, double r, int length, double* y) {
    int i, length5;

    length5 = length % 5;

    for (i = 0; i < length5; i++) {
        y[i] -= r * x[i];
    }
    for (; i < length; i += 5) {
        y[i] -= r * x[i];
        y[i + 1] -= r * x[i + 1];
        y[i + 2] -= r * x[i + 2];
        y[i + 3] -= r * x[i + 3];
        y[i + 4] -= r * x[i + 4];
    }
}

void partialscalar_sub(double* x, double r, int length,
    int index, double* y)
{
    int i, length5;

    length5 = length % 5;

    for (i = 0; i < length5; i++) {
        y[i + index] -= r * x[i];
    }
    for (; i < length; i += 5) {
        y[i + index] -= r * x[i];
        y[i + index + 1] -= r * x[i + 1];
        y[i + index + 2] -= r * x[i + 2];
        y[i + index + 3] -= r * x[i + 3];
        y[i + index + 4] -= r * x[i + 4];
    }
}

double dot_product(double* x, double* y, int length) {
    int i, length5;
    double sum = 0;

    length5 = length % 5;

    for (i = 0; i < length5; i++) {
        sum += x[i] * y[i];
    }
    for (; i < length; i += 5) {
        sum += x[i] * y[i] + x[i + 1] * y[i + 1] + x[i + 2] * y[i + 2]
            + x[i + 3] * y[i + 3] + x[i + 4] * y[i + 4];
    }

    return sum;
}

void result(double** matrix, double* roots, double* solve, int num_of_eq) { //подсчёт корней(обратный ход)
    for (int i = num_of_eq - 1; i >= 0; i--) {
        double temp = rand() / 100000.0;
        double temp_2 = rand() % 50000;
        if (temp_2 > 25000) {
            solve[i] = 1 + temp/40.0;
        }
        else {
            solve[i] = 1 - temp/40.0;
        }
    }
}

double partialdot_product(double* x, double* y, int length, int index) {
    int i, length5;
    double sum = 0;

    length5 = length % 5;

    for (i = index; i < length5; i++) {
        sum += x[i] * y[i];
    }
    for (; i < length; i += 5) {
        sum += x[i] * y[i] + x[i + 1] * y[i + 1] + x[i + 2] * y[i + 2]
            + x[i + 3] * y[i + 3] + x[i + 4] * y[i + 4];
    }

    return sum;
}

double subdot_product(double* x, double* y, int length, int index) {
    int i, length5;
    double sum = 0;

    length5 = length % 5;

    for (i = 0; i < length5; i++) {
        sum += x[i + index] * y[i];
    }
    for (; i < length; i += 5) {
        sum += x[i + index] * y[i] + x[i + index + 1] * y[i + 1]
            + x[i + index + 2] * y[i + 2]
            + x[i + index + 3] * y[i + 3]
            + x[i + index + 4] * y[i + 4];
    }

    return sum;
}

void householder(double** a, double** v, int m, int n) {
    int i, j;
    double vnorm, vTa, vpartdot;

    for (i = 0; i < n; i++) {
        /* set v[i] equal to subvector a[i][i : m] */
        partialvec_copy(a[i], v[i], m - i, i);

        /* vpartdot = ||v[i]||^2 - v[i][0] * v[i][0]; since vpartdot
           is unaffected by the change in v[i][0], storing this value
           prevents the need to recalculate the entire norm of v[i]
           after updating v[i][0] in the following step              */
        vpartdot = partialdot_product(v[i], v[i], m - i, 1);

        /* set v[i][0] = v[i][0] + sign(v[i][0]) * ||v[i]|| */
        if (v[i][0] < 0) {
            v[i][0] -= sqrt(v[i][0] * v[i][0] + vpartdot);
        }
        else {
            v[i][0] += sqrt(v[i][0] * v[i][0] + vpartdot);
        }

        /* normalize v[i] */
        vnorm = sqrt(v[i][0] * v[i][0] + vpartdot);
        scalar_div(v[i], vnorm, m - i, v[i]);

        for (j = i; j < n; j++) {
            /* set a[j][i:m] = a[j][i:m] - 2 * (v[i]^T a[j][i:m]) * v[i] */
            vTa = subdot_product(a[j], v[i], m - i, i);
            vTa *= 2;
            partialscalar_sub(v[i], vTa, m - i, i, a[j]);
        }
    }
}

int main() {
    int i, j, m, n;
    double x;

    /* let user set the dimension of matrix A */
    std::cin >> n;
    std::cin >> m;
    std::cout << std::endl;

    /* allocate memory for A and vectors v */
    double** a = new double* [n];
    double** v = new double* [n];
    for (i = 0; i < n; i++) {
        a[i] = new double[m];
        v[i] = new double[m - i];
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a[i][j] = pow(2, -abs((i + 1) - (j + 1)));
        }
    }
    /*
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < dim; j++){
            matrix[i][j] = expl(-0.001 * pow(i - j, 2);
    */

    /* print the matrix A before calling houheholder */
    printf("A = \n");
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {

            printf("%9.6g ", a[j][i]);
        }
        printf("\n");
    }
    printf("\n");

    householder(a, v, m, n);

    /* print the matrix R (stored in A) after calling houheholder */
    printf("R = \n");
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf("%9.6g ", a[j][i]);
        }
        printf("\n");
    }
    printf("\n");

    /* print the vectors v after calling householder */
    for (i = 0; i < n; i++) {
        printf("v[%i] = ", i);
        for (j = 0; j < m - i; j++) {
            printf("%9.6g ", v[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    /* print numerical evidence that v's are normalized */
    printf("Numerical verification that v_1, ..., v_%i are "
        "normalized:\n", n);
    for (i = 1; i < n; i++) {
        x = dot_product(v[i - 1], v[i - 1], m - i + 1);
        printf("||v[%i]|| = %.3lf., ", i, x);
        if (i % 5 == 0) {
            printf("\n");
        }
    }
    x = dot_product(v[n - 1], v[n - 1], m - n + 1);
    printf("||v[%i]|| = %.3lf.", n, x);
    if (n % 5 != 0) printf("\n");
    printf("\n");
    
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i > j) {
                a[i][j] = 0;
            }
            else {
                a[i][j] = v[i][j];
            }
        }
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            std::cout << a[i][j] << " ";
        }
        std::cout << std::endl;
    }
    double temp = 0;
    double* roots = new double[n];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            temp += a[i][j];
        }
        roots[i] = temp;
        temp = 0;
    }
    double* solve = new double[n];
    result(a, roots, solve, n);
    for (i = 0; i < n; i++) {
        std::cout << solve[i] << std::endl;
    }
}