#define  _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <malloc.h>
#include < math.h >

void main() {
    int N, k; // n-размер матрицы 
    int i, j;

    scanf("%d", &N);

    double** a;
    a = (double**)malloc(N * sizeof(double*)); //    т.е. кол-во строк
    for (i = 0; i < N; i++) {
        a[i] = (double*)malloc(N * sizeof(double)); //    т.е. кол-во столбцов
    }
    double* b;
    b = (double*)malloc(N * sizeof(double)); //    b1-bn
    double* x;
    x = (double*)malloc(N * sizeof(double)); //    x1-xn
    for (i = 0; i < N; i++) {
        x[i] = 0;
    }

    double** r;
    r = (double**)malloc(N * sizeof(double*)); //    т.е. кол-во строк
    for (i = 0; i < N; i++) {
        r[i] = (double*)malloc(N * sizeof(double)); //    т.е. кол-во столбцов
    }
    double* br;
    br = (double*)malloc(N * sizeof(double)); //    b1-bn


    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            r[i][j] = 0;
        }
    }

    for (i = 0; i < N; i++) { // читаю элементы
        for (j = 0; j < N; j++) {
            scanf( "%lf ", &a[i][j]);//  кол-во строк кол-во столбцов
        }
       scanf("%lf ", &b[i]);//  читаю b
    }

    for (int ii = 0; ii < N; ii++) {
        for (int jj = 0; jj < (N); jj++) {
            printf_s("aa %4.3lf ", a[ii][jj]);
        }
        printf_s("%4.3lf ", b[ii]);
        printf_s("\n ");
    }printf_s("\n "); printf_s("\n ");

    double c = 0, s = 0;
    for (i = 0; i < N - 1; i++) {//удаляем х1-хM
        for (k = 1; k < N - i; k++) {
            c = a[i][i] / sqrt(a[i][i] * a[i][i] + a[i + k][i] * a[i + k][i]);
            s = a[i + k][i] / sqrt(a[i][i] * a[i][i] + a[i + k][i] * a[i + k][i]);

            for (j = 0; j < N; j++) {
                r[i][j] = c * a[i][j] + s * a[i + k][j];
                r[i + k][j] = c * a[i + k][j] - s * a[i][j];
            }
            br[i] = c * b[i] + s * b[i + k];
            br[i + k] = c * b[i + k] - s * b[i];

            for (j = 0; j < N; j++) {
                a[i][j] = r[i][j];
                a[i + k][j] = r[i + k][j];

            }
            b[i] = br[i];
            b[i + k] = br[i + k];

            for (int ii = 0; ii < N; ii++) {
                for (int jj = 0; jj < (N); jj++) {
                    printf_s("aa %lf ", a[ii][jj]);
                }
                printf_s("aa %lf ", b[ii]); printf_s("\n ");
            }printf_s("\n "); printf_s("\n ");

        }

    }



    // вывод   матрицы


    for (int ii = 0; ii < N; ii++) {
        for (int jj = 0; jj < (N); jj++) {
            printf_s("aa %lf ", a[ii][jj]);
        }
        printf_s("aa %lf ", b[ii]); printf_s("\n ");
    }printf_s("\n "); printf_s("\n ");


    //cчитаю иксы
    double t = 0;
    for (i = N - 1; i >= 0; i--) {
        t = 0;

        for (int j = i; j < N; j++) {
            t = t + a[i][j] * x[j];
        }
        x[i] = (b[i] - t) / a[i][i];

    }



    //вывод иксов
    printf_s("\n "); printf_s("\n ");

    for (i = 0; i < N; i++) {
        printf("%lf \n", x[i]);
    }


}