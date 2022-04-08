#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#define _CRT_SECURE_NO_WARNINGS 
#pragma warning(disable : 4996)

void diagonal(double** matrix, int k, int num_of_eq, int* changes) {                      //нахождение макс элемента в столбце и проверка его на ноль, если макс = 0, то нахождение минимального и умножение всей строки на -1
	double temp = 0;
	int i_max = k;
	double max = matrix[k][k];

	for (int i = k; i < num_of_eq; i++) {
		if (max < matrix[i][k]) {
			max = matrix[i][k];
			i_max = i;
		}
	}
	if (k != i_max) {
		(*changes)++;
	}

	int min = 0;
	if (max == 0) {
		for (int i = k; i < num_of_eq; i++) {
			if (min > matrix[i][k]) {
				min = matrix[i][k];
				i_max = i;
			}
		}
		for (int i = k; i < num_of_eq; i++) {
			matrix[i_max][i] = (-1) * matrix[i_max][i];
		}
	}

	int j;

	for (j = 0; j <= num_of_eq; j++) {
		temp = matrix[k][j];
		matrix[k][j] = matrix[i_max][j];
		matrix[i_max][j] = temp;
	}
}

void print_matrix(double** matrix, int num_of_eq, int mode) {									//печать матрицы
	if (mode == 0) {
		for (int i = 0; i < num_of_eq; i++) {
			for (int j = 0; j <= num_of_eq; j++) {
				printf("%.1lf ", matrix[i][j]);
			}
			printf("\n");
		}
	}
	else {
		for (int i = 0; i < num_of_eq; i++) {
			for (int j = 0; j < num_of_eq; j++) {
				printf("%.1lf ", matrix[i][j]);
			}
			printf("\n");
		}
	}
}

double determinant(double** matrix, int num_of_eq, int* changes);

void triangle_matrix(double** matrix, int num_of_eq) {                      //привод к триугольной   
	
	int* changes = 0;
	for (int k = 0; k < num_of_eq; k++) {
		for (int i = k + 1; i < num_of_eq; i++) {			
			diagonal(matrix, k, num_of_eq, changes);
			if (matrix[k][k] == 0) {  //если макс = 0 - решение бесконечное множество
				printf("Solution not exist\n");
				exit(0);
			}
			double M = matrix[i][k] / matrix[k][k];
			for (int j = k; j <= num_of_eq; j++) {
				matrix[i][j] -= M * matrix[k][j];
			}
		}
	}
	double det = 0;
	printf("determinant\n");
	det = determinant(matrix, num_of_eq, changes);
	printf("%lf\n", det);
}

void result(double** matrix, double* roots, int num_of_eq) { //подсчёт корней(обратный ход)
	for (int i = num_of_eq - 1; i >= 0; i--) {
		double s = 0; 
		for (int j = i; j < num_of_eq; j++) {
			s = s + matrix[i][j] * matrix[i][num_of_eq];
		}
		roots[i] = (matrix[i][num_of_eq] - s) / matrix[i][i];
	}
}

double determinant(double** matrix, int num_of_eq, int* changes) { //определитель(умножение всех диаг элем-ов

	double det = 1;
	for (int i = 0; i < num_of_eq; i++) {
		det *= matrix[i][i];
	}
	det = det * pow(-1, *changes);
	return det;
}

void matrix_copy(double** matrix, double** copy, int num_of_eq){ //копирование матрицы

	for (int i = 0; i < num_of_eq; i++) {
		for (int j = 0; j < num_of_eq; j++) {
			copy[i][j] = matrix[i][j];
		}
	}
}

void multiple(double** matrix, double** r_matrix, int num_of_eq) { //умножение матрицы

	double** mult = (double**)malloc(num_of_eq * sizeof(double*));
	if (!matrix) {
		exit(0);
	}

	for (int i = 0; i < num_of_eq; i++) {
		mult[i] = (double*)malloc((num_of_eq + 1) * sizeof(double));
		if (!matrix[i]) {
			exit(0);
		}
	}
	int str = 0;
	for (int i = 0; i < num_of_eq; i++) {
		for (int j = 0; j < num_of_eq; j++) {
			mult[i][j] = 0;
			for (int k = 0; k < num_of_eq; k++) {
				mult[i][j] += matrix[i][k] * r_matrix[k][j];
			}
		}
	}
	printf("\nmultpile of matrix and reverse\n"); //вывод умножения данной на обратную к ней(проверк на обратность)
	print_matrix(mult, num_of_eq, 1);
}

void reverse_matrix(double** matrix, double** r_matrix, int num_of_eq) { //нахождение обратной через единичную(когда рядом единичиную приписываешь и приводишь к диаг виду)

	double** E = (double**)malloc(num_of_eq * sizeof(double*));
	if (!matrix) {
		exit(0);
	}

	for (int i = 0; i < num_of_eq; i++) {
		E[i] = (double*)calloc(num_of_eq, sizeof(double));
		if (!E[i]) {
			exit(0);
		}
	}

	double** copy = (double**)malloc(num_of_eq * sizeof(double*));
	if (!matrix) {
		exit(0);
	}

	for (int i = 0; i < num_of_eq; i++) {
		copy[i] = (double*)calloc(num_of_eq, sizeof(double));
		if (!copy[i]) {
			exit(0);
		}
	}

	matrix_copy(matrix, copy, num_of_eq);

	for (int i = 0; i < num_of_eq; i++) {
		E[i][i] = 1;
	}

	double temp;

	for (int k = 0; k < num_of_eq; k++) {
		temp = copy[k][k];
		for (int j = 0; j < num_of_eq; j++) {
			copy[k][j] /= temp;
			E[k][j] /= temp;
		}
		for (int i = k + 1; i < num_of_eq; i++) {
			temp = copy[i][k];
			for (int j = 0; j < num_of_eq; j++) {
				copy[i][j] -= copy[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}
	for (int k = num_of_eq - 1; k > 0; k--)	{
		for (int i = k - 1; i >= 0; i--){
			temp = copy[i][k];
			for (int j = 0; j < num_of_eq; j++) {
				copy[i][j] -= copy[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	matrix_copy(E, r_matrix, num_of_eq);

	print_matrix(r_matrix, num_of_eq, 1);

	multiple(matrix, r_matrix, num_of_eq);
}

void output_res(double* roots, int num_of_eq) { //вывод в файл

	FILE* out = fopen("out.txt", "w");
	for (int i = 0; i < num_of_eq; i++) {
		fprintf(out, "%.3lf\n", roots[i]);
	}
}


void main() {

	int num_of_eq = 0;
	scanf("%d", &num_of_eq);
																					//выделение памяти под матрицу
	double** matrix = (double**)malloc(num_of_eq * sizeof(double*));
	if (!matrix) {
		exit(0);
	}

	for (int i = 0; i < num_of_eq; i++) {
		matrix[i] = (double*)calloc((num_of_eq + 1), sizeof(double));
		if (!matrix[i]) {
			exit(0);
		}
	}
																					//создание памяти для обратной
	double** r_matrix = (double**)malloc(num_of_eq * sizeof(double*));
	if (!matrix) {
		exit(0);
	}

	for (int i = 0; i < num_of_eq; i++) {
		r_matrix[i] = (double*)malloc(num_of_eq * sizeof(double));
		if (!matrix[i]) {
			exit(0);
		}
	}
																					//память для корней
	double* roots = (double*)calloc(num_of_eq , sizeof(double));
	if (!roots) {
		exit(0);
	}

	 /*double temp_i, temp_j;
	for (int i = 0; i < num_of_eq; i++) {
		for (int j = 0; j < num_of_eq; j++) {
			if (j >= i) {
				temp_i = (double)(i + 1);
				temp_j = (double)(j + 1);
				matrix[i][j] = (temp_i * 1.0 )/( 1.0 * temp_j);
			}
			else {
				temp_i = (double)(i + 1);
				temp_j = (double)(j + 1);
				matrix[i][j] = (temp_j * 1.0) / (1.0 * temp_i);
			}
			temp_i = 0;
			temp_j = 0;
		}
	}
	for (int i = 0; i < num_of_eq; i++) {
		matrix[i][num_of_eq] = 1;
	}
	*/

	int temp;
	temp = num_of_eq - 1;
	for(int i = 0; i < num_of_eq; i++) {
		matrix[i][temp] = 1;
		temp--;
	}
	for(int i = 0; i < num_of_eq; i++){
		matrix[i][num_of_eq] = 1;
	}
	
																					//вывод заданной матрицы
	printf("first matrix\n");
	print_matrix(matrix, num_of_eq, 0);
																					//вывод обратной
	printf("\nreverse matrix\n");
	reverse_matrix(matrix, r_matrix, num_of_eq);
																					//вывод треугольной
	printf("\ntriangle matrix\n");
	triangle_matrix(matrix, num_of_eq);
	print_matrix(matrix, num_of_eq, 0);
	printf("\n");
																					//вывод корней в файл
	result(matrix, roots, num_of_eq);
	output_res(roots, num_of_eq);
																					//определитель																	

	for (int i = 0; i < num_of_eq; i++) {                                            //освобождение памяти
		free(matrix[i]);
		free(r_matrix[i]);
	}
	free(matrix);
	free(r_matrix);
	free(roots);
}