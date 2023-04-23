
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

double *ans;
double tmp_num;

void print_matrix(double **arr, int n, double *f);
void print_vector(double *x, int n);
double** copy_matrix(double **a, int n);
double* copy_vector(double *x, int n);
void gen_matrix(int n, double ***arr);
void find_f_vector(double **arr, double **f, int n);

double convergence(double *ans, double **matrix, double *f, int n) {
    double norm = 0., sum;
    int i, j;
    for (i = 0; i < n; i++) {
        sum = 0.;
        for (j = 0; j < n; j++) {
            sum += matrix[i][j] * ans[j];
        }
        sum -= f[i];
        norm += sum * sum;
    }
    return sqrt(norm);
}

int over_relaxation_method(double **matrix, double *f, int n, double omega, double eps) {
    int i, j, step = 0;
    double *pre_ans, sum;
    pre_ans = calloc(n, sizeof(*pre_ans));
    memset(pre_ans, 0, n);
    for (i = 0; i < n; i++) ans[i] = 0;
    while ((!step) || (convergence(ans, matrix, f, n) > eps)) {
        step++;
        for (i = 0; i < n; i++) {
            pre_ans[i] = ans[i];
        }
        for (i = 0; i < n; i++) {
            sum = 0.;
            for (j = 0; j < i; j ++){
                sum += matrix[i][j] * ans[j];
            }
            for (j = i; j < n; j++) {
                sum += matrix[i][j] * pre_ans[j];
            }
            ans[i] = pre_ans[i] + (f[i] - sum) * (omega / matrix[i][i]);
        }
    }
    return step;
}

double **matrix_transpose(double **matrix, int n) {
    int i, j;
    double **res;
    res = calloc(n, sizeof(*res));
    for (i = 0; i < n; i++) {
        res[i] = calloc (n, sizeof(**res));
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            res[i][j] = matrix[j][i];
        }
    }
    return res;
}


double **matrix_multiplication(double **matrix_1, double **matrix_2, int n) {
    int i, j, k;
    double **res, sum;
    res = calloc(n, sizeof(*res));
    for (i = 0; i < n; i++) {
        res[i] = calloc (n, sizeof(**res));
    }
    for (i = 0; i < n; i++) {
        for (k = 0; k < n; k++) {
            sum = 0.;
            for(j = 0; j < n; j++) {
                sum += matrix_1[i][j] * matrix_2[j][k];
            }
            res[i][k] = sum;
        }
    }
    return res;
}

double *matr_to_vec_mult(double **matr, double *vector, int n) {
    double *res, sum;
    int i, j;
    res = calloc(n, sizeof(*res));
    for (i = 0; i < n; i++) {
        sum = 0.;
        for (j = 0; j < n; j++) {
            sum += matr[i][j] * vector[j];
        }
        res[i] = sum;
    }
    return res;
}

void free_all(double **m, double *v, int n);

int main(int argc, char **argv) {
    double tmp_num = 0.1, w, eps = 10e-6;
    printf("Введите способ ввода матриц\n");
    printf("0 - по формуле, или 1 или 2 или 3 - номер примера\n");
    int n = 4, flag;
    scanf("%d", &flag);
    if (flag == 0) n = 20;

    double **matrix = calloc(n, sizeof(*matrix));
    for (int i = 0; i < n; i++) matrix[i] = calloc(n, sizeof(**matrix));
    double *f_vector = calloc(n, sizeof(*f_vector));
    ans = calloc(n, sizeof(*ans));
    if (flag == 0) {
        gen_matrix(n, &matrix);
        find_f_vector(matrix, &f_vector, n);
        printf("\nСгенерированная система:\n");
        printf("Сгенерированный вектор x всегда единичный вектор длины n\n");
        print_matrix(matrix, n, f_vector);
        
    } else {
        FILE *fd;
        if (flag == 1) {
            fd = fopen("/home/uliana/Programm/CHMI/1.txt", "r");
        } else if (flag == 2) {
            fd = fopen("/home/uliana/Programm/CHMI/2.txt", "r");
        } else {
            fd = fopen("/home/uliana/Programm/CHMI/3.txt", "r");
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                fscanf(fd, "%lf", &matrix[i][j]);
            }
            fscanf(fd, "%lf", &f_vector[i]);
        }
        fclose(fd);
    }
    
    
    double **matrix_1 = copy_matrix(matrix, n);
    double *vector_f1 = copy_vector(f_vector, n);

    printf("Приведём матрицу и вектор правой части к симметричному виду:\n");
    matrix_1 = matrix_transpose(matrix_1, n);
    matrix = matrix_multiplication(matrix_1, matrix, n);
    vector_f1 = matr_to_vec_mult(matrix_1, vector_f1, n);
    print_matrix(matrix, n, vector_f1);
    printf ("Заданная точность: %lf\n", eps);
    
    for (w = 0.1; w < 2; w += 0.1) {
        int it = over_relaxation_method(matrix,vector_f1, n, w, eps);
        printf ("w = %.2lf количество итераций - %d\nОтвет:", w,it);
        print_vector(ans, n);
    }

    free_all(matrix_1, vector_f1, n);
    free_all(matrix, f_vector, n);

    return 0;
}
