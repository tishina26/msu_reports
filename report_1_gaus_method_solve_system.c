#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

static double det = 1.0;
double tmp_num;

void print_matrix(double **arr, int n, double *f) {
    for (int i = 0; i < n; i++) {
        printf("|");
        for (int j = 0; j < n; j++) {
            printf("%10.5lf ", arr[i][j]);
        }

        printf(" |  *  |x|  =  ");
        if (f) {
            printf("|%10.5lf |\n", f[i]);
        } else {
            printf("| f |\n");
        }
    }

    printf("\n");

    return;
}

void print_vector(double *x, int n) {
    for (int i = 0; i < n; i++) {
        printf("%10.5lf ", x[i]);
    }

    printf("\n");
}

double** copy_matrix(double **a, int n) {
    double** arr = calloc(n, sizeof(*arr));

    for (int i = 0; i < n; i++) {
        arr[i] = calloc(n, sizeof(*arr[i]));

        for (int j = 0; j < n; j++) {
            arr[i][j] = a[i][j];
        }
    }

    return arr;
}

double* copy_vector(double *x, int n) {
    double* ans = calloc(n, sizeof(*ans));

    for (int i = 0; i < n; i++) {
        ans[i] = x[i];
    }

    return ans;
}

void gen_matrix(int n, double ***arr) {
    int num = 2, mod = 7;
    double c = 0.01;
    tmp_num = 0.0;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                ((*arr)[i][j]) = (double) (((i + 1) * num + j) % mod);
            } else {
                ((*arr)[i][j]) = c;
            }
            if ((*arr)[i][j] > tmp_num) {
                tmp_num = (*arr)[i][j];
            }
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            (*arr)[i][j] /= tmp_num;
        }
    }
}

void find_not_null_elem_1(double** arr, int n, double* f, int i) {
    int k = i;
    while ((k < n - 1) && (arr[k][i] == 0.0)) {
        k++;
    }

    double *tmp = arr[i];
    arr[i] = arr[k];
    arr[k] = tmp;

    if (f) {
        double tmp1 = f[i];
        f[i] = f[k];
        f[k] = tmp1;
    }
    return;
}

void find_not_null_elem_max(double** arr, int n, double* f, int i) {
    int k = i, i_max = i;
    double max = arr[i][i];

    for (int k = i; k < n - 1; k++) {
        if ((arr[k][i] != 0) && (arr[k][i] > max)) {
            max = arr[k][i];
            i_max = k;
        }
    }

    double *tmp = arr[i];
    arr[i] = arr[i_max];
    arr[i_max] = tmp;

    if (f) {
        double tmp1 = f[i];
        f[i] = f[i_max];
        f[i_max] = tmp1;
    }
    return;
}

void make_triangle(double** arr, int n, double* f,
    void find_not_null_elem(double** arr, int n, double* f, int i)) {

    for (int i = 0; i < n - 1; i++) {
        (*find_not_null_elem)(arr, n, f, i);
        double koef = arr[i][i];
        if (koef == 0) {
            printf("Определитель матрицы = 0, пожалуйста, введите другую матрицу\n");
            exit(0);
        }

        if (f) {
            f[i] /= koef;
        }
        for (int j = i; j < n; j++) {
            arr[i][j] /= koef;
        }
        for (int j = i + 1; j < n; j++) {
            koef = arr[j][i];
            if (f) {
                f[j] -= koef * f[i];
            }
            for (int k = i; k < n; k++) {
                arr[j][k] -= koef * arr[i][k];
            }
        }
    }
}

void find_f_vector(double **arr, double **f, int n) {
    for (int i = 0; i < n; i++) {
        double a = 0.0;
        for (int j = 0; j < n; j++) {
            a += arr[i][j];
        }
        (*f)[i] = a;
    }
}

double* GAUSS_method(double **arr, double *f, int n, int flag) {
    double *x = calloc(n, sizeof(*x));
    
    printf("\nПриведенная система:\n");

    if (flag = 1) {
        make_triangle(arr, n, f, find_not_null_elem_1);
    } else {
        make_triangle(arr, n, f, find_not_null_elem_max);
    }

    print_matrix(arr, n, f);

    x[n - 1] = f[n - 1] / arr[n - 1][n - 1];

    for (int i = n - 2; i >= 0; i--) {
        arr[i][i] = f[i];

        for (int j = i + 1; j < n; j++) {
            arr[i][i] -= arr[i][j] * x[j]; 
        }

        x[i] = arr[i][i];
    }

    return x;
}

double determinant(double **matr, int n){
    int i, j, k;
    double max_elem, buf, det = 1, buf_matr;
    for (i = 0; i < n; i++) {
        max_elem = matr[i][i];
        k = i;
        for(j = i + 1; j < n; j++) {
            if (fabs(max_elem) < fabs(matr[j][i])) {
                max_elem = matr[j][i];
                k = j;
            }
        }
        if (k != i) {
            det = det * (-1);
            for (j = i; j < n; j++) {
                buf = matr[i][j];
                matr[i][j] = matr[k][j];
                matr[k][j] = buf;
            }
        }
        det *= matr[i][i];
        buf_matr = matr[i][i];
        for (j = i; j < n; j++) {
            matr[i][j] /= buf_matr;
        }
        for (k = i + 1; k < n; k++) {
            if (matr[k][i] != 0) {
                for (j = i + 1; j < n; j++) {
                    matr[k][j] -= matr[i][j] * matr[k][i];
                }
                matr[k][i] = 0;
            }
        }
    }
    return det;
}

double **inverse_matrix(double **matr, int n){
    double **res, max_elem, buf, buf_matr;
    int i, j, k;
    res = calloc(n, sizeof(*res));
    for (i = 0; i < n; i++) {
        res[i] = calloc (n, sizeof(**res));
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) {
                res[i][j] = 1.0;
            } else {
                res[i][j] = 0.0;
            }
        }
    }
    for (i = 0; i < n; i++) {
        max_elem = matr[i][i];
        k = i;
        for(j = i + 1; j < n; j++) {
            if (fabs(max_elem) < fabs(matr[j][i])) {
                max_elem = matr[j][i];
                k = j;
            }
        }
        if (k != i) {
            for (j = i; j < n; j++) {
                buf = matr[i][j];
                matr[i][j] = matr[k][j];
                matr[k][j] = buf;
            }
            for (j = 0; j < n; j++) {
                buf = res[i][j];
                res[i][j] = res[k][j];
                res[k][j] = buf;
            }
        }
        buf_matr = matr[i][i];
        for (j = i; j < n; j++) {
            matr[i][j] /= buf_matr;
        }
        for (j = 0; j < n; j++) {
            res[i][j] /= buf_matr;
        }
        for (k = i + 1; k < n; k++) {
            if (matr[k][i] != 0) {
                for (j = 0; j < n; j++) {
                    res[k][j] -= res[i][j] * matr[k][i];
                }
                for (j = i + 1; j < n; j++) {
                    matr[k][j] -= matr[i][j] * matr[k][i];
                }
                matr[k][i] = 0;
            }
        }
    }
    for (i = n - 1; i >= 0; i--) {
        for (j = i - 1; j >= 0; j--) {
            for (k = 0; k < n; k++) {
                res[j][k] -= res[i][k] * matr[j][i];
            }
        }
    }
    return res;
}

void free_all(double **m, double *v, int n) {
    for (int i = 0; i < n; i++) {
        free(m[i]);
    }
    free(m);
    free(v);
}

int main(int argc, char **argv) {
    printf("Введите способ ввода матриц\n");
    printf("0 - по формуле, или 1 или 2 или 3 - номер примера\n");
    int n = 4, flag;
    scanf("%d", &flag);
    if (flag == 0) n = 40;

    double **matrix = calloc(n, sizeof(*matrix));
    for (int i = 0; i < n; i++) matrix[i] = calloc(n, sizeof(**matrix));
    double *f_vector = calloc(n, sizeof(*f_vector));
    
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

    double **matrix_1 = nULL;
    double *vector_f1 = nULL;

    matrix_1 = copy_matrix(matrix, n);
    vector_f1 = copy_vector(f_vector, n);

    
    if (flag == 0) {
        printf("\tОпределитель равен %lf * %lf ^ %d\n\n", determinant(matrix_1, n), tmp_num, n);
    } else printf("\tОпределитель равен %lf\n\n", determinant(matrix_1, n));
    free_all(matrix_1, vector_f1, n);
    
    
    matrix_1 = copy_matrix(matrix, n);
    vector_f1 = copy_vector(f_vector, n);
    
    
    matrix_1 = inverse_matrix(matrix_1, n);
    printf("\tОбратная матрица:\n");
    print_matrix(matrix_1,n,nULL);
    free_all(matrix_1, vector_f1, n);
    
    
    
    matrix_1 = copy_matrix(matrix, n);
    vector_f1 = copy_vector(f_vector, n);

    
    
    printf("\n\tВ результате обычного метода Гаусса\n");
    double *ans_x_1 = GAUSS_method(matrix_1, vector_f1, n, 1);
    printf("Вектор х получился\n");
    print_vector(ans_x_1, n);
    free_all(matrix_1, vector_f1, n);

    
    matrix_1 = copy_matrix(matrix, n);
    vector_f1 = copy_vector(f_vector, n);
    
    
    printf("\n\tВ результате метода Гаусса с выбором максимального элемента\n");
    double *ans_x_2 = GAUSS_method(matrix_1, vector_f1, n, 2);
    printf("Вектор х получился\n");
    print_vector(ans_x_2, n);
    free_all(matrix_1, vector_f1, n);

    free_all(matrix, f_vector, n);

    return 0;
}