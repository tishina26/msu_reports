#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>


extern double f1(double x);
extern double f2(double x);
extern double f3(double x);
extern double df1(double x);
extern double df2(double x);
extern double df3(double x);

//_____________________________________
#define a12 0.0
#define b12 1.0
#define a13 -1.9
#define b13 -1.5
#define a23 -1.0
#define b23 0.0

static int count_steps_root;
static int count_steps[3];
static double root_arr[3];
static double integral_arr[3];

static double eps1 = 0.0001;
static double eps2 = 0.0001;
//_____________________________________


// вспомогательная функция, вычисляющая значение F(x) = f(x) - g(x)
double F(double(*f)(double), double(*g)(double), double x) {
	return f(x) - g(x);
}

// for methods choose right or left
int find_case(double(*f)(double), double(*g)(double), 
			  double a, double b) {
	// F = f - g
	// F убывает при F(a) >= 0
	double Fa = F(f, g, a);

	// F выше хорды при 
	// F(0.5a+0.5b) > 0.5 (F(a) + F(b));
	double Fc = F(f, g, (a + b) / 2);
	double cF = 0.5 * (F(f, g, a) + F(f, g, b));

	// случай 1 если F убывает и выше хорды
	if ((Fa > 0) && (Fc > cF)) {
		return 1;
	}

	// случай 1 если F возрастает и ниже хорды
	if ((Fa < 0) && (Fc < cF)) {
		return 1;
	}

	return 2;
}

// метод хорд
void chord_method(double(*f)(double), double(*g)(double),
				 double a, double b, double *new_a, double *new_b) {
	// c = aF(b) - bF(a)    /    F(b) - F(a)
	// case 1 => [c, b]   case 2 => [a, c]
	double c1 = a * F(f, g, b) - b * F(f, g, a);
	double c2 = F(f, g, b) - F(f, g, a);

	if (find_case(f, g, a, b) == 1) {
		(*new_a) = c1 / c2;
	} else {
		(*new_b) = c1 / c2;
	}
}

// метод касательных
void tangent_method(double(*f)(double), double(*g)(double),
				  double(*df)(double), double(*dg)(double), 
				  double a, double b, double *new_a, double *new_b) {
	// c = d - F(d) / dF(d)
	// case 1 => d = b   case 2 => d = a
	int tmp_case = find_case(f, g, a, b);

	double d = b;
	if (tmp_case == 2) {
		d = a;
	}

	double c = d - F(f, g, d) / F(df, dg, d);

	if (tmp_case == 1) {
		(*new_b) = c;
	} else {
		(*new_a) = c;
	}
}

// поиск корня с помощью комбинированного метода
double root(double(*f)(double), double(*g)(double),
			double(*df)(double), double(*dg)(double), 
			double a, double b, double eps) {
	if (F(f, g, a) * F(f, g, b) > 0.0) {
		printf("неверно задан отрезок [a, b]\nпопробуйте другой\n");
		return 0;
	}

	if (b - a < eps1) {
		return (a + b) / 2;
	}

	double new_a = a, new_b = b;  // because must not change after 1 method

	chord_method(f, g, a, b, &new_a, &new_b);
	tangent_method(f, g, df, dg, a, b, &new_a, &new_b);

	count_steps_root++;

	return root(f, g, df, dg, new_a, new_b, eps);
}

double Sympson_integral(double(*f)(double), double a, double b, int n, double* arr, int index) {
	double h = (b - a) / n;
	(arr[0]) = f(a);
	(arr[n]) = f(b);
	double ans = (arr[0]) + (arr[n]);

	for (int i = 1; i < n; i++) {
		// заполняем массив
		if (index == 1) {
			(arr[2 * i]) = f(a + i * h);
		} else {
			(arr[i]) = f(a + i * h);
		}

		ans += (2 + 2 * (i % 2)) * (arr[i]);
	}

	return ans * h / 3;
}

int count = 0;
int n = 90;
double integral(double(*f)(double), double a, double b, double eps) {
	// (Int_n - Int_2n) / 15 <= eps
	// возьмем n = 20 (10 недостаточно)

	// чтобы не перевычислять значения int_n повторня для int_2n
	double* arr_F = (double*)calloc(n*2+1, sizeof(double));

	// последняя цифра нужна, чтобы понять, заполнять массив(1) или нет(0)
	double Int_n = Sympson_integral(f, a, b, n, arr_F, 1);
	double Int_2n =  Sympson_integral(f, a, b, 2 * n, arr_F, 0);

	free(arr_F);
	// 15 - из правила Рунге
	if ((Int_n - Int_2n) / 15 <= eps) {
		return Int_2n;
	}

	if (count >= 15) {
		printf("%d", count);
		printf("слишком большие границы интегрирования, не достигнута требуемая точность\n");
		printf("или неправильно введены границы\n");
		return 0;
	}

	count++;
	n *= 2;

	return integral(f, a, b, eps);
}

// считает площадь под графиком
double area(void) {
	count_steps_root = 0;
	double x1_2 = root(f1, f2, df1, df2, a12, b12, eps1);
	count_steps[0] = count_steps_root;
	count_steps_root = 0;
	double x1_3 = root(f1, f3, df1, df3, a13, b13, eps1);
	count_steps[1] = count_steps_root;
	count_steps_root = 0;
	double x2_3 = root(f2, f3, df2, df3, a23, b23, eps1);
	count_steps[2] = count_steps_root;

	root_arr[0] = x1_2, root_arr[1] = x1_3, root_arr[2] = x2_3;
	integral_arr[0] = integral(f1, x1_3, x1_2, eps2);
	integral_arr[1] = integral(f2, x2_3, x1_2, eps2);
	integral_arr[2] = integral(f3, x1_3, x2_3, eps2);
	return integral_arr[0] - integral_arr[1] - integral_arr[2];
}

void help(void) {
	printf("-constants          \tвывести отрезки для поиска корня и eps\n");
	printf("-test_root         \tвызывает root(index1, index2, a, b, eps) - параметры вводятся вручную\n");
	printf("\t\t\t   index1, index2 - номера функций, у которых надо найти точку пересечения\n");
	printf("\t\t\t   a, b - отрезок, в котором надо искать корень (можете воспользоваться -constants)\n");
	printf("\t\t\t   eps - с какой точностью надо найти корень\n");
	printf("-test_integral     \tвызывает integral(index1, a, b, eps) - параметры вводятся вручную\n");
	printf("\t\t\t   index - номер функции, у которой надо найти интеграл\n");
	printf("\t\t\t   a, b - пределы интегрирования (можете воспользоваться -test_root)\n");
	printf("\t\t\t   eps - с какой точностью надо найти интеграл\n");
	printf("-test_root_auto     \tвызывает root(от 3 разных наборов параметров, подобранных автором)\n");
	printf("\t\t\t   и сравнивает ответ программы с правильным ответом\n");
	printf("-test_integral_auto \tвызывает integral(от 3 разных наборов параметров, подобранных автором)\n");
	printf("\t\t\t   и сравнивает ответ программы с правильным ответом\n");
	printf("-steps             \tвыводит корни и количество шагов, потраченных на их нахождение\n");
}

void constants(void) {
	printf("Границы:\n");
	printf("  1_2: [%lf, %lf]\n", a12, b12);
	printf("  1_3: [%lf, %lf]\n", a13, b13);
	printf("  2_3: [%lf, %lf]\n", a23, b23);
	printf("\nEps1 = %lf\tEps2 = %lf\n", eps1, eps2);
}

void test_root(int i, char* arr[]) {
	i++;
	for (int j = i; j < i + 5; j++) {
		if (arr[j] && ((arr[j][0] >= '0' && arr[j][0] <= '9') || (arr[j][0] == '-' && (arr[j][1] && arr[j][1] >= '0' && arr[j][1] <= '9')))) {
			continue;
		} else {
			printf("введено мало аргументов или не числа\ntest_root(i1, i2, a, b, eps)\n");
		}	return;
	}

	int i1, i2;
	double a, b, eps;
	sscanf(arr[i], "%d", &i1);
	sscanf(arr[i+1], "%d", &i2);
	sscanf(arr[i+2], "%lf", &a);
	sscanf(arr[i+3], "%lf", &b);
	sscanf(arr[i+4], "%lf", &eps);

	if (i1 == i2) {
		printf("Нельзя задавать одинаковые индексы\n");
		return;
	}

	if (i2 == 3 && a <= -2 && b >= -2) {
		printf("Пожалуйста, введите для 3-й функции границы, не включающие число -2\n");
		return;
	}

	double (*inp_fun[])(double) = {0, f1, f2, f3, 0, df1, df2, df3};

	count_steps_root = 0;
	double x = root(inp_fun[i1], inp_fun[i2], inp_fun[i1+4], inp_fun[i2+4], a, b, eps);
	count_steps[i1 * (i2 % 2)] = count_steps_root;
	root_arr[i1 * (i2 % 2)] = x;
	printf("\troot f%d - f%d : %lf\n", i1, i2, x);
}

void test_root_auto(void) {
	char* res[] = {"OK", "NOT ok, check programm"};
	double a, b, eps, root_ans, right_ans;
	int index;
	//________________________________________________________
	a = -5.95, b = 6.68, eps = 0.01;
	printf("Test1: a = %lf, b = %lf, eps = %lf\n", a, b, eps);
	count_steps_root = 0;
	root_ans = root(f1, f2, df1, df2, a, b, eps);
	count_steps[0] = count_steps_root;
	root_arr[0] = root_ans;
	// right_ans = 0.448178;  // 0.448178 [0.1]
	right_ans = (79.0-sqrt(5289.0))/14.0;  // 0.448178 [0.1]
	printf("\tF = f1 - f2 => root = %lf\n", root_ans);
	index = 1;
	if (root_ans - right_ans <= eps && root_ans - right_ans>= -eps) {
		index = 0;
	}
	printf("\t\t%s\n\n", res[index]);
	//________________________________________________________
	a = -1.9, b = 5, eps = 0.001;
	printf("Test2: a = %lf, b = %lf, eps = %lf\n", a, b, eps);
	count_steps_root = 0;
	root_ans = root(f1, f3, df1, df3, a, b, eps);
	count_steps[1] = count_steps_root;
	root_arr[1] = root_ans;
	right_ans = -1.8211369;  // -1.821137 [-1.9, -1]     0 means no answer
	printf("\tF = f1 - f3 => root = %lf\n", root_ans);
	index = 1;
	if (root_ans - right_ans <= eps && root_ans - right_ans >= -eps) {
		index = 0;
	}
	printf("\t\t%s\n\n", res[index]);
	//________________________________________________________
	a = -1.999, b = 1000, eps = 0.0001;
	printf("Test3: a = %lf, b = %lf, eps = %lf\n", a, b, eps);
	count_steps_root = 0;
	root_ans = root(f2, f3, df2, df3, a, b, eps);
	count_steps[2] = count_steps_root;
	root_arr[2] = root_ans;
	right_ans = (-7+sqrt(37))/6;  // -0.152873 [-1, 0]
	printf("\tF = f2 - f3 => root = %lf\n", root_ans);
	index = 1;
	if (root_ans - right_ans <= eps && root_ans - right_ans>= -eps) {
		index = 0;
	}
	printf("\t\t%s\n\n", res[index]);
	//________________________________________________________
}

void test_integral(int i, char* arr[]) {
	i++;
	for (int j = i; j < i + 4; j++) {
		if (arr[j] && ((arr[j][0] >= '0' && arr[j][0] <= '9') || (arr[j][0] == '-' && (arr[j][1] && arr[j][1] >= '0' && arr[j][1] <= '9')))) {
			continue;
		} else {
			printf("аргументы введены неверно\ntest_integral(index, a, b, eps)\n");
		}	return;
	}

	int index;
	double a, b, eps;
	sscanf(arr[i], "%d", &index);
	if (index != 1 && index != 2 && index != 3) {
		printf("Неправильно введен индекс функции (попробуйте 1, 2 или 3)\ntest_integral(index, a, b, eps)\n");
		return;
	}
	sscanf(arr[i+1], "%lf", &a);
	sscanf(arr[i+2], "%lf", &b);
	sscanf(arr[i+3], "%lf", &eps);

	if (index == 3 && a <= -2 && b >= -2) {
		printf("Пожалуйста, введите для 3-й функции границы, не включающие число -2\n");
		return;
	}

	double (*inp_fun[])(double) = {0, f1, f2, f3};
	double integ = integral(inp_fun[index], a, b, eps);
	printf("\tintegral f%d [%lf; %lf]: %lf\n", index, a, b, integ);
	integral_arr[index - 1] = integ;
}

void test_integral_auto(void) {
	char* res[] = {"OK", "NOT ok, check programm"};
	double a, b, eps, integral_ans, right_ans;
	int index;
	//________________________________________________________
	a = 0.0, b = 12345.0, eps = 0.001;
	printf("Test1: a = %lf, b = %lf, eps = %lf\n", a, b, eps);
	integral_ans = integral(f1, a, b, eps);
	integral_arr[0] = integral_ans;
	right_ans = 219420339550.875;  // 8.322250
	printf("\tf1 [%lf; %lf] => integral = %lf\n", a, b, integral_ans);
	index = 1;
	if (integral_ans - right_ans <= eps && integral_ans - right_ans >= -eps) {
		index = 0;
	}
	printf("\t\t%s\n\n", res[index]);
	//________________________________________________________
	a = 7.95, b = 99, eps = 0.001;
	printf("Test2: a = %lf, b = %lf, eps = %lf\n", a, b, eps);
	integral_ans = integral(f2, a, b, eps);
	integral_arr[1] = integral_ans;
	right_ans = 14697.74625;  // 0.867290
	printf("\tf2 [%lf; %lf] => integral = %lf\n", a, b, integral_ans);
	index = 1;
	if (integral_ans - right_ans <= eps && integral_ans - right_ans >= -eps) {
		index = 0;
	}
	printf("\t\t%s\n\n", res[index]);
	//________________________________________________________
	a = 9.584, b = 51.453, eps = 0.01;
	printf("Test3: a = %lf, b = %lf, eps = %lf\n", a, b, eps);
	integral_ans = integral(f3, a, b, eps);
	integral_arr[2] = integral_ans;
	right_ans = log(53453.0/11584.0); // 2.334854
	printf("\tf3 [%lf; %lf] => integral = %lf\n", a, b, integral_ans);
	index = 1;
	if (integral_ans - right_ans <= eps && integral_ans - right_ans >= -eps) {
		index = 0;
	}
	printf("\t\t%s\n\n", res[index]);
	//________________________________________________________
}

void steps(void) {
	printf("Результаты приведены для последних поисков корней\n");
	printf("f1 - f2 = 0 => root = %lf;   steps = %d\n", root_arr[0], count_steps[0]);
	printf("f1 - f3 = 0 => root = %lf;   steps = %d\n", root_arr[1], count_steps[1]);
	printf("f2 - f3 = 0 => root = %lf;   steps = %d\n", root_arr[2], count_steps[2]);
	
	printf("Результаты после поиска площади:\n");
	area();
	printf("f1 - f2 = 0 => root = %lf;   steps = %d\n", root_arr[0], count_steps[0]);
	printf("f1 - f3 = 0 => root = %lf;   steps = %d\n", root_arr[1], count_steps[1]);
	printf("f2 - f3 = 0 => root = %lf;   steps = %d\n", root_arr[2], count_steps[2]);

	printf("\n\nдополнительные данные:\n\n");
	printf("f1 - f2 = 0 => (x; y) = (%lf; %lf)\n", root_arr[0], f1(root_arr[0]));
	printf("f1 - f3 = 0 => (x; y) = (%lf; %lf)\n", root_arr[1], f1(root_arr[1]));
	printf("f2 - f3 = 0 => (x; y) = (%lf; %lf)\n\n", root_arr[2], f3(root_arr[2]));

	printf("integral f1 [root; root] = %lf\n", integral_arr[0]);
	printf("integral f2 [root; root] = %lf\n", integral_arr[1]);
	printf("integral f3 [root; root] = %lf\n", integral_arr[2]);

}

void hello(void) {
	printf("\nФайл make предусматривает make all_flags компиляцию со всеми флагами\n\n");
	printf("_______________________________________________________________\n");
	printf("| Функции:                                                    |\n");
	printf("|     f1 = 0.35x^2-0.95x+2.7                                  |\n");
	printf("|     f2 = 3x+1                                               |\n");
	printf("|     f3 = 1/(x+2)                                            |\n");
	printf("| Флаг '-help' покажет, какие флаги предусмотрены программой  |\n");
	printf("|_____________________________________________________________|\n");
}

int main(int argc, char *argv[]) {
	static char* flags[7] = {"-help", "-constants", "-test_root", "-test_integral", "-test_root_auto", "-test_integral_auto", "-steps"};
	hello();

	printf("\n\n_______________________________________________\n");
	printf("    Answer: S = %.6lf |\n", area());
	printf("f1 - f2 = 0 => root = %lf;   steps = %d\n", root_arr[0], count_steps[0]);
	printf("f1 - f3 = 0 => root = %lf;   steps = %d\n", root_arr[1], count_steps[1]);
	printf("f2 - f3 = 0 => root = %lf;   steps = %d\n", root_arr[2], count_steps[2]);
	printf("_______________________________________________\n\n\n");

	if (argc > 1) {
		for(int i = 1; i < argc; i++) {
			if (strcmp(argv[i], flags[0]) == 0) {
				help();
			} else if (strcmp(argv[i], flags[1]) == 0) {
				constants();
			} else if (strcmp(argv[i], flags[2]) == 0) {
				test_root(i, argv);
				i += 5;
			} else if (strcmp(argv[i], flags[3]) == 0) {
				test_integral(i, argv);
				i += 4;
			} else if (strcmp(argv[i], flags[4]) == 0) {
				test_root_auto();
			} else if (strcmp(argv[i], flags[5]) == 0) {
				test_integral_auto();
			} else if (strcmp(argv[i], flags[6]) == 0) {
				steps();
			} else {
				printf("нет такого флага %s\n", argv[i]);
				printf("или введено неверное количество входных данных\n");
				printf("попробуйте '-help'\n");
			}
		}
	}

	return 0;
}