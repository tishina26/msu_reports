#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double p(double x) {
	return 2;
}
double q(double x) {
	return (- 1 / x);
}
double f(double x) {
	return 3;
}

int main(int argc, char *argv[]) {
	double ab[] = {0.2, 0.5};
	double a = ab[0], b = ab[1];

	int n;
	printf("Введите частоту разбиения: \n");
	scanf("%d", &n);
	
	double h = (b - a) / n;
	
	double data1[] = {1, 0, 2};
	double data2[] = {0.5, -1, 1};
	double sigma1 = data1[0], gamma1 = data1[1], delta1 = data1[2];
	double sigma2 = data2[0], gamma2 = data2[1], delta2 = data2[2];

	double alpha[n], beta[n];
	alpha[0] = - gamma1 / (sigma1 * h - gamma1);
	beta[0] = delta1 * h / (sigma1 * h - gamma1);
	
	double x = a + h;
	
	for (int i = 1; i < n; i++) {
		double B = 1 / (h * h) + p(x) / (2 * h);
		double A = 1 / (h * h) - p(x) / (2 * h);
		double C = 2 / (h * h) - q(x);
		alpha[i] = B / (C - A * alpha[i - 1]);
		beta[i] = (A * beta[i - 1] - f(x)) / (C - A * alpha[i - 1]);
		x += h;
	}

	double y = (delta2 * h + gamma2 * beta[n - 1]) / 
				(sigma2 * h + gamma2 * (1 - alpha[n - 1]));
	
	printf("Точки полученные методом прогонки решения краевой задачи:\n");
	
	for (int i = n - 1; i >= 0; i--) {
		printf("(%.3lf; %.3lf)\n", x, y);
		x -= h;
		y = alpha[i] * y + beta[i];
	}
	printf("(%.3lf; %.3lf)\n", x, y);

	return 0;
}