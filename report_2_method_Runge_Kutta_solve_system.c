#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double f(double x, double y) {
    return (x - x * x) * y;
}
double f1(double x, double u, double v) {
	return cos(u + 1.1 * v) + 2.1;
}
double f2(double x, double u, double v) {
	return 1.1 / (x + 2.1 * u * u) + x + 1;
}

int main(int argc, char **argv) {
	double a, b;
	printf("Задайте отрезок, на котором строить график\n");
	scanf("%lf %lf", &a, &b);

	int n;
	printf("Задайте частоту разбивки для корректного применения методов\n");
	scanf("%d", &n);

	double h = (b - a) / n;
	
	// Алгоритм Рунге-Кутта второго порядка точности для одного уравнения

	double x0 = 0, y0 = 1;
	
	printf("\nАлгоритм R_K_2 для 1 уравнения дает точки:\n");
	
	double x_i = x0, y_i = y0;

	for (int i = 0; i < n; i++) {
		y_i = y_i + h/2. * (f(x_i,y_i) + f(x_i + h, y_i + h*f(x_i,y_i)));
		x_i += h;
		printf("(%.3lf;\t%.3lf)\n", x_i, y_i);
	}

	// Алгоритм Рунге-Кутта четвертого порядка точности для одного уравнения
	printf("\nАлгоритм R_K_4 для 1 уравнения дает точки:\n");

	x_i = x0; y_i = y0;

	for (int i = 0; i < n; i++) {
		double k1 = h*f(x_i,y_i);
		double k2 = h*f(x_i + h/2.,y_i + k1/2.);
		double k3 = h*f(x_i + h/2.,y_i + k2/2.);
		double k4 = h*f(x_i + h,y_i + k3);
		x_i += h;
		y_i = y_i + (k1 + 2*k2 + 2*k3 + k4) / 6.;
		printf("(%.3lf;\t%.3lf)\n", x_i, y_i);
	}
	
	// Алгоритм Рунге-Кутта второго порядка точности для системы из 2 уравнений
	
	x0 = 0;
	double y10 = 1, y20 = 0.05;
	
	printf("\nАлгоритм R_K_2 для системы уравнений дает точки:\n");
	
	x_i = x0;
	double u_i = y10, v_i = y20, u, v;

	for (int i = 0; i < n; i++) {
		double du = u_i + h*f1(x_i,u_i,v_i);
		double dv = v_i + h*f2(x_i,u_i,v_i);
		u = u_i + (f1(x_i,u_i,v_i) + f1(x_i + h,du,dv))*h/2;
		v = v_i + (f2(x_i,u_i,v_i) + f2(x_i + h,du,dv))*h/2;
		x_i += h;
		u_i = u; v_i = v;
		printf("(%.3lf,\t%.3lf,\t%.3lf)\n", x_i, u_i, v_i);
	}

	//Алгоритм Рунге-Кутта четвертого порядка точности для системы из 2 уравнений
	printf("\nАлгоритм R_K_4 для системы уравнений дает точки:\n");

	x_i = x0, u_i = y10, v_i = y20;
	
	for (int i = 0; i < n; i++) {
		double k1 = h*f1(x_i,u_i,v_i);
		double l1 = h*f2(x_i,u_i,v_i);
		double k2 = h*f1(x_i + h/2.,u_i + k1/2., v_i + l1/2.);
		double l2 = h*f2(x_i + h/2.,u_i + k1/2., v_i + l1/2.);
		double k3 = h*f1(x_i + h/2.,u_i + k2/2., v_i + l2/2.);
		double l3 = h*f2(x_i + h/2.,u_i + k2/2., v_i + l2/2.);
		double k4 = h*f1(x_i + h,u_i + k3, v_i + l3);
		double l4 = h*f2(x_i + h,u_i + k3, v_i + l3);
		x_i += h;
		u_i = u_i + (k1 + 2*k2 + 2*k3 + k4) / 6.;
		v_i = v_i + (l1 + 2*l2 + 2*l3 + l4) / 6.;
		printf("(%.3lf,\t%.3lf,\t%.3lf)\n", x_i, u_i, v_i);
	}
	return 0;
}