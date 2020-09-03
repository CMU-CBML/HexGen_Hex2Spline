#include "PDEdata.h"
#include <cmath>

#define PI 3.141592654

double f_source(const array<double, 3>& x)
{
	//double tmp = f_source_1(x);
	double tmp = f_source_5(x);
	//double tmp = f_source_6(x);
	//double tmp = f_source_7(x);
	//double tmp = f_source_8(x);
	//double tmp = f_source_9(x);
	//double tmp = f_source_10(x);
	return tmp;
}

double f_source_1(const array<double, 3>& x)
{
	return 0.;
}

double f_source_2(const array<double, 3>& x)//analysis domain is a unit cube
{
	double tmp1 = -1600.*exp(-654. + 800.*x[0] * (1. - x[0]) + 800.*x[1] * (1. - x[1]) + 800.*x[2] * (1. - x[2])) / 3.;
	double tmp2 = exp((7. - 20.*x[0])*(7. - 20.*x[0]) + (7. - 20.*x[1])*(7. - 20.*x[1]) + (7. - 20.*x[2])*(7. - 20.*x[2]));
	double tmp3 = exp((13. - 20.*x[0])*(13. - 20.*x[0]) + (13. - 20.*x[1])*(13. - 20.*x[1]) + (13. - 20.*x[2])*(13. - 20.*x[2]));
	double tmp4 = 1011. - 1040.*x[0] + 800.*x[0] * x[0] - 1040.*x[1] + 800.*x[1] * x[1] - 1040.*x[2] + 800.*x[2] * x[2];
	double tmp5 = 291. - 560.*x[0] + 800.*x[0] * x[0] - 560.*x[1] + 800.*x[1] * x[1] - 560.*x[2] + 800.*x[2] * x[2];
	return tmp1*(tmp2*tmp4 + tmp3*tmp5);
}

double f_source_3(const array<double, 3>& x)//analysis domain is a unit cube
{
	double tmp1 = -1600.*exp(-300. + 400.*x[0] * (1. - x[0]) + 400.*x[1] * (1. - x[1]) + 400.*x[2] * (1. - x[2])) / 3.;
	double tmp2 = 597. - 800.*x[0] * (1. - x[0]) - 800.*x[1] * (1. - x[1]) - 800.*x[2] * (1. - x[2]);
	return tmp1*tmp2;
}

double f_source_4(const array<double, 3>& x)//analysis domain is a unit cube
{
	double tmp0 = 1. + 20.4124*x[0] + 20.4124*x[1] - 40.8248*x[2];
	double tmp1 = 1. / cosh(tmp0);
	return 5000.*tmp1*tmp1*tanh(tmp0);
}

double f_source_5(const array<double, 3>& x)//analysis domain is arbitrary
{
	double tmp = -10.*(x[0] + x[1] + x[2]);
	return tmp;
}

double f_source_6(const array<double, 3>& x)//analysis domain is arbitrary
{
	double tmp = 2.*(x[0] * (1. - x[0])*x[1] * (1. - x[1]) + x[0] * (1. - x[0])*x[2] * (1. - x[2]) + x[2] * (1. - x[2])*x[1] * (1. - x[1]));
	return tmp;
}

double f_source_7(const array<double, 3>& x)
{
	double acoef, bcoef, nmpl[3], xorg[3];
	//double tmp0 = acoef*(nmpl[0] * (x[0] - range[0][0]) + nmpl[1] * (x[1] - range[1][0]) + nmpl[2] * (x[2] - range[2][0]));
	double tmp0 = acoef*(nmpl[0] * (x[0] - xorg[0]) + nmpl[1] * (x[1] - xorg[1]) + nmpl[2] * (x[2] - xorg[2]));
	double tmp1 = 1. / cosh(tmp0);
	return bcoef*tmp1*tmp1*tanh(tmp0);
}

double f_source_8(const array<double, 3>& x)
{
	double range[3][2], dmlen[3], acoef, nmpl[3];
	double x1[3] = { (x[0] - range[0][0]) / dmlen[0], (x[1] - range[1][0]) / dmlen[1], (x[2] - range[2][0]) / dmlen[2] };
	double len2[3] = { dmlen[0] * dmlen[0], dmlen[1] * dmlen[1], dmlen[2] * dmlen[2] };
	double tmp0 = -acoef / exp(acoef*((x1[0] - nmpl[0])*(x1[0] - nmpl[0]) + (x1[1] - nmpl[1])*(x1[1] - nmpl[1]) + (x1[2] - nmpl[2])*(x1[2] - nmpl[2])));
	double tmp1 = (-2. + acoef)*(1. / len2[0] + 1. / len2[1] + 1. / len2[2]);
	double tmp2 = 4.*acoef*(x1[0] * (x1[0] - 1.) / len2[0] + x1[1] * (x1[1] - 1.) / len2[1] + x1[2] * (x1[2] - 1.) / len2[2]);
	return tmp0*(tmp1 + tmp2);
}

double f_source_9(const array<double, 3>& x)
{
	//double PI(3.141592654);
	return 3.*PI*PI*sin(PI*x[0]) * sin(PI*x[1]) * sin(PI*x[2]);
}

double f_source_10(const array<double, 3>& x)
{
	return -exp((x[0] + x[1] + x[2]) / 3.) / 3.;
}



double exact_sol(const array<double, 3>& x)
{
	//double tmp = exact_sol_1(x);
	//double tmp = exact_sol_1y(x);
	//double tmp = exact_sol_1z(x);
	double tmp = exact_sol_5(x);
	//double tmp = exact_sol_6(x);
	//double tmp = exact_sol_7(x);
	//double tmp = exact_sol_8(x);
	//double tmp = exact_sol_9(x);
	//double tmp = exact_sol_10(x);
	return tmp;
}

double exact_sol_1(const array<double, 3>& x)
{
	return x[0];
}

double exact_sol_1y(const array<double, 3>& x)
{
	return x[1];
}

double exact_sol_1z(const array<double, 3>& x)
{
	return x[2];
}

double exact_sol_2(const array<double, 3>& x)
{
	double tmp1 = 1. / exp((20.*x[0] - 7.)*(20.*x[0] - 7.) + (20.*x[1] - 7.)*(20.*x[1] - 7.) + (20.*x[2] - 7.)*(20.*x[2] - 7.));
	double tmp2 = 1. / exp((20.*x[0] - 13.)*(20.*x[0] - 13.) + (20.*x[1] - 13.)*(20.*x[1] - 13.) + (20.*x[2] - 13.)*(20.*x[2] - 13.));
	return 2.*(tmp1 + tmp2) / 3.;
}

double exact_sol_3(const array<double, 3>& x)
{
	return 2. / (3.*exp((20.*x[0] - 10.)*(20.*x[0] - 10.) + (20.*x[1] - 10.)*(20.*x[1] - 10.) + (20.*x[2] - 10.)*(20.*x[2] - 10.)));
}

double exact_sol_4(const array<double, 3>& x)
{
	return tanh(1. - 50.*(-0.408248*x[0] - 0.408248*x[1] + 0.816497*x[2]));
}

double exact_sol_5(const array<double, 3>& x)
{
	double tmp = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2])*(x[0] + x[1] + x[2]) + x[0] * x[1] * x[2];
	return tmp;
}

double exact_sol_6(const array<double, 3>& x)
{
	double tmp = x[0] * (1. - x[0])*x[1] * (1. - x[1])*x[2] * (1. - x[2]);;
	return tmp;
}

double exact_sol_7(const array<double, 3>& x)
{
	double xorg[3], acoef, nmpl[3];
	//double x1[3] = { (x[0] - range[0][0]), (x[1] - range[1][0]), (x[2] - range[2][0]) };
	double x1[3] = { (x[0] - xorg[0]), (x[1] - xorg[1]), (x[2] - xorg[2]) };
	return tanh(acoef*(nmpl[0] * x1[0] + nmpl[1] * x1[1] + nmpl[2] * x1[2]));
}

double exact_sol_8(const array<double, 3>& x)
{
	double range[3][2], acoef, nmpl[3], dmlen[3];
	double x1[3] = { (x[0] - range[0][0]) / dmlen[0], (x[1] - range[1][0]) / dmlen[1], (x[2] - range[2][0]) / dmlen[2] };
	double tmp = 1. / exp(acoef*((x1[0] - nmpl[0])*(x1[0] - nmpl[0]) + (x1[1] - nmpl[1])*(x1[1] - nmpl[1]) + (x1[2] - nmpl[2])*(x1[2] - nmpl[2])));
	return 1. / exp(acoef*((x1[0] - nmpl[0])*(x1[0] - nmpl[0]) + (x1[1] - nmpl[1])*(x1[1] - nmpl[1]) + (x1[2] - nmpl[2])*(x1[2] - nmpl[2])));
}

double exact_sol_9(const array<double, 3>& x)
{
	//double PI(3.141592654);
	return sin(PI*x[0]) * sin(PI*x[1]) * sin(PI*x[2]);
}

double exact_sol_10(const array<double, 3>& x)
{
	return exp((x[0] + x[1] + x[2]) / 3.);
}



void grad_sol(const array<double, 3>& x, array<double,3>& del_u)
{
	//grad_sol_1(x, del_u);
	//grad_sol_1y(x, del_u);
	//grad_sol_1z(x, del_u);
	grad_sol_5(x, del_u);
	//grad_sol_9(x, del_u);
	//grad_sol_10(x, del_u);
}

void grad_sol_1(const array<double, 3>& x, array<double, 3>& del_u)
{
	del_u[0] = 1.;
	del_u[1] = 0.;
	del_u[2] = 0.;
}

void grad_sol_1y(const array<double, 3>& x, array<double, 3>& del_u)
{
	del_u[0] = 0.;
	del_u[1] = 1.;
	del_u[2] = 0.;
}

void grad_sol_1z(const array<double, 3>& x, array<double, 3>& del_u)
{
	del_u[0] = 0.;
	del_u[1] = 0.;
	del_u[2] = 1.;
}

void grad_sol_5(const array<double, 3>& x, array<double, 3>& del_u)
{
	double sum2(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]), sum(x[0] + x[1] + x[2]);
	del_u[0] = sum2 + x[1] * x[2] + 2.*x[0] * sum;
	del_u[1] = sum2 + x[0] * x[2] + 2.*x[1] * sum;
	del_u[2] = sum2 + x[0] * x[1] + 2.*x[2] * sum;
}

void grad_sol_9(const array<double, 3>& x, array<double, 3>& del_u)
{
	del_u[0] = PI*cos(PI*x[0])*sin(PI*x[1])*sin(PI*x[2]);
	del_u[1] = PI*sin(PI*x[0])*cos(PI*x[1])*sin(PI*x[2]);
	del_u[2] = PI*sin(PI*x[0])*sin(PI*x[1])*cos(PI*x[2]);
}

void grad_sol_10(const array<double, 3>& x, array<double, 3>& del_u)
{
	del_u[0] = exp((x[0] + x[1] + x[2]) / 3.) / 3.;
	del_u[1] = exp((x[0] + x[1] + x[2]) / 3.) / 3.;
	del_u[2] = exp((x[0] + x[1] + x[2]) / 3.) / 3.;
}