#ifndef PDE_DATA_H
#define PDE_DATA_H

#include <array>

using namespace std;

double f_source(const array<double, 3>& pt);
double f_source_1(const array<double, 3>& pt);
double f_source_2(const array<double, 3>& pt);
double f_source_3(const array<double, 3>& pt);
double f_source_4(const array<double, 3>& pt);
double f_source_5(const array<double, 3>& pt);
double f_source_6(const array<double, 3>& pt);
double f_source_7(const array<double, 3>& pt);
double f_source_8(const array<double, 3>& pt);
double f_source_9(const array<double, 3>& pt);
double f_source_10(const array<double, 3>& pt);

double exact_sol(const array<double, 3>& x);
double exact_sol_1(const array<double, 3>& x);
double exact_sol_1y(const array<double, 3>& x);
double exact_sol_1z(const array<double, 3>& x);
double exact_sol_2(const array<double, 3>& x);
double exact_sol_3(const array<double, 3>& x);
double exact_sol_4(const array<double, 3>& x);
double exact_sol_5(const array<double, 3>& x);
double exact_sol_6(const array<double, 3>& x);
double exact_sol_7(const array<double, 3>& x);
double exact_sol_8(const array<double, 3>& x);
double exact_sol_9(const array<double, 3>& x);
double exact_sol_10(const array<double, 3>& x);

void grad_sol(const array<double, 3>& x, array<double, 3>& del_u);
void grad_sol_1(const array<double, 3>& x, array<double, 3>& del_u);
void grad_sol_1y(const array<double, 3>& x, array<double, 3>& del_u);
void grad_sol_1z(const array<double, 3>& x, array<double, 3>& del_u);
void grad_sol_5(const array<double, 3>& x, array<double, 3>& del_u);
void grad_sol_9(const array<double, 3>& x, array<double, 3>& del_u);
void grad_sol_10(const array<double, 3>& x, array<double, 3>& del_u);

#endif