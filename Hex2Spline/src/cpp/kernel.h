#ifndef KERNEL_H
#define KERNEL_H

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include "KnotInsertion.h"
//#include "TTSP_2D.h"
#include "TTSP_3D.h"
#include "Elasiticity_3D.h"
#include "Laplace.h"
#include "LeastSquare.h"
#include "atlstr.h" 
#include <sstream>
#include <iomanip>
using namespace std;

class kernel
{
public:
	void run();

	void run_complex();
	void run_complex_fit();
	void run_complex_glb();
	void run_Bspline();
	void run_Bspline_fit();

	//convergence 3D study
	void run_conv();//solid T-splines
	void run_coupling();//generalized Bezier
	void run_coupling_Bezier();//explicitly convert spline to Bezier
	void run_coupling_comp();//comparison using B-splines-like

	void run_UB_adapt_fit();//built without open knot vectors, problem
	void run_THS3D();//built with open knot vectors
	void run_HS3D();//without truncation, rationalization instead

	void run_C0C1Bezier();//C1 functions with necessary C0 Bezier functions
	void run_C0Bezier();//All C0 Bezier functions for comparison
	void run_C02iBezier();//improvement of run_C0Bezier

	void run_UBsplines();//B-splines-like to handle unstructured
	void run_adaptive_fit();

	void run_benchmark();
	void run_leastsquare();
	void run_leastsquare_1();
	void run_pipeline();
	void output_err(string fn, const vector<int>& dof, const vector<double>& err);
	void output_err(string fn, const vector<double>& dof, const vector<double>& err);

	void run_MeshQualityImprove();
	//Abaqus with GEM
	void run_AbaqusGEM();

	// For User Friendly
	void run_coupling_comp(int flag_sharp, int ngref, int nlref, double tol_sharp, string fn);//comparison using B-splines-like
	void run_MeshQualityImprove(int mode, int flag_sharp, double tol_sharp, int opt_par1, double opt_par2, string fn);

//private:
//	TruncatedTspline_3D tt3;
//	LinearElasticity le;
//	Laplace lap;
};

#endif