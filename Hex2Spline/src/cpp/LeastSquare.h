#ifndef LEASTSQUARE_H
#define LEASTSQUARE_H

#include <vector>
#include <array>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "BasicDataStructure.h"
//#include "TTSP_3D.h"

using namespace std;
using namespace Eigen;

class LeastSquare
{
public:
	LeastSquare();
	void GaussInfo(int ng=3);
	void SetProblem(const vector<int>& IDBC_in, const vector<double>& gh_in);
	void SetBoundary();
	void BasisFunction(double u, double v, double w, const vector<array<double, 3>>& pt, double Nx[64], double dNdx[64][3], double& detJ);
	void BasisFunctionMap(double Nx0[64], const vector<vector<double>>& cmat, vector<double>& Nx);
	void ElementMatrix(double Nx[64], double detJ, double EK[64][64]);
	void ElementMatrix(const vector<double>& Nx, vector<vector<double>>& EK);
	void ElementMatrix_Laplace(double dNdx[64][3], double detJ, double EK[64][64]);
	void ElementForce(double Nx[64], double Jmod, double Fb, double EF[64]);
	void ElementForce(const vector<double>& Nx, double Fb, vector<double>& EF);
	void Assembly(double EK[64][64], double EF[64], const vector<vector<double>>& cmat, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF);
	void Assembly(const vector<vector<double>>& EK, const vector<double>& EF, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF);
	void Assembly_Fitting(double EK[64][64], double EF[64], const vector<vector<double>>& cmat, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF);
	void BuildLinearSystem(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF);
	void BuildLinearSystem_Fitting(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF);
	void BuildLinearSystem_Fitting_C012(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF);
	void Solver(SparseMatrix<double>& GK, VectorXd& GF);

	void Para2Phys(double u, double v, double w, const vector<array<double, 3>>& bzpt, double Nx[64], array<double,3>& pt);
	//void Para2Phys(double u, double v, const double cpt[25][3], double pt[3]);
	void DispCal(double u,double v,double w,const BezierElement3D& bzel,double& disp,double& detJ);
	void ElementError(const BezierElement3D& bzel, double& L2, double& H1);
	void VisualizeVTK(const vector<BezierElement3D>& bzmesh, string fn);
	void VisualizeVTK_1(const vector<BezierElement3D>& bzmesh, string fn);//with smooth boundary
	void VisualizeVTK_2(const vector<BezierElement3D>& bzmesh, string fn);//with smooth boundary, remove elemenets
	void VisualizeBoundarySurface(const vector<BezierElement3D>& bzmesh, const vector<array<double, 3>>& cpts, string fn);
	void VisualizeError(const vector<BezierElement3D>& bzmesh, const vector<double>& err, string fn);
	//void ElementErrorEstimate(int eid,double& err,double& area,vector<int>& pid);
	void ErrorEstimate(const vector<BezierElement3D>& bzmesh, vector<double>& err);
	void Run(const vector<BezierElement3D>& bzmesh, string fn, vector<double>& err);
	void Run_Fitting(const vector<BezierElement3D>& bzmesh, string fn, vector<double>& sol);
	void Run_Fitting_C012(const vector<BezierElement3D>& bzmesh, string fn, vector<double>& sol);//blended B-splines

	double ElementErrorEstimate(const BezierElement3D& bzel);//based on residual and bubble function, return energy norm (square)
	void BasisFunctionMore(double u, double v, double w, const BezierElement3D& bzel, double pt[3], vector<double>& Nx, vector<array<double, 3>>& dNdx, double& detJ, double dtdx1[3][3]);
	void BubbleFunction(double u, double v, double w, double dtdx[3][3], double& f, double dfdx[3]);
	void DispFunction(const vector<int>& IEN, const vector<array<double, 3>>& dNdx, double dfdx[3]);
	void TotalErrorEstimate(const vector<BezierElement3D>& bzmesh, vector<double>& err);

	void OutputMatrix(const SparseMatrix<double>& GK, string fn);

	void SetProblem_SurfaceFitting(const vector<int>& IDBC_in);
	void Run_SurfaceFitting(const vector<BezierElement3D>& bzmesh, string fn, vector<array<double,3>>& sol);
	void BuildLinearSystem_SurfaceFitting(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK, vector<VectorXd>& GF);
	void ElementForce_SurfaceFitting(const vector<double>& Nx, double Fb[3], vector<vector<double>>& EF);
	void Para2Phys_SurfaceFitting(double u, double v, double w, const vector<array<double, 3>>& cnpt, array<double, 3>& pt);
	void Assembly_SurfaceFitting(const vector<vector<double>>& EK, const vector<vector<double>>& EF, const vector<int>& IEN, SparseMatrix<double>& GK, vector<VectorXd>& GF);
	void Solver_SurfaceFitting(SparseMatrix<double>& GK, vector<VectorXd>& GF);

	//void BasisFunction_TSP(double u, double v, const BezierElement3D& bzel, vector<double>& Nx, vector<array<double,2>>& dNdx, double& detJ);
	//void ElementMatrix_TSP(const vector<array<double,2>>& dNdx, double detJ, vector<vector<double>>& EK);
	//void ElementForce_TSP(const vector<double>& Nx, double detJ, double Fb, vector<double>& EF);
	//void Assembly_TSP(const vector<vector<double>>& EK, const vector<double>& EF, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF);
	//void BuildLinearSystem_TSP(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF);
	//void Para2Phys_TSP(double u, double v, const BezierElement3D& bzel, array<double,3>& pt);
	//void DispCal_TSP(double u,double v,const BezierElement3D& bzel,double& disp,double& detJ);
	//void Quantityh_TSP(double u,double v,const BezierElement3D& bzel,double val[3],double& detJ);

	void FindRange(const vector<BezierElement3D>& bzmesh);
	//void setProblem_ptestHO(const TruncatedTspline_3D& tts);
	double f_source_ls(const array<double, 3>& pt);
	//double f_source_1_ls(const array<double,3>& pt);
	//double f_source_2_ls(const array<double, 3>& pt);
	//double f_source_3_ls(const array<double, 3>& pt);
	//double f_source_4_ls(const array<double, 3>& pt);
	//double f_source_5_ls(const array<double, 3>& pt);
	//double f_source_6_ls(const array<double, 3>& pt);
	double f_source_7_ls(const array<double, 3>& pt);
	//double f_source_8_ls(const array<double, 3>& pt);
	double exact_sol_ls(const array<double, 3>& x);
	//double exact_sol_1_ls(const array<double, 3>& x);
	//double exact_sol_2_ls(const array<double, 3>& x);
	//double exact_sol_3_ls(const array<double, 3>& x);
	//double exact_sol_4_ls(const array<double, 3>& x);
	//double exact_sol_5_ls(const array<double, 3>& x);
	//double exact_sol_6_ls(const array<double, 3>& x);
	double exact_sol_7_ls(const array<double, 3>& x);
	//double exact_sol_8_ls(const array<double, 3>& x);

	void GetRemoveRegion(double xy[3][2]);
	void GetEqParameter(double xy[3][2], double nrm[3], double a);
	
private:
	vector<double> Gpt;
	vector<double> wght;
	vector<int> IDBC;
	vector<int> BCList;
	vector<double> gh;
	vector<array<double,3>> gh3;
	vector<double> uh;
	vector<array<double,3>> uh3;
	int npt;
	int neq;
	double range[3][2];
	double rmv[3][2];

	//used for solution
	double nmpl[3];
	double acoef;
	double bcoef;
	double dmlen[3];
	double xorg[3];
};

//double LDomainSolution(double x,double y);
//void LDomainSolutionDeriv(double x, double y, double& u, double& ux, double& uy);

//double ptestHO_sol_1(double x, double y);
//double ptestHO_sol_2(double x, double y);
//double ptestHO_sol_3(double x, double y);

//double f_source_1(const array<double,3>& x);//u=x+y, f=0
//double f_source_2(double x, double y);//u=x^2+y^2, f=-4
//double f_source_3(double x, double y);//u=x^3+y^3, f=-6(x+y)

#endif