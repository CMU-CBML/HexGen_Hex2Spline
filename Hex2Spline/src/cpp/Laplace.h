#ifndef LAPLACE_H
#define LAPLACE_H

#include <vector>
#include <array>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "BasicDataStructure.h"
#include "TTSP_3D.h"

using namespace std;
using namespace Eigen;

class Laplace
{
public:
	Laplace();
	void GaussInfo(int ng=3);
	void SetProblem(const vector<int>& IDBC_in, const vector<double>& gh_in);
	void SetBoundary();
	void BasisFunction(double u, double v, double w, const vector<array<double, 3>>& pt, double Nx[64], double dNdx[64][3], double& detJ);
	void ElementMatrix(double dNdx[64][3], double detJ, double EK[64][64]);
	void ElementForce(double Nx[64], double Jmod, double Fb, double EF[64]);
	void Assembly(double EK[64][64], double EF[64], const vector<vector<double>>& cmat, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF);
	void BuildLinearSystem(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF);
	void Solver(SparseMatrix<double>& GK, VectorXd& GF);

	void Para2Phys(double u, double v, double w, const vector<array<double, 3>>& bzpt, double Nx[64], array<double,3>& pt);
	//void Para2Phys(double u, double v, const double cpt[25][3], double pt[3]);
	void DispCal(double u,double v,double w,const BezierElement3D& bzel,double& disp,double& detJ);
	void ElementError(const BezierElement3D& bzel, double& L2, double& H1);

	void VisualizeLSDYNA(const vector<BezierElement3D>& bzmesh, string fn);

	void VisualizeVTK(const vector<BezierElement3D>& bzmesh, string fn);
	void VisualizeVTK_1(const vector<BezierElement3D>& bzmesh, string fn);//with smooth boundary
	void VisualizeVTK_2(const vector<BezierElement3D>& bzmesh, string fn);//with smooth boundary, remove elemenets
	void VisualizeError(const vector<BezierElement3D>& bzmesh, const vector<double>& err, string fn);
	void VisualizeError_2(const vector<BezierElement3D>& bzmesh, const vector<double>& err, string fn);//remove elemenets
	//void ElementErrorEstimate(int eid,double& err,double& area,vector<int>& pid);
	void ErrorEstimate(const vector<BezierElement3D>& bzmesh, vector<double>& err);
	void Run(const vector<BezierElement3D>& bzmesh, string fn, vector<double>& err);

	double ElementErrorEstimate(const BezierElement3D& bzel);//based on residual and bubble function, return energy norm (square)
	void BasisFunctionMore(double u, double v, double w, const BezierElement3D& bzel, double pt[3], vector<double>& Nx, vector<array<double, 3>>& dNdx, double& detJ, double dtdx1[3][3]);
	void BubbleFunction(double u, double v, double w, double dtdx[3][3], double& f, double dfdx[3]);
	void DispFunction(const vector<int>& IEN, const vector<array<double, 3>>& dNdx, double dfdx[3]);
	void TotalErrorEstimate(const vector<BezierElement3D>& bzmesh, vector<double>& err);

	void OutputMatrix(const SparseMatrix<double>& GK, string fn);

	void BasisFunction_1(double u, double v, double w, const BezierElement3D& bzel, vector<double>& Nx, vector<array<double, 3>>& dNdx, double& detJ);
	void ElementMatrix_1(const vector<array<double, 3>>& dNdx, double detJ, vector<vector<double>>& EK);
	void ElementForce_1(const vector<double>& Nx, double Jmod, double Fb, vector<double>& EF);
	void Assembly_1(const vector<vector<double>>& EK, const vector<double>& EF, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF);
	void BuildLinearSystem_1(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF);
	void Assembly_Dense(const vector<vector<double>>& EK, const vector<double>& EF, const vector<int>& IEN, MatrixXd& GK, VectorXd& GF);
	void BuildLinearSystem_Dense(const vector<BezierElement3D>& bzmesh, MatrixXd& GK, VectorXd& GF);
	void Solver_Dense(MatrixXd& GK, VectorXd& GF);

	void FindRange(const vector<BezierElement3D>& bzmesh);
	//void setProblem_ptestHO(const TruncatedTspline_3D& tts);

	//coupling
	void Set_Problem(const vector<int>& IDBC_in);
	void Run_Coupling(const vector<BezierElement3D>& bzmesh, string fn);
	void Solver_Coupling(SparseMatrix<double>& GK, VectorXd& GF);
	void VisualizeDisp_Coupling(const vector<BezierElement3D>& bzmesh, string fn);
	void ErrorEstimate_Coupling(const vector<BezierElement3D>& bzmesh, vector<double>& err);
	void VisualizeError_Coupling(const vector<BezierElement3D>& bzmesh, const vector<double>& err, string fn);

	void Assembly_Coupling(const vector<vector<double>>& EK, const vector<double>& EF, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF);
	void BuildLinearSystem_Coupling(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF);
	void InitializeSparseMatrix(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK);
	void BasisFunction_IGA(double u, double v, double w, const BezierElement3D& bzel, vector<double>& Nx, vector<array<double, 3>>& dNdx, double& detJ);
	
	void BasisFunction_RKPM(double u, double v, double w, const BezierElement3D& bzel, vector<double>& Nx, vector<array<double, 3>>& dNdx, double& detJ, array<double,3>& x);
	void RKPM(double x[3], const BezierElement3D& bzel, vector<double>& f, vector<array<double,3>>& df);
	void Monomial(double x[3], int deg, vector<double>& ply, vector<array<double, 3>>& dp);
	void KernelFunction(double x[3], double a, double& kf, double dkf[3]);
	void KernelFunction_2(double x[3], double a, double& kf, double dkf[3]);
	void MMatrix(const vector<double>& ply, const vector<array<double, 3>>& dp, double kf, double dkf[3], MatrixXd& m, MatrixXd& dmdx, MatrixXd& dmdy, MatrixXd& dmdz);
	void ElementMatrix_RKPM(int start, const vector<array<double, 3>>& dNdx, double detJ, vector<vector<double>>& EK);
	void ElementForce_RKPM(int start, const vector<double>& Nx, double Jmod, double Fb, vector<double>& EF);

	void GetFaceInfo(double u, double v, const double pts[16][3], double& detJ, double nm[3]);
	void ElementMatrix_RKPM_Interface(int mfid, const vector<double>& Nx, const vector<array<double, 3>>& dNdx, double nm[3], double beta, double detJ, vector<vector<double>>& EK);
	void NitscheInterface(const BezierElement3D& bzel, double beta, vector<vector<double>>& EK, vector<double>& EF);
	void ElementMatrix_RKPM_Boundary(int mfid, const vector<double>& Nx, const vector<array<double, 3>>& dNdx, double nm[3], double beta, double detJ, vector<vector<double>>& EK);
	void ElementForce_RKPM_Boundary(int mfid, const vector<double>& Nx, const vector<array<double, 3>>& dNdx, double nm[3], double beta, double detJ, double Fb, vector<double>& EF);
	void NitscheBoundary_RKPM(const BezierElement3D& bzel, double beta, vector<vector<double>>& EK, vector<double>& EF);

	void ElementMatrix_IGA_Boundary(int mfid, const vector<double>& Nx, const vector<array<double, 3>>& dNdx, double nm[3], double beta, double detJ, vector<vector<double>>& EK);
	void ElementForce_IGA_Boundary(int mfid, const vector<double>& Nx, const vector<array<double, 3>>& dNdx, double nm[3], double beta, double detJ, double Fb, vector<double>& EF);
	void NitscheBoundary_IGA(const BezierElement3D& cell, double beta, vector<vector<double>>& EK, vector<double>& EF);

	void DispCal_Coupling(double u, double v, double w, const BezierElement3D& bzel, double pt[3], double& disp, double& detJ);
	void ElementError_Coupling(const BezierElement3D& bzel, double& L2, double& H1);
	void OutputError(const vector<double>& err, string fn);
	void CheckReproduce(const vector<BezierElement3D>& bzmesh, string fn);


	//coupling with Bezier
	void Set_Problem_Bezier(const vector<int>& IDBC_in, const vector<double>& gh_in);
	void Run_Coupling_Bezier(const vector<BezierElement3D>& bzmesh, string fn);
	void BasisFunction_Bezier(double u, double v, double w, const vector<array<double, 3>>& pt, vector<double>& Nx, vector<array<double, 3>>& dNdx, double& detJ);
	void BasisFunction_BZAll(double u, double v, double w, const BezierElement3D& bzel, vector<double>& Nx, vector<array<double, 3>>& dNdx, double& detJ);
	void InitializeSparseMatrix_Bezier(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK);
	void BuildLinearSystem_Coupling_Bezier(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF);
	void Assembly_Coupling_Bezier(const vector<vector<double>>& EK, const vector<double>& EF, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF);
	void NitscheBoundary_Bezier(const BezierElement3D& cell, double beta, vector<vector<double>>& EK, vector<double>& EF);
	void ElementMatrix_Bezier_Boundary(const vector<double>& Nx, const vector<array<double, 3>>& dNdx, double nm[3], double beta, double detJ, vector<vector<double>>& EK);
	void ElementForce_Bezier_Boundary(const vector<double>& Nx, const vector<array<double, 3>>& dNdx, double nm[3], double beta, double detJ, double Fb, vector<double>& EF);
	void NitscheInterface_Bezier(const BezierElement3D& cell, double beta, vector<vector<double>>& EK, vector<double>& EF);
	void ElementMatrix_Bezier_Interface(int mfid, const vector<double>& Nx, const vector<array<double, 3>>& dNdx, double nm[3], double beta, double detJ, vector<vector<double>>& EK);
	void VisualizeDisp_Coupling_Bezier(const vector<BezierElement3D>& bzmesh, string fn);
	void DispCal_Coupling_Bezier(double u, double v, double w, const BezierElement3D& bzel, double pt[3], double& disp, double dudx[3], double& detJ);
	void ErrorEstimate_Coupling_Bezier(const vector<BezierElement3D>& bzmesh, vector<vector<double>>& err);
	void ElementError_Coupling_Bezier(const BezierElement3D& bzel, double& L2, double& H1);
	void OutputError_All(const vector<vector<double>>& err, string fn);
	void Solver_Coupling_Bezier(SparseMatrix<double>& GK, VectorXd& GF);
	void Solver_Matlab(SparseMatrix<double>& GK, VectorXd& GF);
	void Solver_MKL(SparseMatrix<double>& GK, VectorXd& GF);
	void Write2MKL_SparseMat(SparseMatrix<double>& GK, VectorXd& GF, string fn);
	void Write2Matlab_SparseMat(SparseMatrix<double>& GK, VectorXd& GF, string fn);
	void ReadSolution(string fn);
	void VisualizeError_BezierCouple(const vector<BezierElement3D>& bzmesh, const vector<double>& err, string fn);

	//for comparison
	void NitscheBoundary_Spline(const BezierElement3D& cell, double beta, vector<vector<double>>& EK, vector<double>& EF);


	double f_source_lap(const array<double, 3>& pt);
	//double f_source_1(const array<double,3>& pt);
	//double f_source_2(const array<double, 3>& pt);
	//double f_source_3(const array<double, 3>& pt);
	//double f_source_4(const array<double, 3>& pt);
	//double f_source_5(const array<double, 3>& pt);
	//double f_source_6(const array<double, 3>& pt);
	double f_source_7_lap(const array<double, 3>& pt);
	//double f_source_8(const array<double, 3>& pt);
	//double f_source_9(const array<double, 3>& pt);
	double exact_sol_lap(const array<double, 3>& x);
	//double exact_sol_1(const array<double, 3>& x);
	//double exact_sol_2(const array<double, 3>& x);
	//double exact_sol_3(const array<double, 3>& x);
	//double exact_sol_4(const array<double, 3>& x);
	//double exact_sol_5(const array<double, 3>& x);
	//double exact_sol_6(const array<double, 3>& x);
	double exact_sol_7_lap(const array<double, 3>& x);
	//double exact_sol_8(const array<double, 3>& x);
	//double exact_sol_9(const array<double, 3>& x);

	void GetRemoveRegion(double xy[3][2]);
	void GetEqParameter(double xy[3][2], double nrm[3], double a);

	void PipelineTmp(const vector<BezierElement3D>& bzmesh, string fn);
	
private:
	vector<double> Gpt;
	vector<double> wght;
	vector<int> IDBC;
	vector<int> BCList;
	vector<double> gh;
	vector<double> uh;
	VectorXd u0;
	long long int npt;//not enough, INT_MAX=32767
	long long int neq;
	double range[3][2];
	double rmv[3][2];

	int nglev;
	double hemax;

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