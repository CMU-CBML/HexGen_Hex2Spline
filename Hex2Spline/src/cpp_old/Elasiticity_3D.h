#ifndef ELASITICITY_3D_H
#define ELASITICITY_3D_H

#include <vector>
#include <array>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "BasicDataStructure.h"
//#include "Truncated_Tspline.h"
#include "mkl_pardiso.h"
#include "mkl_types.h"

using namespace std;
using namespace Eigen;

class LinearElasticity
{
public:
	LinearElasticity();
	void SetMaterialProp(double Ey=1., double v=0.3);
	void GaussInfo(int ng=3);
	void SetProblem(const vector<int>& IDBC_in, const vector<double>& gh_in);
	void SetBoundary();
	void BasisFunction(double u, double v, double w, const vector<array<double, 3>>& pt, double Nx[64], double dNdx[64][3], double& detJ);
	void ElementMatrix(double dNdx[64][3], double detJ, double EK[192][192]);
	void ElementForce(double Nx[16], double Jmod, double Fb[2], double EF[32]);
	void Assembly(double EK[192][192], double EF[192], const vector<vector<double>>& cmat, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF);
	void BuildLinearSystem(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF);
	void Solver(SparseMatrix<double>& GK, VectorXd& GF);
	void Para2Phys(double u, double v, const double cpt[16][3], double pt[3]);
	void DispStrainStress(double u,double v,double w, const BezierElement3D& bzel,double disp[3],double se[6],double ss[6]);
	//void StressCalculate(double u,double v,const BezierElement& bzel,double stress[6]);
	void VisualizeVTK(const vector<BezierElement3D>& bzmesh, string fn);
	//void VisualizeError(string fn);
	//void ElementErrorEstimate(int eid,double& err,double& area,vector<int>& pid);
	//void TotalErrorEstimate();
	void Run(const vector<BezierElement3D>& bzmesh, string fn);

	void OutputMatrix(SparseMatrix<double>& GK, string fn);

	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void Set_Problem_Bezier(const vector<int>& IDBC_in, const vector<double>& gh_in);
	void Run_Coupling_Bezier(const vector<BezierElement3D>& bzmesh, string fn);
	void InitializeSparseMatrix_Bezier(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK);
	void BuildLinearSystem_Coupling_Bezier(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF);
	void BasisFunction_IGA(double u, double v, double w, const BezierElement3D& bzel, vector<double>& Nx, vector<array<double, 3>>& dNdx, double& detJ);
	void ElementMatrix(const vector<array<double, 3>>& dNdx, double detJ, vector<vector<double>>& EK);
	void Assembly_Coupling_Bezier(const vector<vector<double>>& EK, const vector<double>& EF, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF);
	void Solver_Coupling_Bezier(SparseMatrix<double>& GK, VectorXd& GF);
	void Solver_MKL(SparseMatrix<double>& GK, VectorXd& GF);
	void VisualizeDisp_Coupling_Bezier(const vector<BezierElement3D>& bzmesh, string fn);
	void DispCal_Coupling_Bezier(const vector<int>& IEN, const vector<double>& Nx, array<double, 3>& disp);
	void StressCal_Coupling_Bezier(const vector<int>& IEN, const vector<array<double,3>>& dNdx, array<double, 7>& ss);
	double MaxVonMises(const vector<BezierElement3D>& bzmesh, vector<double>& vs);
	void ReadSolutionAbaqus(string fn);
	void WriteDispAbaqus(string fn);
	void ErrorEvaluation(const vector<BezierElement3D>& bzmesh);//patch test

	//sparse matrix using MKL
	void Run_MKL(const vector<BezierElement3D>& bzmesh, string fn);
	void InitializeSparseMatrix_MKL(const vector<BezierElement3D>& bzmesh, 
		vector<MKL_INT>& ja, vector<MKL_INT>& ia, vector<double>& a);
	void BuildLinearSystem_MKL(const vector<BezierElement3D>& bzmesh, 
		vector<MKL_INT>& ja, vector<MKL_INT>& ia, vector<double>& a, vector<double>& b);
	void Assembly_MKL(const vector<vector<double>>& EK, const vector<double>& EF, const vector<int>& IEN,
		vector<MKL_INT>& Acol, vector<MKL_INT>& Arow, vector<double>& Aval, vector<double>& GF);
	void Solver_MKL_Direct(vector<MKL_INT>& ja, vector<MKL_INT>& ia, vector<double>& a, vector<double>& b);

	void Run_UB_AdaptFit(const vector<BezierElement3D>& bzmesh, string fn, vector<double>& vs, double& vsmax);
	
private:
	vector<double> Gpt;
	vector<double> wght;
	double lambda;//Lame parameters
	double mu;
	vector<int> IDBC;
	vector<int> BCList;
	vector<int> BC_Neumann;
	vector<double> gh;
	vector<double> uh;
	int npt;
	int neq;
	int dim;

	vector<array<double, 3>> pts;
};

//void ExactDisp_PlateHole(double x, double y, double u[2]);
//void ExactStress_PlateHole(double x, double y, double ss[3]);

#endif