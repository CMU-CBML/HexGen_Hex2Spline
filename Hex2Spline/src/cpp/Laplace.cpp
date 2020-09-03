#include "Laplace.h"
#include "PDEdata.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <ctime>
#include <iomanip>
//#include "Matlab_Solver_wap.h"
//#include <Eigen/IterativeLinearSolvers>
//#include <unsupported/Eigen/IterativeSolvers>
#include "mkl_pardiso.h"
#include "mkl_types.h"

#define EIGEN_DONT_PARALLELIZE

using namespace std;

typedef unsigned int uint;
typedef long long int llint;

Laplace::Laplace()
{
	range[0][0] = 0.; range[0][1] = 1.;
	range[1][0] = 0.; range[1][1] = 1.;
	range[2][0] = 0.; range[2][1] = 1.;

	dmlen[0] = range[0][1] - range[0][0];
	dmlen[1] = range[1][1] - range[1][0];
	dmlen[2] = range[2][1] - range[2][0];
	nmpl[0] = -0.408248; nmpl[1] = -0.408248; nmpl[2] = 0.816497;
	acoef = 50.;
	xorg[0] = 0.; xorg[1] = 0.; xorg[2] = 0.;

	nglev = 0;
}

void Laplace::GaussInfo(int ng)
{
	Gpt.clear();
	wght.clear();
	switch(ng)
	{
	case 2:
		{
			Gpt.resize(ng);
			wght.resize(ng);
			Gpt[0]=0.2113248654051871;			Gpt[1]=0.7886751345948129;
			wght[0]=1.;			wght[1]=1.;
			break;
		}
	case 3:
		{
			Gpt.resize(ng);
			wght.resize(ng);
			Gpt[0]=0.1127016653792583;			Gpt[1]=0.5;			Gpt[2]=0.8872983346207417;
			wght[0]=0.5555555555555556;			wght[1]=0.8888888888888889;			wght[2]=0.5555555555555556;
			break;
		}
	case 4:
		{
			Gpt.resize(ng);
			wght.resize(ng);
			Gpt[0]=0.06943184420297371;			Gpt[1]=0.33000947820757187;			Gpt[2]=0.6699905217924281;			Gpt[3]=0.9305681557970262;
			wght[0]=0.3478548451374539;			wght[1]=0.6521451548625461;			wght[2]=0.6521451548625461;			wght[3]=0.3478548451374539;
			break;
		}
	case 5:
	{
			  Gpt.resize(ng);
			  wght.resize(ng);
			  Gpt[0] = 0.046910077030668;			Gpt[1] = 0.2307653449471585;			Gpt[2] = 0.5;			Gpt[3] = 0.7692346550528415;  Gpt[4] = 0.953089922969332;
			  wght[0] = 0.2369268850561891;			wght[1] = 0.4786286704993665;			wght[2] = 0.5688888888888889;			wght[3] = 0.4786286704993665; wght[4] = 0.2369268850561891;
			  break;
	}
	default:
		{
			Gpt.resize(2);
			wght.resize(2);
			Gpt[0]=0.2113248654051871;			Gpt[1]=0.7886751345948129;
			wght[0]=1.;			wght[1]=1.;
			break;
		}
	}
}

void Laplace::SetProblem(const vector<int>& IDBC_in, const vector<double>& gh_in)
{
	npt = IDBC_in.size();
	IDBC = IDBC_in;
	gh = gh_in;
	neq = 0;
	//int pid(0);
	for (unsigned int i = 0; i < IDBC.size(); i++)
	{
		if (IDBC[i] != -1) neq++;
	}

	//uh.resize(npt);
	//for (int i = 0; i<npt; i++)
	//{
	//	if (IDBC[i] != -1)
	//		uh[i] = 0.;
	//	else
	//		uh[i] = gh[i];
	//}
}

void Laplace::BasisFunction(double u, double v, double w, const vector<array<double, 3>>& pt, double Nx[64], double dNdx[64][3], double& detJ)
{
	double Nu[4] = { (1. - u)*(1. - u)*(1. - u), 3.*(1. - u)*(1. - u)*u, 3.*(1. - u)*u*u, u*u*u };
	double Nv[4] = { (1. - v)*(1. - v)*(1. - v), 3.*(1. - v)*(1. - v)*v, 3.*(1. - v)*v*v, v*v*v };
	double Nw[4] = { (1. - w)*(1. - w)*(1. - w), 3.*(1. - w)*(1. - w)*w, 3.*(1. - w)*w*w, w*w*w };
	double dNdu[4] = { -3.*(1. - u)*(1. - u), 3. - 12.*u + 9.*u*u, 3.*(2. - 3.*u)*u, 3.*u*u };
	double dNdv[4] = { -3.*(1. - v)*(1. - v), 3. - 12.*v + 9.*v*v, 3.*(2. - 3.*v)*v, 3.*v*v };
	double dNdw[4] = { -3.*(1. - w)*(1. - w), 3. - 12.*w + 9.*w*w, 3.*(2. - 3.*w)*w, 3.*w*w };
	double dNdt[64][3];
	int i, j, k, a, b, loc(0);
	for (i = 0; i<4; i++)
	{
		for (j = 0; j<4; j++)
		{
			for (k = 0; k < 4; k++)
			{
				Nx[loc] = Nu[k] * Nv[j] * Nw[i];
				dNdt[loc][0] = dNdu[k] * Nv[j] * Nw[i];
				dNdt[loc][1] = Nu[k] * dNdv[j] * Nw[i];
				dNdt[loc][2] = Nu[k] * Nv[j] * dNdw[i];
				loc++;
			}
		}
	}
	Matrix3d dxdt = Matrix3d::Zero();
	loc = 0;
	for (i = 0; i<4; i++)
	{
		for (j = 0; j<4; j++)
		{
			for (k = 0; k < 4; k++)
			{
				for (a = 0; a<3; a++)
				{
					for (b = 0; b<3; b++)
					{
						dxdt(a, b) += pt[loc][a] * dNdt[loc][b];
					}
				}
				loc++;
			}
		}
	}
	//cout << dxdt << "\n\n";
	Matrix3d dtdx = dxdt.inverse();
	for (i = 0; i<64; i++)
	{
		dNdx[i][0] = dNdt[i][0] * dtdx(0, 0) + dNdt[i][1] * dtdx(1, 0) + dNdt[i][2] * dtdx(2, 0);
		dNdx[i][1] = dNdt[i][0] * dtdx(0, 1) + dNdt[i][1] * dtdx(1, 1) + dNdt[i][2] * dtdx(2, 1);
		dNdx[i][2] = dNdt[i][0] * dtdx(0, 2) + dNdt[i][1] * dtdx(1, 2) + dNdt[i][2] * dtdx(2, 2);
	}
	detJ = dxdt.determinant();
	detJ = 0.125*detJ;

	//double tmp = detJ;
	////cout << "detJ: " << detJ << "\n";
	//detJ = dxdt(0, 0)*(dxdt(1, 1)*dxdt(2, 2) - dxdt(1, 2)*dxdt(2, 1)) - dxdt(0, 1)*(dxdt(1, 0)*dxdt(2, 2) - dxdt(1, 2)*dxdt(2, 0))
	//	+ dxdt(0, 2)*(dxdt(1, 0)*dxdt(2, 1) - dxdt(1, 1)*dxdt(2, 0));
	//Matrix3d dtx = Matrix3d::Zero();
	//dtx << (dxdt(1, 1)*dxdt(2, 2) - dxdt(1, 2)*dxdt(2, 1)) / detJ, (dxdt(0, 2)*dxdt(2, 1) - dxdt(0, 1)*dxdt(2, 2)) / detJ, (dxdt(0, 1)*dxdt(1, 2) - dxdt(0, 2)*dxdt(1, 1)) / detJ,
	//	(dxdt(1, 2)*dxdt(2, 0) - dxdt(1, 0)*dxdt(2, 2)) / detJ, (dxdt(0, 0)*dxdt(2, 2) - dxdt(0, 2)*dxdt(2, 0)) / detJ, (dxdt(0, 2)*dxdt(1, 0) - dxdt(0, 0)*dxdt(1, 2)) / detJ,
	//	(dxdt(1, 0)*dxdt(2, 1) - dxdt(1, 1)*dxdt(2, 0)) / detJ, (dxdt(0, 1)*dxdt(2, 0) - dxdt(0, 0)*dxdt(2, 1)) / detJ, (dxdt(0, 0)*dxdt(1, 1) - dxdt(0, 1)*dxdt(1, 0)) / detJ;
	//Matrix3d mt = dtx - dtdx;
	//cout << mt << "\n";
	////detJ *= 0.125;
	////cout << "detJ diff: " << fabs(detJ-tmp) << "\n";
	//getchar();
}

void Laplace::ElementMatrix(double dNdx[64][3], double detJ, double EK[64][64])
{
	int i, j;
//#pragma omp parallel for
	for (i = 0; i<64; i++)
	{
		for (j = 0; j<64; j++)
		{
			EK[i][j] += (dNdx[i][0] * dNdx[j][0] + dNdx[i][1] * dNdx[j][1] + dNdx[i][2] * dNdx[j][2])*detJ;
		}
	}
}

void Laplace::ElementForce(double Nx[64], double detJ, double Fb, double EF[64])
{
	for (int i = 0; i<64; i++)
	{
		EF[i] += Nx[i] * detJ*Fb;
	}
}

void Laplace::Assembly(double EK1[64][64], double EF1[64], const vector<vector<double>>& cmat, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF)
{
	//unsigned int i, j, A, B;
	int i, j, A, B;
	vector<vector<double>> EK(IEN.size(), vector<double>(IEN.size(), 0.));
	vector<double> EF(IEN.size(), 0.);
	//cout << "before mapping\n";
//#pragma omp parallel for
	for (i = 0; i<IEN.size(); i++)
	{
		for (j = 0; j<IEN.size(); j++)
		{
			for (A = 0; A<64; A++)
			{
				for (B = 0; B<64; B++)
				{
					EK[i][j] += cmat[i][A] * cmat[j][B] * EK1[A][B];
				}
			}
		}
	}
//#pragma omp barrier
	//cout << "after mapping\n";
//#pragma omp parallel for
	for (i = 0; i < IEN.size(); i++)
	{
		for (A = 0; A < 64; A++)
		{
			EF[i] += cmat[i][A] * EF1[A];
		}
	}
//#pragma omp barrier
	//cout << "before global\n";
	for (i = 0; i<IEN.size(); i++)
	{
		for (j = 0; j<IEN.size(); j++)
		{
			A = IEN[i]; B = IEN[j];
			if (A >= IDBC.size() || B >= IDBC.size())
			{
				cout << "out of range!\n";
				cout << "A: " << A << "/" << IDBC.size() << "\n";
				cout << "B: " << B << "/" << IDBC.size() << "\n";
				getchar();
			}
			if (IDBC[A] != -1 && IDBC[B] != -1)
			{
				GK.coeffRef(IDBC[A], IDBC[B]) += EK[i][j];
			}
			else if (IDBC[A] != -1 && IDBC[B] == -1)
			{
				GF(IDBC[A]) -= EK[i][j] * gh[B];
			}
		}
	}
	for (i = 0; i < IEN.size(); i++)
	{
		A = IEN[i];
		if (IDBC[A] != -1)
		{
			GF(IDBC[A]) += EF[i];
		}
	}
	//cout << "after global\n";
}

void Laplace::BuildLinearSystem(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF)
{
	int e, i, j, k, a, b;
	double EK[64][64];
	double EF[64];
	double Nx[64];
	double dNdx[64][3];
	double detJ, Fb(0.);
	array<double,3> pt;
	cout << "npt: " << npt << "\n";
	cout << "neq: " << neq << "\n";
	cout << "# of eles: " << bzmesh.size() << "\n";
	for (e = 0; e<bzmesh.size(); e++)
	{
		//cout << "element id: " << e << "/" << bzmesh.size() << "\n";
		if (e != 0 && e % 100 == 0)
		{
			cout << e << " ";
		}
//#pragma omp parallel for
		for (i = 0; i<64; i++)
		{
			EF[i] = 0.;
			for (j = 0; j<64; j++)
			{
				EK[i][j] = 0.;
			}
		}
//#pragma omp barrier
		//cout << "before loop\n";
		for (i = 0; i<Gpt.size(); i++)
		{
			for (j = 0; j<Gpt.size(); j++)
			{
				for (k = 0; k < Gpt.size(); k++)
				{
					BasisFunction(Gpt[i], Gpt[j], Gpt[k], bzmesh[e].pts, Nx, dNdx, detJ);
					detJ = wght[i] * wght[j] * wght[k] * detJ;
					ElementMatrix(dNdx, detJ, EK);

					Para2Phys(Gpt[i], Gpt[j], Gpt[k], bzmesh[e].pts, Nx, pt);
					Fb = f_source(pt);
					ElementForce(Nx, detJ, Fb, EF);
				}
			}
		}
		//cout << "after loop\n";
		Assembly(EK, EF, bzmesh[e].cmat, bzmesh[e].IEN, GK, GF);
		//cout << "after assembly\n";
	}

	//cout << GF << "\n";
	//getchar();
}

void Laplace::Solver(SparseMatrix<double>& GK, VectorXd& GF)
{
	cout << "Solving linear system...\n";

	SimplicialLDLT<SparseMatrix<double>> solver;
	VectorXd sol = solver.compute(GK).solve(GF);

	//ConjugateGradient<SparseMatrix<double>, Lower|Upper> solver;
	////double tol(double(neq)*10.e-16);
	////double tol(1.e-6);
	////cout << "tol before: " << solver.tolerance() << "\n";
	////solver.setTolerance(tol);
	////cout << "tol after: " << solver.tolerance() << "\n";
	////solver.setMaxIterations(300);
	//solver.compute(GK);
	//VectorXd x0 = VectorXd::Ones(neq);
	//VectorXd sol = solver.solve(GF);
	//VectorXd sol = solver.solveWithGuess(GF, x0);
	//cout << "\n# iterations: " << solver.iterations() << '\n';
	//cout << "estimated error: " << solver.error() << '\n';
	//getchar();

	uh.resize(npt);
	for(int i=0; i<npt; i++)
	{
		if (IDBC[i] != -1)
			uh[i]=sol(IDBC[i]);
		else
			uh[i]=gh[i];
	}

	cout << "Done solving!\n";

	//cout << GF << "\n";
	//cout << sol << "\n";
	//cout << GF.maxCoeff() << "\n";
	//cout << GF.minCoeff() << "\n";
	//getchar();

	//uh.resize(npt);
	//for (int i = 0; i<npt; i++)
	//{
	//	if (IDBC[i] != -1)
	//		uh[i] = 0.;
	//	else
	//		uh[i] = gh[i];
	//}
}

void Laplace::Solver_Dense(MatrixXd& GK, VectorXd& GF)
{
	//SimplicialLDLT<SparseMatrix<double>> solver;
	VectorXd sol = GK.ldlt().solve(GF);

	//ConjugateGradient<SparseMatrix<double>, Lower|Upper> solver;
	//solver.compute(GK);
	//VectorXd sol = solver.solve(GF);
	//cout << "\n# iterations: " << solver.iterations() << '\n';
	//cout << "estimated error: " << solver.error() << '\n';
	//getchar();

	uh.resize(npt);
	for (int i = 0; i<npt; i++)
	{
		if (IDBC[i] != -1)
			uh[i] = sol(IDBC[i]);
		else
			uh[i] = gh[i];
	}

	//uh.resize(npt);
	//for (int i = 0; i<npt; i++)
	//{
	//	if (IDBC[i] != -1)
	//		uh[i] = 0.;
	//	else
	//		uh[i] = gh[i];
	//}
}

void Laplace::Para2Phys(double u, double v, double w, const vector<array<double, 3>>& bzpt, double Nx[64], array<double, 3>& pt)
{
	pt[0] = 0.; pt[1] = 0.; pt[2] = 0.;
	for (int i = 0; i < 64; i++)
	{
		pt[0] += bzpt[i][0] * Nx[i];
		pt[1] += bzpt[i][1] * Nx[i];
		pt[2] += bzpt[i][2] * Nx[i];
	}
}

void Laplace::DispCal(double u,double v,double w,const BezierElement3D& bzel,double& disp,double& detJ)
{
	double Nx[64];
	double dNdx[64][3];
	BasisFunction(u,v,w,bzel.pts,Nx,dNdx,detJ);
	double uloc[64];
	unsigned int i,j;
	for(i=0;i<64;i++)
	{
		uloc[i]=0.;
		for(j=0;j<bzel.IEN.size();j++)
		{
			uloc[i]+=bzel.cmat[j][i]*uh[bzel.IEN[j]];
		}
	}
	//displacement
	disp=0.;
	for(i=0;i<64;i++)
	{
		disp+=Nx[i]*uloc[i];
	}
}

void Laplace::ElementError(const BezierElement3D& bzel, double& L2, double& H1)
{
	L2=0.; H1=0.;
	uint i, j, k;
	double pt1[3], disp, detJ, sol_e;
	array<double, 3> pt;
	for (i = 0; i < Gpt.size(); i++)
	{
		for (j = 0; j < Gpt.size(); j++)
		{
			for (k = 0; k < Gpt.size(); k++)
			{
				//BasisFunction(Gpt[i], Gpt[j], Gpt[k], bzel.pts, Nx, dNdx, detJ);
				bzel.Para2Phys(Gpt[i], Gpt[j], Gpt[k], pt1);
				pt[0] = pt1[0]; pt[1] = pt1[1]; pt[2] = pt1[2];
				DispCal(Gpt[i], Gpt[j], Gpt[k], bzel, disp, detJ);
				sol_e = exact_sol_lap(pt);
				//sol_e = exact_sol(pt);
				L2 += (wght[i] * wght[j] * wght[k] * detJ*(disp - sol_e)*(disp - sol_e));
			}
		}
	}

	//for(unsigned int gid1=0;gid1<Gpt.size();gid1++)
	//{
	//	for(unsigned int gid2=0;gid2<Gpt.size();gid2++)
	//	{
	//		double us, ue, detJ, val[3], ux, uy;
	//		array<double,3> x;
	//		//DispCal(Gpt[gid1],Gpt[gid2],bzel,us,detJ);
	//		//Para2Phys(Gpt[gid1],Gpt[gid2],bzel.pts,x);
	//		//DispCal_TSP(Gpt[gid1],Gpt[gid2],bzel,us,detJ);
	//		Quantityh_TSP(Gpt[gid1],Gpt[gid2],bzel,val,detJ);
	//		Para2Phys_TSP(Gpt[gid1],Gpt[gid2],bzel,x);
	//		//ue=LDomainSolution(x[0],x[1]);
	//		//LDomainSolutionDeriv(x[0],x[1],ue,ux,uy);
	//		//ue=ptestHO_sol_1(x[0],x[1]);
	//		//ue=ptestHO_sol_2(x[0],x[1]);
	//		ue=ptestHO_sol_3(x[0],x[1]);
	//		//L2+=wght[gid1]*wght[gid2]*detJ*(us-ue)*(us-ue);
	//		L2+=wght[gid1]*wght[gid2]*detJ*(val[0]-ue)*(val[0]-ue);
	//		//H1+=wght[gid1]*wght[gid2]*detJ*((val[1]-ux)*(val[1]-ux)+(val[2]-uy)*(val[2]-uy));
	//	}
	//}
}

void Laplace::ErrorEstimate(const vector<BezierElement3D>& bzmesh, vector<double>& err)
{
	err.clear();
	err.resize(bzmesh.size(),0.);
#pragma omp parallel for
	for (int e = 0; e < bzmesh.size(); e++)
	{
		double L2, H1;
		ElementError(bzmesh[e],L2,H1);
		err[e] = L2;
	}
}

double Laplace::ElementErrorEstimate(const BezierElement3D& bzel)
{
	double err_energy(0.), left(0.), right(0.);
	double detJ, dtdx[3][3], pt[3];
	vector<double> Nx;
	vector<array<double, 3>> dNdx;
	double w, dwdx[3], dudx[3], fb(0.);
	array<double, 3> pt1;
	uint i, j, k;
	for (i = 0; i < Gpt.size(); i++)
	{
		for (j = 0; j < Gpt.size(); j++)
		{
			for (k = 0; k < Gpt.size(); k++)
			{
				BasisFunctionMore(Gpt[i], Gpt[j], Gpt[k], bzel, pt, Nx, dNdx, detJ, dtdx);
				BubbleFunction(Gpt[i], Gpt[j], Gpt[k], dtdx, w, dwdx);
				DispFunction(bzel.IEN, dNdx, dudx);
				detJ = wght[i] * wght[j] * wght[k] * detJ;
				left += detJ*(dwdx[0] * dwdx[0] + dwdx[1] * dwdx[1] + dwdx[2] * dwdx[2]);
				right -= detJ*(dudx[0] * dwdx[0] + dudx[1] * dwdx[1] + dudx[2] * dwdx[2]);
				//pt1[0] = pt[0]; pt1[1] = pt[1]; pt1[2] = pt[2];
				//fb=f_source_3(pt1);
				//right += detJ*fb*w;
			}
		}
	}
	if (left != 0.)
	{
		double coef=right/left;
		err_energy = coef*coef*left;
	}
	else
	{
		cerr << "wrong left!\n";
		getchar();
	}
	return err_energy;
}

void Laplace::BasisFunctionMore(double u, double v, double w, const BezierElement3D& bzel, double pt[3], vector<double>& Nx, vector<array<double, 3>>& dNdx, double& detJ, double dtdx1[3][3])
{
	double Nu[4] = { (1. - u)*(1. - u)*(1. - u), 3.*(1. - u)*(1. - u)*u, 3.*(1. - u)*u*u, u*u*u };
	double Nv[4] = { (1. - v)*(1. - v)*(1. - v), 3.*(1. - v)*(1. - v)*v, 3.*(1. - v)*v*v, v*v*v };
	double Nw[4] = { (1. - w)*(1. - w)*(1. - w), 3.*(1. - w)*(1. - w)*w, 3.*(1. - w)*w*w, w*w*w };
	double dNdu[4] = { -3.*(1. - u)*(1. - u), 3. - 12.*u + 9.*u*u, 3.*(2. - 3.*u)*u, 3.*u*u };
	double dNdv[4] = { -3.*(1. - v)*(1. - v), 3. - 12.*v + 9.*v*v, 3.*(2. - 3.*v)*v, 3.*v*v };
	double dNdw[4] = { -3.*(1. - w)*(1. - w), 3. - 12.*w + 9.*w*w, 3.*(2. - 3.*w)*w, 3.*w*w };
	double Nt[64];
	double dNdt[64][3], dNdx0[64][3];
	int i, j, k, a, b, loc(0);
	pt[0] = 0.; pt[1] = 0.; pt[2] = 0.;
	for (i = 0; i<4; i++)
	{
		for (j = 0; j<4; j++)
		{
			for (k = 0; k < 4; k++)
			{
				Nt[loc] = Nu[k] * Nv[j] * Nw[i];
				pt[0] += bzel.pts[loc][0] * Nt[loc];
				pt[1] += bzel.pts[loc][1] * Nt[loc];
				pt[2] += bzel.pts[loc][2] * Nt[loc];
				dNdt[loc][0] = dNdu[k] * Nv[j] * Nw[i];
				dNdt[loc][1] = Nu[k] * dNdv[j] * Nw[i];
				dNdt[loc][2] = Nu[k] * Nv[j] * dNdw[i];
				loc++;
			}
		}
	}
	Matrix3d dxdt = Matrix3d::Zero();
	loc = 0;
	for (i = 0; i<4; i++)
	{
		for (j = 0; j<4; j++)
		{
			for (k = 0; k < 4; k++)
			{
				for (a = 0; a<3; a++)
				{
					for (b = 0; b<3; b++)
					{
						dxdt(a, b) += bzel.pts[loc][a] * dNdt[loc][b];
					}
				}
				loc++;
			}
		}
	}
	//cout << dxdt << "\n\n";
	Matrix3d dtdx = dxdt.inverse();
	for (i = 0; i<64; i++)
	{
		dNdx0[i][0] = dNdt[i][0] * dtdx(0, 0) + dNdt[i][1] * dtdx(1, 0) + dNdt[i][2] * dtdx(2, 0);
		dNdx0[i][1] = dNdt[i][0] * dtdx(0, 1) + dNdt[i][1] * dtdx(1, 1) + dNdt[i][2] * dtdx(2, 1);
		dNdx0[i][2] = dNdt[i][0] * dtdx(0, 2) + dNdt[i][1] * dtdx(1, 2) + dNdt[i][2] * dtdx(2, 2);
	}
	detJ = dxdt.determinant();
	detJ = 0.125*detJ;
	Nx.clear(); dNdx.clear();
	Nx.resize(bzel.IEN.size(),0.); dNdx.resize(bzel.IEN.size());
	for (i = 0; i < bzel.cmat.size(); i++)
	{
		dNdx[i][0] = 0.; dNdx[i][1] = 0.; dNdx[i][2] = 0.;
		for (j = 0; j < bzel.cmat[i].size(); j++)
		{
			Nx[i] += bzel.cmat[i][j] * Nt[j];
			dNdx[i][0] += bzel.cmat[i][j] * dNdx0[j][0];
			dNdx[i][1] += bzel.cmat[i][j] * dNdx0[j][1];
			dNdx[i][2] += bzel.cmat[i][j] * dNdx0[j][2];
		}
	}
	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			dtdx1[i][j] = dtdx(i,j);
		}
	}
}

void Laplace::BubbleFunction(double u, double v, double w, double dtdx[3][3], double& f, double dfdx[3])
{
	double wt[3] = { u*(1. - u), v*(1. - v), w*(1. - w)};
	double dwdt[3] = { 1. - 2.*u, 1. - 2.*v, 1. - 2.*w};
	f = wt[0]*wt[1]*wt[2];
	double dfdt[3] = { dwdt[0] * wt[1] * wt[2], wt[0] * dwdt[1] * wt[2], wt[0] * wt[1] * dwdt[2] };
	for (int i = 0; i < 3; i++)
	{
		dfdx[i] = dfdt[0] * dtdx[0][i] + dfdt[1] * dtdx[1][i] + dfdt[2] * dtdx[2][i];
	}
}

void Laplace::DispFunction(const vector<int>& IEN, const vector<array<double, 3>>& dNdx, double dfdx[3])
{
	dfdx[0] = 0.; dfdx[1] = 0.; dfdx[2] = 0.;
	for (uint i = 0; i < IEN.size(); i++)
	{
		dfdx[0] += uh[IEN[i]] * dNdx[i][0];
		dfdx[1] += uh[IEN[i]] * dNdx[i][1];
		dfdx[2] += uh[IEN[i]] * dNdx[i][2];
	}
}

void Laplace::TotalErrorEstimate(const vector<BezierElement3D>& bzmesh, vector<double>& err)
{
	double enorm;
	err.clear();
	err.resize(bzmesh.size(), 0.);
	for (uint e = 0; e < bzmesh.size(); e++)
	{
		enorm=ElementErrorEstimate(bzmesh[e]);
		err[e] = enorm;
	}
}

void Laplace::VisualizeVTK(const vector<BezierElement3D>& bzmesh, string fn)
{
	vector<array<double, 3>> spt;//s means sample
	vector<double> sdisp;
	//vector<array<double, 3>> sdisp_err;
	//vector<array<double, 3>> sse;
	//vector<array<double, 3>> sss;
	vector<array<int, 8>> sele;
	vector<double> errL2;
	vector<array<double, 3>> lpt;//visulize parameter lines
	vector<array<int, 2>> led;//line connectivity
	int ns(4), ecount(0);
	vector<double> su(ns);
	for (int i = 0; i<ns; i++)
	{
		su[i] = double(i) / (double(ns) - 1.);
	}

	for (unsigned int e = 0; e<bzmesh.size(); e++)
	{
		//double errtmp;
		//ElementError(bzmesh[e],errtmp);
		//errL2.push_back(sqrt(errtmp));
		int loc(0);
		int pstart = spt.size();
		for (int a = 0; a<ns; a++)
		{
			for (int b = 0; b<ns; b++)
			{
				for (int c = 0; c < ns; c++)
				{
					double pt1[3];
					double disp, detmp;
					bzmesh[e].Para2Phys(su[c], su[b], su[a], pt1);
					//DispCal(su[c], su[b], su[a], bzmesh[e], disp, detmp);
					array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
					spt.push_back(pt);
					//disp = exact_sol(pt);
					//sdisp.push_back(disp);
					if (a == 0 || a == ns - 1 || b == 0 || b == ns - 1)
					{
						lpt.push_back(pt);
					}
				}
			}
		}
		int nns[2] = { ns*ns*ns, ns*ns };
		for (int a = 0; a < ns - 1; a++)
		{
			for (int b = 0; b < ns - 1; b++)
			{
				for (int c = 0; c < ns - 1; c++)
				{
					array<int, 8> el;
					el[0] = pstart + a * nns[1] + b * ns + c;
					el[1] = pstart + a * nns[1] + b * ns + c + 1;
					el[2] = pstart + a * nns[1] + (b + 1)*ns + c + 1;
					el[3] = pstart + a * nns[1] + (b + 1)*ns + c;
					el[4] = pstart + (a + 1)*nns[1] + b * ns + c;
					el[5] = pstart + (a + 1)*nns[1] + b * ns + c + 1;
					el[6] = pstart + (a + 1)*nns[1] + (b + 1)*ns + c + 1;
					el[7] = pstart + (a + 1)*nns[1] + (b + 1)*ns + c;
					sele.push_back(el);
				}
			}
		}
		int lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(su[a], 0., 0., pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(su[a], 1., 0., pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(su[a], 0., 1., pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(su[a], 1., 1., pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(0., su[a], 0., pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(1., su[a], 0., pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(0., su[a], 1., pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(1., su[a], 1., pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(0., 0., su[a], pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(1., 0., su[a], pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(0., 1., su[a], pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(1., 1., su[a], pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		/*for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> lc;
			lc[0] = ecount * 4 * (ns - 1) + a;
			lc[1] = ecount * 4 * (ns - 1) + a + 1;
			led.push_back(lc);
			lc[0] = ecount * 4 * (ns - 1) + 3 * ns - 4 + a;
			lc[1] = ecount * 4 * (ns - 1) + 3 * ns - 4 + a + 1;
			led.push_back(lc);
		}
		for (int a = 0; a < ns - 2; a++)
		{
			array<int, 2> lc;
			lc[0] = ecount * 4 * (ns - 1) + ns + 2 * a;
			lc[1] = ecount * 4 * (ns - 1) + ns + 2 * a + 2;
			led.push_back(lc);
			lc[0] = ecount * 4 * (ns - 1) + ns + 2 * a - 1;
			lc[1] = ecount * 4 * (ns - 1) + ns + 2 * a + 1;
			led.push_back(lc);
		}
		array<int, 2> lc1;
		lc1[0] = ecount * 4 * (ns - 1);
		lc1[1] = ecount * 4 * (ns - 1) + ns;
		led.push_back(lc1);
		lc1[0] = ecount * 4 * (ns - 1) + 3 * ns - 5;
		lc1[1] = ecount * 4 * (ns - 1) + 4 * ns - 5;
		led.push_back(lc1);
		ecount++;*/
	}

	string fname = fn + "_disp.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i<spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 9 * sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "8 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3]
				<< " " << sele[i][4] << " " << sele[i][5] << " " << sele[i][6] << " " << sele[i][7] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "12\n";
		}
		//fout<<"\nPOINT_DATA "<<sdisp.size()<<"\nVECTORS disp float\n";
		//for(i=0;i<sdisp.size();i++)
		//{
		//	fout << sdisp[i][0] << " " << sdisp[i][1] << " " << sdisp[i][2] << "\n";
		//}
		//fout << "\nPOINT_DATA " << sse.size() << "\nVECTORS strain float\n";
		//for (i = 0; i<sse.size(); i++)
		//{
		//	fout << sse[i][0] << " " << sse[i][1] << " " << sse[i][2] << "\n";
		//}
		//fout<<"POINT_DATA "<<sss.size()<<"\nVECTORS stress float\n";
		//for(i=0;i<sss.size();i++)
		//{
		//	fout<<sss[i][0]<<" "<<sss[i][1]<<" "<<sss[i][2]<<"\n";
		//}
		//fout<<"POINT_DATA "<<sdisp_err.size()<<"\nVECTORS disp_err float\n";
		//for(i=0;i<sdisp_err.size();i++)
		//{
		//	fout<<sdisp_err[i][0]<<" "<<sdisp_err[i][1]<<" 0\n";
		//}

		//fout<<"POINT_DATA "<<sdisp.size()<<"\nSCALARS u float 1\nLOOKUP_TABLE default\n";
		//for(uint i=0;i<sdisp.size();i++)
		//{
		//	fout<<sdisp[i]<<"\n";
		//}

		//fout<<"\nCELL_DATA "<<errL2.size()<<"\nSCALARS Error float 1\nLOOKUP_TABLE default\n";
		//for(i=0;i<errL2.size();i++)
		//{
		//	fout<<errL2[i]<<"\n";
		//	//fout<<eles[i].type<<"\n";
		//}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	string fname1(fn+"-lines.vtk");
	ofstream fout1;
	fout1.open(fname1.c_str());
	if(fout1.is_open())
	{
		fout1<<"# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout1<<"POINTS "<<lpt.size()<<" float\n";
		for(i=0;i<lpt.size();i++)
		{
			fout1<<lpt[i][0]<<" "<<lpt[i][1]<<" "<<lpt[i][2]<<"\n";
		}
		fout1<<"\nCELLS "<<led.size()<<" "<<3*led.size()<<'\n';
		for(i=0;i<led.size();i++)
		{
			fout1<<"2 "<<led[i][0]<<" "<<led[i][1]<<'\n';
		}
		fout1<<"\nCELL_TYPES "<<led.size()<<'\n';
		for(i=0;i<led.size();i++)
		{
			fout1<<"3\n";
		}
		fout1.close();
	}
	else
	{
		cout<<"Cannot open "<<fname1<<"!\n";
	}
}

void Laplace::VisualizeVTK_1(const vector<BezierElement3D>& bzmesh, string fn)
{
	vector<array<double, 3>> spt;//sample points
	vector<double> sdisp;
	//vector<array<double, 3>> sdisp_err;
	//vector<array<double, 3>> sse;
	//vector<array<double, 3>> sss;
	vector<array<int, 8>> sele;
	vector<double> errL2;
	vector<array<double, 3>> lpt;//visulize parameter lines
	vector<array<int, 2>> led;//line connectivity
	//int ecount(0);
	//vector<double> su;
	//for (int i = 0; i<ns; i++)
	//{
	//	su[i] = double(i) / (double(ns) - 1.);
	//}

	for (unsigned int e = 0; e<bzmesh.size(); e++)
	{
		//double errtmp;
		//ElementError(bzmesh[e],errtmp);
		//errL2.push_back(sqrt(errtmp));

		//if (bzmesh[e].type == 1)
		{
			int ns(2);
			if (bzmesh[e].type == 1) ns = 5;
			vector<double> su(ns);
			for (int i = 0; i<ns; i++)
			{
				su[i] = double(i) / (double(ns) - 1.);
			}

			int loc(0);
			int pstart = spt.size();
			for (int a = 0; a<ns; a++)
			{
				for (int b = 0; b<ns; b++)
				{
					for (int c = 0; c < ns; c++)
					{
						double pt1[3];
						double disp, detmp;
						bzmesh[e].Para2Phys(su[c], su[b], su[a], pt1);
						DispCal(su[c], su[b], su[a], bzmesh[e], disp, detmp);
						array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
						spt.push_back(pt);
						sdisp.push_back(disp);
					}
				}
			}
			int nns[2] = { ns*ns*ns, ns*ns };
			for (int a = 0; a<ns - 1; a++)
			{
				for (int b = 0; b<ns - 1; b++)
				{
					for (int c = 0; c < ns - 1; c++)
					{
						array<int, 8> el;
						el[0] = pstart + a*nns[1] + b*ns + c;
						el[1] = pstart + a*nns[1] + b*ns + c + 1;
						el[2] = pstart + a*nns[1] + (b + 1)*ns + c + 1;
						el[3] = pstart + a*nns[1] + (b + 1)*ns + c;
						el[4] = pstart + (a + 1)*nns[1] + b*ns + c;
						el[5] = pstart + (a + 1)*nns[1] + b*ns + c + 1;
						el[6] = pstart + (a + 1)*nns[1] + (b + 1)*ns + c + 1;
						el[7] = pstart + (a + 1)*nns[1] + (b + 1)*ns + c;
						sele.push_back(el);
					}
				}
			}
			//edges
			int lstart = lpt.size();
			for (int a = 0; a < ns; a++)
			{
				double pt1[3];
				bzmesh[e].Para2Phys(su[a], 0., 0., pt1);
				array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
				lpt.push_back(pt);
			}
			for (int a = 0; a < ns - 1; a++)
			{
				array<int, 2> ed = { lstart + a, lstart + a + 1 };
				led.push_back(ed);
			}
			lstart = lpt.size();
			for (int a = 0; a < ns; a++)
			{
				double pt1[3];
				bzmesh[e].Para2Phys(su[a], 1., 0., pt1);
				array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
				lpt.push_back(pt);
			}
			for (int a = 0; a < ns - 1; a++)
			{
				array<int, 2> ed = { lstart + a, lstart + a + 1 };
				led.push_back(ed);
			}
			lstart = lpt.size();
			for (int a = 0; a < ns; a++)
			{
				double pt1[3];
				bzmesh[e].Para2Phys(su[a], 0., 1., pt1);
				array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
				lpt.push_back(pt);
			}
			for (int a = 0; a < ns - 1; a++)
			{
				array<int, 2> ed = { lstart + a, lstart + a + 1 };
				led.push_back(ed);
			}
			lstart = lpt.size();
			for (int a = 0; a < ns; a++)
			{
				double pt1[3];
				bzmesh[e].Para2Phys(su[a], 1., 1., pt1);
				array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
				lpt.push_back(pt);
			}
			for (int a = 0; a < ns - 1; a++)
			{
				array<int, 2> ed = { lstart + a, lstart + a + 1 };
				led.push_back(ed);
			}
			lstart = lpt.size();
			for (int a = 0; a < ns; a++)
			{
				double pt1[3];
				bzmesh[e].Para2Phys(0., su[a], 0., pt1);
				array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
				lpt.push_back(pt);
			}
			for (int a = 0; a < ns - 1; a++)
			{
				array<int, 2> ed = { lstart + a, lstart + a + 1 };
				led.push_back(ed);
			}
			lstart = lpt.size();
			for (int a = 0; a < ns; a++)
			{
				double pt1[3];
				bzmesh[e].Para2Phys(1., su[a], 0., pt1);
				array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
				lpt.push_back(pt);
			}
			for (int a = 0; a < ns - 1; a++)
			{
				array<int, 2> ed = { lstart + a, lstart + a + 1 };
				led.push_back(ed);
			}
			lstart = lpt.size();
			for (int a = 0; a < ns; a++)
			{
				double pt1[3];
				bzmesh[e].Para2Phys(0., su[a], 1., pt1);
				array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
				lpt.push_back(pt);
			}
			for (int a = 0; a < ns - 1; a++)
			{
				array<int, 2> ed = { lstart + a, lstart + a + 1 };
				led.push_back(ed);
			}
			lstart = lpt.size();
			for (int a = 0; a < ns; a++)
			{
				double pt1[3];
				bzmesh[e].Para2Phys(1., su[a], 1., pt1);
				array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
				lpt.push_back(pt);
			}
			for (int a = 0; a < ns - 1; a++)
			{
				array<int, 2> ed = { lstart + a, lstart + a + 1 };
				led.push_back(ed);
			}
			lstart = lpt.size();
			for (int a = 0; a < ns; a++)
			{
				double pt1[3];
				bzmesh[e].Para2Phys(0., 0., su[a], pt1);
				array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
				lpt.push_back(pt);
			}
			for (int a = 0; a < ns - 1; a++)
			{
				array<int, 2> ed = { lstart + a, lstart + a + 1 };
				led.push_back(ed);
			}
			lstart = lpt.size();
			for (int a = 0; a < ns; a++)
			{
				double pt1[3];
				bzmesh[e].Para2Phys(1., 0., su[a], pt1);
				array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
				lpt.push_back(pt);
			}
			for (int a = 0; a < ns - 1; a++)
			{
				array<int, 2> ed = { lstart + a, lstart + a + 1 };
				led.push_back(ed);
			}
			lstart = lpt.size();
			for (int a = 0; a < ns; a++)
			{
				double pt1[3];
				bzmesh[e].Para2Phys(0., 1., su[a], pt1);
				array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
				lpt.push_back(pt);
			}
			for (int a = 0; a < ns - 1; a++)
			{
				array<int, 2> ed = { lstart + a, lstart + a + 1 };
				led.push_back(ed);
			}
			lstart = lpt.size();
			for (int a = 0; a < ns; a++)
			{
				double pt1[3];
				bzmesh[e].Para2Phys(1., 1., su[a], pt1);
				array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
				lpt.push_back(pt);
			}
			for (int a = 0; a < ns - 1; a++)
			{
				array<int, 2> ed = { lstart + a, lstart + a + 1 };
				led.push_back(ed);
			}
		}
	}

	string fname = fn + "_disp.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i<spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 9 * sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "8 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3]
				<< " " << sele[i][4] << " " << sele[i][5] << " " << sele[i][6] << " " << sele[i][7] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "12\n";
		}
		//fout<<"\nPOINT_DATA "<<sdisp.size()<<"\nVECTORS disp float\n";
		//for(i=0;i<sdisp.size();i++)
		//{
		//	fout << sdisp[i][0] << " " << sdisp[i][1] << " " << sdisp[i][2] << "\n";
		//}
		//fout << "\nPOINT_DATA " << sse.size() << "\nVECTORS strain float\n";
		//for (i = 0; i<sse.size(); i++)
		//{
		//	fout << sse[i][0] << " " << sse[i][1] << " " << sse[i][2] << "\n";
		//}
		//fout<<"POINT_DATA "<<sss.size()<<"\nVECTORS stress float\n";
		//for(i=0;i<sss.size();i++)
		//{
		//	fout<<sss[i][0]<<" "<<sss[i][1]<<" "<<sss[i][2]<<"\n";
		//}
		//fout<<"POINT_DATA "<<sdisp_err.size()<<"\nVECTORS disp_err float\n";
		//for(i=0;i<sdisp_err.size();i++)
		//{
		//	fout<<sdisp_err[i][0]<<" "<<sdisp_err[i][1]<<" 0\n";
		//}

		fout << "POINT_DATA " << sdisp.size() << "\nSCALARS err float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<sdisp.size(); i++)
		{
			fout << sdisp[i] << "\n";
		}

		//fout<<"\nCELL_DATA "<<errL2.size()<<"\nSCALARS Error float 1\nLOOKUP_TABLE default\n";
		//for(i=0;i<errL2.size();i++)
		//{
		//	fout<<errL2[i]<<"\n";
		//	//fout<<eles[i].type<<"\n";
		//}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	string fname1(fn+"-lines.vtk");
	ofstream fout1;
	fout1.open(fname1.c_str());
	if(fout1.is_open())
	{
		fout1<<"# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout1<<"POINTS "<<lpt.size()<<" float\n";
		for(i=0;i<lpt.size();i++)
		{
			fout1<<lpt[i][0]<<" "<<lpt[i][1]<<" "<<lpt[i][2]<<"\n";
		}
		fout1<<"\nCELLS "<<led.size()<<" "<<3*led.size()<<'\n';
		for(i=0;i<led.size();i++)
		{
			fout1<<"2 "<<led[i][0]<<" "<<led[i][1]<<'\n';
		}
		fout1<<"\nCELL_TYPES "<<led.size()<<'\n';
		for(i=0;i<led.size();i++)
		{
			fout1<<"3\n";
		}
		fout1.close();
	}
	else
	{
		cout<<"Cannot open "<<fname1<<"!\n";
	}
}

void Laplace::VisualizeVTK_2(const vector<BezierElement3D>& bzmesh, string fn)
{
	vector<array<double, 3>> spt;//sample points
	vector<double> sdisp;
	//vector<array<double, 3>> sdisp_err;
	//vector<array<double, 3>> sse;
	//vector<array<double, 3>> sss;
	vector<array<int, 8>> sele;
	vector<double> errL2;
	vector<array<double, 3>> lpt;//visulize parameter lines
	vector<array<int, 2>> led;//line connectivity
	//int ecount(0);
	//vector<double> su;
	//for (int i = 0; i<ns; i++)
	//{
	//	su[i] = double(i) / (double(ns) - 1.);
	//}

	int corner[8] = {0,3,15,12,48,51,63,60};

	for (unsigned int e = 0; e<bzmesh.size(); e++)
	{
		//double errtmp;
		//ElementError(bzmesh[e],errtmp);
		//errL2.push_back(sqrt(errtmp));

		double ctp[3] = { 0., 0., 0. };
		for (int j = 0; j < 8; j++)
		{
			ctp[0] += bzmesh[e].pts[corner[j]][0];
			ctp[1] += bzmesh[e].pts[corner[j]][1];
			ctp[2] += bzmesh[e].pts[corner[j]][2];
		}
		ctp[0] /= 8.; ctp[1] /= 8.; ctp[2] /= 8.;
		if (ctp[0] >= rmv[0][0] && ctp[0] <= rmv[0][1] && ctp[1] >= rmv[1][0] && ctp[1] <= rmv[1][1] && ctp[2] >= rmv[2][0] && ctp[2] <= rmv[2][1])
		{
			continue;
		}
		else
		{
			int ns(2);
			if (bzmesh[e].type == 1) ns = 5;
			vector<double> su(ns);
			for (int i = 0; i<ns; i++)
			{
				su[i] = double(i) / (double(ns) - 1.);
			}

			int loc(0);
			int pstart = spt.size();
			for (int a = 0; a<ns; a++)
			{
				for (int b = 0; b<ns; b++)
				{
					for (int c = 0; c < ns; c++)
					{
						double pt1[3];
						double disp, detmp;
						bzmesh[e].Para2Phys(su[c], su[b], su[a], pt1);
						DispCal(su[c], su[b], su[a], bzmesh[e], disp, detmp);
						array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
						//disp = exact_sol(pt);//tmp
						spt.push_back(pt);
						sdisp.push_back(disp);
					}
				}
			}
			int nns[2] = { ns*ns*ns, ns*ns };
			for (int a = 0; a<ns - 1; a++)
			{
				for (int b = 0; b<ns - 1; b++)
				{
					for (int c = 0; c < ns - 1; c++)
					{
						array<int, 8> el;
						el[0] = pstart + a*nns[1] + b*ns + c;
						el[1] = pstart + a*nns[1] + b*ns + c + 1;
						el[2] = pstart + a*nns[1] + (b + 1)*ns + c + 1;
						el[3] = pstart + a*nns[1] + (b + 1)*ns + c;
						el[4] = pstart + (a + 1)*nns[1] + b*ns + c;
						el[5] = pstart + (a + 1)*nns[1] + b*ns + c + 1;
						el[6] = pstart + (a + 1)*nns[1] + (b + 1)*ns + c + 1;
						el[7] = pstart + (a + 1)*nns[1] + (b + 1)*ns + c;
						sele.push_back(el);
					}
				}
			}
			//edges
			int lstart = lpt.size();
			for (int a = 0; a < ns; a++)
			{
				double pt1[3];
				bzmesh[e].Para2Phys(su[a], 0., 0., pt1);
				array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
				lpt.push_back(pt);
			}
			for (int a = 0; a < ns - 1; a++)
			{
				array<int, 2> ed = { lstart + a, lstart + a + 1 };
				led.push_back(ed);
			}
			lstart = lpt.size();
			for (int a = 0; a < ns; a++)
			{
				double pt1[3];
				bzmesh[e].Para2Phys(su[a], 1., 0., pt1);
				array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
				lpt.push_back(pt);
			}
			for (int a = 0; a < ns - 1; a++)
			{
				array<int, 2> ed = { lstart + a, lstart + a + 1 };
				led.push_back(ed);
			}
			lstart = lpt.size();
			for (int a = 0; a < ns; a++)
			{
				double pt1[3];
				bzmesh[e].Para2Phys(su[a], 0., 1., pt1);
				array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
				lpt.push_back(pt);
			}
			for (int a = 0; a < ns - 1; a++)
			{
				array<int, 2> ed = { lstart + a, lstart + a + 1 };
				led.push_back(ed);
			}
			lstart = lpt.size();
			for (int a = 0; a < ns; a++)
			{
				double pt1[3];
				bzmesh[e].Para2Phys(su[a], 1., 1., pt1);
				array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
				lpt.push_back(pt);
			}
			for (int a = 0; a < ns - 1; a++)
			{
				array<int, 2> ed = { lstart + a, lstart + a + 1 };
				led.push_back(ed);
			}
			lstart = lpt.size();
			for (int a = 0; a < ns; a++)
			{
				double pt1[3];
				bzmesh[e].Para2Phys(0., su[a], 0., pt1);
				array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
				lpt.push_back(pt);
			}
			for (int a = 0; a < ns - 1; a++)
			{
				array<int, 2> ed = { lstart + a, lstart + a + 1 };
				led.push_back(ed);
			}
			lstart = lpt.size();
			for (int a = 0; a < ns; a++)
			{
				double pt1[3];
				bzmesh[e].Para2Phys(1., su[a], 0., pt1);
				array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
				lpt.push_back(pt);
			}
			for (int a = 0; a < ns - 1; a++)
			{
				array<int, 2> ed = { lstart + a, lstart + a + 1 };
				led.push_back(ed);
			}
			lstart = lpt.size();
			for (int a = 0; a < ns; a++)
			{
				double pt1[3];
				bzmesh[e].Para2Phys(0., su[a], 1., pt1);
				array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
				lpt.push_back(pt);
			}
			for (int a = 0; a < ns - 1; a++)
			{
				array<int, 2> ed = { lstart + a, lstart + a + 1 };
				led.push_back(ed);
			}
			lstart = lpt.size();
			for (int a = 0; a < ns; a++)
			{
				double pt1[3];
				bzmesh[e].Para2Phys(1., su[a], 1., pt1);
				array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
				lpt.push_back(pt);
			}
			for (int a = 0; a < ns - 1; a++)
			{
				array<int, 2> ed = { lstart + a, lstart + a + 1 };
				led.push_back(ed);
			}
			lstart = lpt.size();
			for (int a = 0; a < ns; a++)
			{
				double pt1[3];
				bzmesh[e].Para2Phys(0., 0., su[a], pt1);
				array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
				lpt.push_back(pt);
			}
			for (int a = 0; a < ns - 1; a++)
			{
				array<int, 2> ed = { lstart + a, lstart + a + 1 };
				led.push_back(ed);
			}
			lstart = lpt.size();
			for (int a = 0; a < ns; a++)
			{
				double pt1[3];
				bzmesh[e].Para2Phys(1., 0., su[a], pt1);
				array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
				lpt.push_back(pt);
			}
			for (int a = 0; a < ns - 1; a++)
			{
				array<int, 2> ed = { lstart + a, lstart + a + 1 };
				led.push_back(ed);
			}
			lstart = lpt.size();
			for (int a = 0; a < ns; a++)
			{
				double pt1[3];
				bzmesh[e].Para2Phys(0., 1., su[a], pt1);
				array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
				lpt.push_back(pt);
			}
			for (int a = 0; a < ns - 1; a++)
			{
				array<int, 2> ed = { lstart + a, lstart + a + 1 };
				led.push_back(ed);
			}
			lstart = lpt.size();
			for (int a = 0; a < ns; a++)
			{
				double pt1[3];
				bzmesh[e].Para2Phys(1., 1., su[a], pt1);
				array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
				lpt.push_back(pt);
			}
			for (int a = 0; a < ns - 1; a++)
			{
				array<int, 2> ed = { lstart + a, lstart + a + 1 };
				led.push_back(ed);
			}
		}
	}

	string fname = fn + "_disp.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i<spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 9 * sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "8 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3]
				<< " " << sele[i][4] << " " << sele[i][5] << " " << sele[i][6] << " " << sele[i][7] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "12\n";
		}
		//fout<<"\nPOINT_DATA "<<sdisp.size()<<"\nVECTORS disp float\n";
		//for(i=0;i<sdisp.size();i++)
		//{
		//	fout << sdisp[i][0] << " " << sdisp[i][1] << " " << sdisp[i][2] << "\n";
		//}
		//fout << "\nPOINT_DATA " << sse.size() << "\nVECTORS strain float\n";
		//for (i = 0; i<sse.size(); i++)
		//{
		//	fout << sse[i][0] << " " << sse[i][1] << " " << sse[i][2] << "\n";
		//}
		//fout<<"POINT_DATA "<<sss.size()<<"\nVECTORS stress float\n";
		//for(i=0;i<sss.size();i++)
		//{
		//	fout<<sss[i][0]<<" "<<sss[i][1]<<" "<<sss[i][2]<<"\n";
		//}
		//fout<<"POINT_DATA "<<sdisp_err.size()<<"\nVECTORS disp_err float\n";
		//for(i=0;i<sdisp_err.size();i++)
		//{
		//	fout<<sdisp_err[i][0]<<" "<<sdisp_err[i][1]<<" 0\n";
		//}

		fout << "POINT_DATA " << sdisp.size() << "\nSCALARS err float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<sdisp.size(); i++)
		{
			fout << sdisp[i] << "\n";
		}

		//fout<<"\nCELL_DATA "<<errL2.size()<<"\nSCALARS Error float 1\nLOOKUP_TABLE default\n";
		//for(i=0;i<errL2.size();i++)
		//{
		//	fout<<errL2[i]<<"\n";
		//	//fout<<eles[i].type<<"\n";
		//}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	string fname1(fn + "-lines.vtk");
	ofstream fout1;
	fout1.open(fname1.c_str());
	if (fout1.is_open())
	{
		fout1 << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout1 << "POINTS " << lpt.size() << " float\n";
		for (i = 0; i<lpt.size(); i++)
		{
			fout1 << lpt[i][0] << " " << lpt[i][1] << " " << lpt[i][2] << "\n";
		}
		fout1 << "\nCELLS " << led.size() << " " << 3 * led.size() << '\n';
		for (i = 0; i<led.size(); i++)
		{
			fout1 << "2 " << led[i][0] << " " << led[i][1] << '\n';
		}
		fout1 << "\nCELL_TYPES " << led.size() << '\n';
		for (i = 0; i<led.size(); i++)
		{
			fout1 << "3\n";
		}
		fout1.close();
	}
	else
	{
		cout << "Cannot open " << fname1 << "!\n";
	}
}

void Laplace::VisualizeError(const vector<BezierElement3D>& bzmesh, const vector<double>& err, string fn)
{
	int cn[8] = {0,3,15,12,48,51,63,60};
	string fname = fn + "_err.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i,j;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << 8*bzmesh.size() << " float\n";
		for (i = 0; i<bzmesh.size(); i++)
		{
			for (j = 0; j < 8; j++)
			{
				fout << bzmesh[i].pts[cn[j]][0] << " " << bzmesh[i].pts[cn[j]][1] << " " << bzmesh[i].pts[cn[j]][2] << "\n";
			}
		}
		fout << "\nCELLS " << bzmesh.size() << " " << 9 * bzmesh.size() << '\n';
		for (i = 0; i<bzmesh.size(); i++)
		{
			fout << "8 " << 8 * i << " " << 8 * i + 1 << " " << 8 * i + 2 << " " << 8 * i + 3
				<< " " << 8 * i + 4 << " " << 8 * i + 5 << " " << 8 * i + 6 << " " << 8 * i + 7 << '\n';
		}
		fout << "\nCELL_TYPES " << bzmesh.size() << '\n';
		for (i = 0; i<bzmesh.size(); i++)
		{
			fout << "12\n";
		}
		//fout << "POINT_DATA " << sdisp.size() << "\nSCALARS err float 1\nLOOKUP_TABLE default\n";
		//for (uint i = 0; i<sdisp.size(); i++)
		//{
		//	fout << sdisp[i] << "\n";
		//}
		fout<<"\nCELL_DATA "<<err.size()<<"\nSCALARS Error float 1\nLOOKUP_TABLE default\n";
		for(i=0;i<err.size();i++)
		{
			fout<<sqrt(err[i])<<"\n";
			//fout<<eles[i].type<<"\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	//string fname1(fn+"-lines.vtk");
	//ofstream fout1;
	//fout1.open(fname1.c_str());
	//if(fout1.is_open())
	//{
	//	fout1<<"# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
	//	fout1<<"POINTS "<<lpt.size()<<" float\n";
	//	for(i=0;i<lpt.size();i++)
	//	{
	//		fout1<<lpt[i][0]<<" "<<lpt[i][1]<<" "<<lpt[i][2]<<"\n";
	//	}
	//	fout1<<"\nCELLS "<<led.size()<<" "<<3*led.size()<<'\n';
	//	for(i=0;i<led.size();i++)
	//	{
	//		fout1<<"2 "<<led[i][0]<<" "<<led[i][1]<<'\n';
	//	}
	//	fout1<<"\nCELL_TYPES "<<led.size()<<'\n';
	//	for(i=0;i<led.size();i++)
	//	{
	//		fout1<<"3\n";
	//	}
	//	fout1.close();
	//}
	//else
	//{
	//	cout<<"Cannot open "<<fname1<<"!\n";
	//}
}

void Laplace::VisualizeError_BezierCouple(const vector<BezierElement3D>& bzmesh, const vector<double>& err, string fn)
{
	int cn[8] = { 0, 3, 15, 12, 48, 51, 63, 60 };
	string fname1 = fn + "_err_reg.vtk";
	string fname2 = fn + "_err_irr.vtk";
	ofstream fout;
	fout.open(fname1.c_str());
	unsigned int i, j;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		uint nel(0);
		for (i = 0; i < bzmesh.size(); i++)
		{
			//if (bzmesh[i].bzflag == 0) nel++;
			if (bzmesh[i].type != 1) nel++;
		}
		fout << "POINTS " << 8 * nel << " float\n";
		for (i = 0; i<bzmesh.size(); i++)
		{
			//if (bzmesh[i].bzflag == 0)
			if (bzmesh[i].type != 1)
			{
				for (j = 0; j < 8; j++)
				{
					fout << bzmesh[i].pts[cn[j]][0] << " " << bzmesh[i].pts[cn[j]][1] << " " << bzmesh[i].pts[cn[j]][2] << "\n";
				}
			}
		}
		fout << "\nCELLS " << nel << " " << 9 * nel << '\n';
		for (i = 0; i<nel; i++)
		{
			fout << "8 " << 8 * i << " " << 8 * i + 1 << " " << 8 * i + 2 << " " << 8 * i + 3
				<< " " << 8 * i + 4 << " " << 8 * i + 5 << " " << 8 * i + 6 << " " << 8 * i + 7 << '\n';
		}
		fout << "\nCELL_TYPES " << nel << '\n';
		for (i = 0; i<nel; i++)
		{
			fout << "12\n";
		}
		//fout << "POINT_DATA " << sdisp.size() << "\nSCALARS err float 1\nLOOKUP_TABLE default\n";
		//for (uint i = 0; i<sdisp.size(); i++)
		//{
		//	fout << sdisp[i] << "\n";
		//}
		fout << "\nCELL_DATA " << nel << "\nSCALARS Error float 1\nLOOKUP_TABLE default\n";
		for (i = 0; i<err.size(); i++)
		{
			//if (bzmesh[i].bzflag == 0)
			if (bzmesh[i].type != 1)
			{
				fout << sqrt(err[i]) << "\n";
				//fout<<eles[i].type<<"\n";
			}
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname1 << "!\n";
	}

	fout.open(fname2.c_str());
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		uint nel(0);
		for (i = 0; i < bzmesh.size(); i++)
		{
			//if (bzmesh[i].bzflag != 0) nel++;
			if (bzmesh[i].type == 1) nel++;
		}
		fout << "POINTS " << 8 * nel << " float\n";
		for (i = 0; i<bzmesh.size(); i++)
		{
			//if (bzmesh[i].bzflag != 0)
			if (bzmesh[i].type == 1)
			{
				for (j = 0; j < 8; j++)
				{
					fout << bzmesh[i].pts[cn[j]][0] << " " << bzmesh[i].pts[cn[j]][1] << " " << bzmesh[i].pts[cn[j]][2] << "\n";
				}
			}
		}
		fout << "\nCELLS " << nel << " " << 9 * nel << '\n';
		for (i = 0; i<nel; i++)
		{
			fout << "8 " << 8 * i << " " << 8 * i + 1 << " " << 8 * i + 2 << " " << 8 * i + 3
				<< " " << 8 * i + 4 << " " << 8 * i + 5 << " " << 8 * i + 6 << " " << 8 * i + 7 << '\n';
		}
		fout << "\nCELL_TYPES " << nel << '\n';
		for (i = 0; i<nel; i++)
		{
			fout << "12\n";
		}
		//fout << "POINT_DATA " << sdisp.size() << "\nSCALARS err float 1\nLOOKUP_TABLE default\n";
		//for (uint i = 0; i<sdisp.size(); i++)
		//{
		//	fout << sdisp[i] << "\n";
		//}
		fout << "\nCELL_DATA " << nel << "\nSCALARS Error float 1\nLOOKUP_TABLE default\n";
		for (i = 0; i<err.size(); i++)
		{
			//if (bzmesh[i].bzflag != 0)
			if (bzmesh[i].type == 1)
			{
				fout << sqrt(err[i]) << "\n";
				//fout<<eles[i].type<<"\n";
			}
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname2 << "!\n";
	}
}

void Laplace::Run(const vector<BezierElement3D>& bzmesh, string fn, vector<double>& err)
{
	GaussInfo(4);
	//FindRange(bzmesh);
	//VisualizeVTK(bzmesh, fn);
	SparseMatrix<double> GK(neq,neq);
	//MatrixXd GK(neq, neq);
	VectorXd GF(neq);
	//GK.setZero();
	GF.setZero();
	cout<<"Building linear system...\n";
	//BuildLinearSystem(bzmesh,GK,GF);
	BuildLinearSystem_1(bzmesh, GK, GF);
	//BuildLinearSystem_Dense(bzmesh, GK, GF);
	//OutputMatrix(GK,fn);
	Solver(GK,GF);
	//Solver_Dense(GK, GF);
	//VisualizeVTK(bzmesh,fn);
	VisualizeVTK_1(bzmesh, fn);//with smoother boundary
	//VisualizeVTK_2(bzmesh, fn);//smooth boundary, remove elements
	//ErrorEstimate(bzmesh,err);//with exact solution
	//TotalErrorEstimate(bzmesh, err);//without exact solution
	//VisualizeError(bzmesh, err, fn);
	cout << "\nDone run Laplace...\n";
}

void Laplace::BasisFunction_1(double u, double v, double w, const BezierElement3D& bzel, vector<double>& Nx, vector<array<double, 3>>& dNdx, double& detJ)
{
	double Nx0[64];
	double dNdx0[64][3];
	BasisFunction(u, v, w, bzel.pts, Nx0, dNdx0, detJ);
//#pragma omp parallel for
	for (int i = 0; i < bzel.IEN.size(); i++)
	{
		Nx[i] = 0.;
		dNdx[i][0] = 0.; dNdx[i][1] = 0.; dNdx[i][2] = 0.;
		for (int j = 0; j < 64; j++)
		{
			Nx[i] += bzel.cmat[i][j] * Nx0[j];
			dNdx[i][0] += bzel.cmat[i][j] * dNdx0[j][0];
			dNdx[i][1] += bzel.cmat[i][j] * dNdx0[j][1];
			dNdx[i][2] += bzel.cmat[i][j] * dNdx0[j][2];
		}
	}
}

void Laplace::ElementMatrix_1(const vector<array<double, 3>>& dNdx, double detJ, vector<vector<double>>& EK)
{
	int i, j;
//#pragma omp parallel for
	for (i = 0; i < dNdx.size(); i++)
	{
		for (j = 0; j < dNdx.size(); j++)
		{
			EK[i][j] += (dNdx[i][0] * dNdx[j][0] + dNdx[i][1] * dNdx[j][1] + dNdx[i][2] * dNdx[j][2])*detJ;
		}
	}
}

void Laplace::ElementForce_1(const vector<double>& Nx, double Jmod, double Fb, vector<double>& EF)
{
	for (uint i = 0; i < Nx.size(); i++)
	{
		EF[i] += Nx[i] * Jmod * Fb;
	}
}

void Laplace::Assembly_1(const vector<vector<double>>& EK, const vector<double>& EF, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF)
{
	//unsigned int i, j, A, B;
	int i, j, A, B;
	//cout << "before global\n";
	for (i = 0; i<IEN.size(); i++)
	{
		for (j = 0; j<IEN.size(); j++)
		{
			A = IEN[i]; B = IEN[j];
			//if (A >= IDBC.size() || B >= IDBC.size())
			//{
			//	cout << "out of range!\n";
			//	cout << "A: " << A << "/" << IDBC.size() << "\n";
			//	cout << "B: " << B << "/" << IDBC.size() << "\n";
			//	getchar();
			//}
			if (IDBC[A] != -1 && IDBC[B] != -1)
			{
				GK.coeffRef(IDBC[A], IDBC[B]) += EK[i][j];
			}
			else if (IDBC[A] != -1 && IDBC[B] == -1)
			{
				GF(IDBC[A]) -= EK[i][j] * gh[B];
			}
		}
	}
	for (i = 0; i < IEN.size(); i++)
	{
		A = IEN[i];
		if (IDBC[A] != -1)
		{
			GF(IDBC[A]) += EF[i];
		}
	}
	//cout << "after global\n";
}

void Laplace::Assembly_Coupling(const vector<vector<double>>& EK, const vector<double>& EF, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF)
{
	int i, j, A, B;
	for (i = 0; i<IEN.size(); i++)
	{
		for (j = 0; j<IEN.size(); j++)
		{
			A = IEN[i]; B = IEN[j];
			if (A != -1 && B != -1)
			{
				if (IDBC[A] != -1 && IDBC[B] != -1)
				{
					GK.coeffRef(IDBC[A], IDBC[B]) += EK[i][j];
				}
			}
			//else if (IDBC[A] != -1 && IDBC[B] == -1)
			//{
			//	GF(IDBC[A]) -= EK[i][j] * gh[B];
			//}
		}
	}
	for (i = 0; i < IEN.size(); i++)
	{
		A = IEN[i];
		if (A != -1 && IDBC[A] != -1)
		{
			GF(IDBC[A]) += EF[i];
		}
	}
}

void Laplace::Assembly_Dense(const vector<vector<double>>& EK, const vector<double>& EF, const vector<int>& IEN, MatrixXd& GK, VectorXd& GF)
{
	//unsigned int i, j, A, B;
	int i, j, A, B;
	//cout << "before global\n";
	for (i = 0; i<IEN.size(); i++)
	{
		for (j = 0; j<IEN.size(); j++)
		{
			A = IEN[i]; B = IEN[j];
			//if (A >= IDBC.size() || B >= IDBC.size())
			//{
			//	cout << "out of range!\n";
			//	cout << "A: " << A << "/" << IDBC.size() << "\n";
			//	cout << "B: " << B << "/" << IDBC.size() << "\n";
			//	getchar();
			//}
			if (IDBC[A] != -1 && IDBC[B] != -1)
			{
				GK(IDBC[A], IDBC[B]) += EK[i][j];
			}
			else if (IDBC[A] != -1 && IDBC[B] == -1)
			{
				GF(IDBC[A]) -= EK[i][j] * gh[B];
			}
		}
	}
	for (i = 0; i < IEN.size(); i++)
	{
		A = IEN[i];
		if (IDBC[A] != -1)
		{
			GF(IDBC[A]) += EF[i];
		}
	}
	//cout << "after global\n";
}

void Laplace::BuildLinearSystem_1(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF)
{
	cout << "Initialize sparse matrix...\n";
	cout << "# of eles: " << bzmesh.size() << "\n";
	int nresv(6 * 64);
	if (nresv < neq)
	{
		GK.reserve(VectorXi::Constant(neq, nresv));
	}
	else
	{
		GK.reserve(VectorXi::Constant(neq, neq));
	}
	for (int e = 0; e < bzmesh.size(); e++)
	{
		//if (e != 0 && e % 500 == 0)
		//{
		//	cout << e << " ";
		//}
		for (int i = 0; i < bzmesh[e].IEN.size(); i++)
		{
			for (int j = 0; j < bzmesh[e].IEN.size(); j++)
			{
				int A(bzmesh[e].IEN[i]), B(bzmesh[e].IEN[j]);
				if (IDBC[A] != -1 && IDBC[B] != -1)
				{
					GK.coeffRef(IDBC[A], IDBC[B]) = 0.;
				}
			}
		}
	}
	GK.makeCompressed();
	cout << "done initializing\n";
	cout << "npt: " << npt << "\n";
	cout << "neq: " << neq << "\n";
	cout << "# of eles: " << bzmesh.size() << "\n";
	//clock_t t0 = clock();
#pragma omp parallel for
	for (int e = 0; e<bzmesh.size(); e++)
	{
		if (e != 0 && e % 500 == 0)
		{
			cout << e << " ";
		}
		vector<vector<double>> EK(bzmesh[e].IEN.size(), vector<double>(bzmesh[e].IEN.size(), 0.));
		vector<double> EF(bzmesh[e].IEN.size(), 0.);
		vector<double> Nx(bzmesh[e].IEN.size());
		vector<array<double, 3>> dNdx(bzmesh[e].IEN.size());
		double detJ, Fb;
		array<double, 3> pt;
		for (int i = 0; i<Gpt.size(); i++)
		{
			for (int j = 0; j<Gpt.size(); j++)
			{
				for (int k = 0; k < Gpt.size(); k++)
				{
					BasisFunction_1(Gpt[i], Gpt[j], Gpt[k], bzmesh[e], Nx, dNdx, detJ);
					detJ = wght[i] * wght[j] * wght[k] * detJ;
					ElementMatrix_1(dNdx, detJ, EK);

					bzmesh[e].Para2Phys(Gpt[i], Gpt[j], Gpt[k], pt.data());
					//Fb = f_source(pt);
					Fb = f_source_lap(pt);
					ElementForce_1(Nx, detJ, Fb, EF);
				}
			}
		}
		Assembly_1(EK, EF, bzmesh[e].IEN, GK, GF);
	}
	cout << "Done assembly!\n";
	//clock_t t1 = clock();
	//double elaps = double(t1 - t0) / CLOCKS_PER_SEC;
	//cout << "\ntime of basis: " << elaps << "\n";
}

void Laplace::BuildLinearSystem_Dense(const vector<BezierElement3D>& bzmesh, MatrixXd& GK, VectorXd& GF)
{
	cout << "npt: " << npt << "\n";
	cout << "neq: " << neq << "\n";
	cout << "# of eles: " << bzmesh.size() << "\n";
	clock_t t0 = clock();
#pragma omp parallel for
	for (int e = 0; e<bzmesh.size(); e++)
	{
		if (e != 0 && e % 100 == 0)
		{
			cout << e << " ";
		}
		vector<vector<double>> EK(bzmesh[e].IEN.size(), vector<double>(bzmesh[e].IEN.size(), 0.));
		vector<double> EF(bzmesh[e].IEN.size(), 0.);
		vector<double> Nx(bzmesh[e].IEN.size());
		vector<array<double, 3>> dNdx(bzmesh[e].IEN.size());
		double detJ, Fb;
		array<double, 3> pt;
		for (int i = 0; i<Gpt.size(); i++)
		{
			for (int j = 0; j<Gpt.size(); j++)
			{
				for (int k = 0; k < Gpt.size(); k++)
				{
					BasisFunction_1(Gpt[i], Gpt[j], Gpt[k], bzmesh[e], Nx, dNdx, detJ);
					detJ = wght[i] * wght[j] * wght[k] * detJ;
					ElementMatrix_1(dNdx, detJ, EK);

					bzmesh[e].Para2Phys(Gpt[i], Gpt[j], Gpt[k], pt.data());
					Fb = f_source(pt);
					ElementForce_1(Nx, detJ, Fb, EF);
				}
			}
		}
		Assembly_Dense(EK, EF, bzmesh[e].IEN, GK, GF);
	}
	clock_t t1 = clock();
	double elaps = double(t1 - t0) / CLOCKS_PER_SEC;
	cout << "\ntime of basis: " << elaps << "\n";
	//getchar();
}

void Laplace::OutputMatrix(const SparseMatrix<double>& GK, string fn)
{
	string fname = fn + "_mat.txt";
	ofstream fout;
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << GK << "\n";
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	string fn1 = fn + "_ratio.txt";
	fout.open(fn1.c_str());
	if (fout.is_open())
	{
		double ratio = double(GK.nonZeros()) / double(neq*neq);
		fout << neq << " " << GK.nonZeros() << "\n";
		fout << neq*neq << " " << ratio << "\n";
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fn1 << "!\n";
	}

	//cout << "# of nonzeros: " << GK.nonZeros() << "\n";
}

void Laplace::FindRange(const vector<BezierElement3D>& bzmesh)
{
	double tmp[3][2] = { { 1.e5, -1.e5 }, { 1.e5, -1.e5 }, { 1.e5, -1.e5 } };
	for (uint i = 0; i < bzmesh.size(); i++)
	{
		for (int j = 0; j < 64; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				if (bzmesh[i].pts[j][k] < tmp[k][0]) tmp[k][0] = bzmesh[i].pts[j][k];
				if (bzmesh[i].pts[j][k] > tmp[k][1]) tmp[k][1] = bzmesh[i].pts[j][k];
			}
		}
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			range[i][j] = tmp[i][j];
		}
	}
}









void Laplace::BuildLinearSystem_Coupling(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF)
{
	InitializeSparseMatrix(bzmesh, GK);
	//find max h
	double hmax(0.), dst;
	int bzop[4][2] = { { 0, 63 }, { 3, 60 }, { 15, 48 }, { 12, 51 } };
	for (int e = 0; e < bzmesh.size(); e++)
	{
		for (int i = 0; i<4; i++)
		{
			dst = (bzmesh[e].pts[bzop[i][0]][0] - bzmesh[e].pts[bzop[i][1]][0])*(bzmesh[e].pts[bzop[i][0]][0] - bzmesh[e].pts[bzop[i][1]][0]) +
				(bzmesh[e].pts[bzop[i][0]][1] - bzmesh[e].pts[bzop[i][1]][1])*(bzmesh[e].pts[bzop[i][0]][1] - bzmesh[e].pts[bzop[i][1]][1]) +
				(bzmesh[e].pts[bzop[i][0]][2] - bzmesh[e].pts[bzop[i][1]][2])*(bzmesh[e].pts[bzop[i][0]][2] - bzmesh[e].pts[bzop[i][1]][2]);
			if (dst>hmax) hmax = dst;
		}
		
	}
	hmax = sqrt(hmax);
	double beta(1.e4 / hmax);
	cout << "hmax: " << hmax << "\n";
	hemax = hmax;

	//assembly
	//int nglev(3);
	int nic(pow(2., double(nglev)));//# of integration cells in each direction
	double clen(pow(2., double(-nglev)));
	double clen3(clen*clen*clen);
	cout << "# of eles: " << bzmesh.size() << "\n";
#pragma omp parallel for
	for (int e = 0; e<bzmesh.size(); e++)
	{
		if (e != 0 && e % 100 == 0) cout << e << " ";
		//cout << e << " ";
		vector<vector<double>> EK(bzmesh[e].IEN.size(), vector<double>(bzmesh[e].IEN.size(), 0.));//IEN is different for different types of elements
		vector<double> Nx(bzmesh[e].IEN.size()), EF(bzmesh[e].IEN.size(), 0.);
		vector<array<double, 3>> dNdx(bzmesh[e].IEN.size());
		double detJ, Fb;
		array<double, 3> pt;
		if (bzmesh[e].type == 0)
		{
			for (int i = 0; i<Gpt.size(); i++)
			{
				for (int j = 0; j<Gpt.size(); j++)
				{
					for (int k = 0; k < Gpt.size(); k++)
					{
						BasisFunction_IGA(Gpt[i], Gpt[j], Gpt[k], bzmesh[e], Nx, dNdx, detJ);
						detJ = wght[i] * wght[j] * wght[k] * detJ;
						ElementMatrix_1(dNdx, detJ, EK);

						bzmesh[e].Para2Phys(Gpt[i], Gpt[j], Gpt[k], pt.data());
						Fb = f_source(pt);
						ElementForce_1(Nx, detJ, Fb, EF);
					}
				}
			}
			//boundary condition
			//NitscheBoundary_IGA(bzmesh[e], beta, EK, EF);
		}
		else//irregular and boundary
		{
			double gps[3] = { 0., 0., 0. };
			int icu, icv, icw, i, j, k;
			for (icw = 0; icw < nic; icw++)
			{
				for (icv = 0; icv < nic; icv++)
				{
					for (icu = 0; icu < nic; icu++)
					{
						for (k = 0; k<Gpt.size(); k++)
						{
							for (j = 0; j<Gpt.size(); j++)
							{
								for (i = 0; i < Gpt.size(); i++)
								{
									gps[0] = (Gpt[i] + double(icu))*clen;
									gps[1] = (Gpt[j] + double(icv))*clen;
									gps[2] = (Gpt[k] + double(icw))*clen;
									BasisFunction_RKPM(gps[0], gps[1], gps[2], bzmesh[e], Nx, dNdx, detJ, pt);//only this function is different
									detJ = wght[i] * wght[j] * wght[k] * detJ * clen3;
									ElementMatrix_RKPM(bzmesh[e].cmat.size(), dNdx, detJ, EK);

									//bzmesh[e].Para2Phys(Gpt[i], Gpt[j], Gpt[k], pt.data());
									Fb = f_source(pt);
									ElementForce_RKPM(bzmesh[e].cmat.size(), Nx, detJ, Fb, EF);
								}
							}
						}
					}
				}
			}
			//Nitsche's method coupling interface
			//NitscheInterface(bzmesh[e], beta, EK, EF);

			//boundary condition
			NitscheBoundary_RKPM(bzmesh[e], beta, EK, EF);
		}

		Assembly_Coupling(EK, EF, bzmesh[e].IEN, GK, GF);
	}
	cout << "Done assembly!\n";

	//ofstream fout;
	////fout.open("../io/ptest3/mat.txt");
	//fout.open("../io/ptest3/gf.txt");
	//if (fout.is_open())
	//{
	//	fout << GF << "\n";
	//	fout.close();
	//}
	//cout << "done output matrix!\n";
	//getchar();
}

void Laplace::InitializeSparseMatrix(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK)
{
	cout << "Initialize sparse matrix...\n";
	cout << "# of eles: " << bzmesh.size() << "\n";
	int nresv(6 * 64);
	if (nresv < neq)
	{
		GK.reserve(VectorXi::Constant(neq, nresv));
	}
	else
	{
		GK.reserve(VectorXi::Constant(neq, neq));
	}
	for (int e = 0; e < bzmesh.size(); e++)
	{
		if (e != 0 && e % 500 == 0)
		{
			cout << e << " ";
		}
		for (int i = 0; i < bzmesh[e].IEN.size(); i++)
		{
			for (int j = 0; j < bzmesh[e].IEN.size(); j++)
			{
				int A(bzmesh[e].IEN[i]), B(bzmesh[e].IEN[j]);
				if (A != -1 && B != -1)
				{
					if (IDBC[A] != -1 && IDBC[B] != -1)
					{
						GK.coeffRef(IDBC[A], IDBC[B]) = 0.;
					}
				}
			}
		}
	}
	GK.makeCompressed();
	cout << "done initializing\n";
	cout << "npt: " << npt << "\n";
	cout << "neq: " << neq << "\n";
}

void Laplace::BasisFunction_IGA(double u, double v, double w, const BezierElement3D& bzel, vector<double>& Nx, vector<array<double, 3>>& dNdx, double& detJ)
{
	double Nx0[64];
	double dNdx0[64][3];
	BasisFunction(u, v, w, bzel.pts, Nx0, dNdx0, detJ);
	//#pragma omp parallel for
	for (int i = 0; i < bzel.cmat.size(); i++)
	{
		Nx[i] = 0.;
		dNdx[i][0] = 0.; dNdx[i][1] = 0.; dNdx[i][2] = 0.;
		for (int j = 0; j < 64; j++)
		{
			if (bzel.cmat[i][j] != 0.)
			{
				Nx[i] += bzel.cmat[i][j] * Nx0[j];
				dNdx[i][0] += bzel.cmat[i][j] * dNdx0[j][0];
				dNdx[i][1] += bzel.cmat[i][j] * dNdx0[j][1];
				dNdx[i][2] += bzel.cmat[i][j] * dNdx0[j][2];
			}
		}
	}
}

void Laplace::BasisFunction_RKPM(double u, double v, double w, const BezierElement3D& bzel, vector<double>& Nx, vector<array<double, 3>>& dNdx, double& detJ, array<double,3>& x)
{
	vector<double> Nxg(bzel.cmat.size());//two many operations
	vector<array<double, 3>> dNdxg(bzel.cmat.size());//too many operations
	BasisFunction_IGA(u, v, w, bzel, Nxg, dNdxg, detJ);

	bzel.Para2Phys(u, v, w, x.data());
	vector<double> f;
	vector<array<double, 3>> df;
	RKPM(x.data(), bzel, f, df);

	for (uint i = 0; i < Nxg.size(); i++)
	{
		Nx[i] = Nxg[i];
		dNdx[i][0] = dNdxg[i][0]; dNdx[i][1] = dNdxg[i][1]; dNdx[i][2] = dNdxg[i][2];
	}
	for (uint i = 0; i < f.size(); i++)
	{
		Nx[Nxg.size() + i] = f[i];
		dNdx[Nxg.size() + i][0] = df[i][0];
		dNdx[Nxg.size() + i][1] = df[i][1];
		dNdx[Nxg.size() + i][2] = df[i][2];
	}
}

void Laplace::RKPM(double x[3], const BezierElement3D& bzel, vector<double>& phi, vector<array<double, 3>>& dphi)
{
	double x0[2] = { 0., 0. };
	int deg(3);
	vector<double> ply, ply0;
	vector<array<double, 3>> dp, dp0;
	//double kf, dkf[2];
	vector<double> kf(bzel.mp.size());
	vector<array<double, 3>> dkf(bzel.mp.size());
	Monomial(x0, deg, ply0, dp0);
	MatrixXd m = MatrixXd::Zero(ply0.size(), ply0.size());
	MatrixXd dmdx = MatrixXd::Zero(ply0.size(), ply0.size());
	MatrixXd dmdy = MatrixXd::Zero(ply0.size(), ply0.size());
	MatrixXd dmdz = MatrixXd::Zero(ply0.size(), ply0.size());
	//int nnz(0);
	//cout << "before loop\n";
	for (uint i = 0; i < bzel.mp.size(); i++)
	{
		double x1[3] = { x[0] - bzel.mp[i][0], x[1] - bzel.mp[i][1], x[2] - bzel.mp[i][2] };
		Monomial(x1, deg, ply, dp);
		KernelFunction(x1, bzel.ra[i], kf[i], dkf[i].data());
		MMatrix(ply, dp, kf[i], dkf[i].data(), m, dmdx, dmdy, dmdz);
		//if (kf[i] != 0.) nnz++;
	}
	//cout << "after loop\n";
	//cout << "nnz: " << nnz << "\n";
	//cout << m << "\n";
	//cout << m.determinant() << "\n";
	//getchar();

	PartialPivLU<MatrixXd> dec(m);
	VectorXd p0 = VectorXd::Zero(ply0.size()); p0(0) = 1.;
	VectorXd mi_p0 = dec.solve(p0);

	VectorXd dmix_p0_tmp = dmdx*mi_p0;
	VectorXd dmix_p0 = dec.solve(dmix_p0_tmp);
	VectorXd dmiy_p0_tmp = dmdy*mi_p0;
	VectorXd dmiy_p0 = dec.solve(dmiy_p0_tmp);
	VectorXd dmiz_p0_tmp = dmdz*mi_p0;
	VectorXd dmiz_p0 = dec.solve(dmiz_p0_tmp);

	//MatrixXd mi = m.inverse();
	//MatrixXd mtmp1 = -mi*dmdx;
	//MatrixXd dmidx = mtmp1*mi;
	//MatrixXd mtmp2 = -mi*dmdy;
	//MatrixXd dmidy = mtmp2*mi;
	//MatrixXd mtmp3 = -mi*dmdz;
	//MatrixXd dmidz = mtmp3*mi;
	//cout << "mi:\n";
	//cout << mi << "\n";
	//cout << "dmidx:\n";
	//cout << dmidx << "\n";
	//cout << "dmidy:\n";
	//cout << dmidy << "\n";
	//getchar();

	phi.clear();
	dphi.clear();
	phi.resize(bzel.mp.size(), 0.);
	dphi.resize(bzel.mp.size());
	for (uint pid = 0; pid < bzel.mp.size(); pid++)
	{
		double x1[3] = { x[0] - bzel.mp[pid][0], x[1] - bzel.mp[pid][1], x[2] - bzel.mp[pid][2] };
		Monomial(x1, deg, ply, dp);
		//KernelFunction(x1, bzel.ra[pid], kf[pid], dkf[pid].data());//don't need to be called again

		phi[pid] = 0.;
		for (uint i = 0; i < ply.size(); i++)
		{
			phi[pid] += ply[i] * mi_p0(i);
		}
		//not yet multiplied by kf

		//dphi[pid][0] = 0.; dphi[pid][1] = 0.; dphi[pid][2] = 0.;
		double term[2] = { 0., 0. };
		for (uint i = 0; i < ply.size(); i++)
		{
			term[0] += dp[i][0] * mi_p0(i);
			term[1] += ply[i] * dmix_p0(i);
		}
		dphi[pid][0] = (term[0] - term[1])*kf[pid] + phi[pid] * dkf[pid][0];
		term[0] = 0.; term[1] = 0.;
		for (uint i = 0; i < ply.size(); i++)
		{
			term[0] += dp[i][1] * mi_p0(i);
			term[1] += ply[i] * dmiy_p0(i);
		}
		dphi[pid][1] = (term[0] - term[1])*kf[pid] + phi[pid] * dkf[pid][1];
		term[0] = 0.; term[1] = 0.;
		for (uint i = 0; i < ply.size(); i++)
		{
			term[0] += dp[i][2] * mi_p0(i);
			term[1] += ply[i] * dmiz_p0(i);
		}
		dphi[pid][2] = (term[0] - term[1])*kf[pid] + phi[pid] * dkf[pid][2];

		phi[pid] *= kf[pid];

		//vector<double> pa(ply.size());
		//for (uint i = 0; i < ply.size(); i++)
		//{
		//	pa[i] = ply0[i];
		//}
		//vector<double> tmp1(ply.size());
		//phi[pid] = 0.;
		//for (uint i = 0; i < ply.size(); i++)
		//{
		//	tmp1[i] = 0.;
		//	for (uint j = 0; j < ply.size(); j++)
		//	{
		//		tmp1[i] += mi(i, j)*pa[j];
		//	}
		//	phi[pid] += ply[i] * tmp1[i];
		//}
		////phi *= kf;

		//dphi[pid][0] = 0.; dphi[pid][1] = 0.; dphi[pid][2] = 0.;
		//double tmp2, r1(0.), r2(0.);
		//for (uint i = 0; i < ply.size(); i++)
		//{
		//	tmp2 = 0.;
		//	for (uint j = 0; j < ply.size(); j++)
		//	{
		//		tmp2 += dmidx(i, j)*pa[j];
		//	}
		//	r1 += dp[i][0] * tmp1[i];
		//	r2 += ply[i] * tmp2;
		//}
		//dphi[pid][0] = (r1 + r2)*kf[pid] + phi[pid] * dkf[pid][0];
		//r1 = 0.; r2 = 0.;
		//for (uint i = 0; i < ply.size(); i++)
		//{
		//	tmp2 = 0.;
		//	for (uint j = 0; j < ply.size(); j++)
		//	{
		//		tmp2 += dmidy(i, j)*pa[j];
		//	}
		//	r1 += dp[i][1] * tmp1[i];
		//	r2 += ply[i] * tmp2;
		//}
		//dphi[pid][1] = (r1 + r2)*kf[pid] + phi[pid] * dkf[pid][1];
		//r1 = 0.; r2 = 0.;
		//for (uint i = 0; i < ply.size(); i++)
		//{
		//	tmp2 = 0.;
		//	for (uint j = 0; j < ply.size(); j++)
		//	{
		//		tmp2 += dmidz(i, j)*pa[j];
		//	}
		//	r1 += dp[i][2] * tmp1[i];
		//	r2 += ply[i] * tmp2;
		//}
		//dphi[pid][2] = (r1 + r2)*kf[pid] + phi[pid] * dkf[pid][2];
		//phi[pid] *= kf[pid];
	}
}

void Laplace::Monomial(double x[3], int deg, vector<double>& ply, vector<array<double, 3>>& dp)
{
	ply.clear();
	dp.clear();
	if (deg == 1)
	{
		ply.resize(4);
		dp.resize(4);
		ply[0] = 1.; dp[0][0] = 0.; dp[0][1] = 0.; dp[0][2] = 0.;

		ply[1] = x[0]; ply[2] = x[1]; ply[3] = x[2];
		dp[1][0] = 1.; dp[1][1] = 0.; dp[1][2] = 0.;
		dp[2][0] = 0.; dp[2][1] = 1.; dp[2][2] = 0.;
		dp[3][0] = 0.; dp[3][1] = 0.; dp[3][2] = 1.;
	}
	else if (deg == 2)
	{
		ply.resize(10);
		dp.resize(10);
		ply[0] = 1.; dp[0][0] = 0.; dp[0][1] = 0.; dp[0][2] = 0.;

		ply[1] = x[0]; ply[2] = x[1]; ply[3] = x[2];
		ply[4] = x[0] * x[0]; ply[5] = x[1] * x[1]; ply[6] = x[2] * x[2];
		ply[7] = x[0] * x[1]; ply[8] = x[0] * x[2]; ply[9] = x[1] * x[2];

		dp[1][0] = 1.; dp[1][1] = 0.; dp[1][2] = 0.;  
		dp[2][0] = 0.; dp[2][1] = 1.; dp[2][2] = 0.;
		dp[3][0] = 0.; dp[3][1] = 0.; dp[3][2] = 1.;
		dp[4][0] = 2.*x[0]; dp[4][1] = 0.; dp[4][2] = 0.;
		dp[5][0] = 0.; dp[5][1] = 2.*x[1]; dp[5][2] = 0.;
		dp[6][0] = 0.; dp[6][1] = 0.; dp[6][2] = 2.*x[2];
		dp[7][0] = x[1]; dp[7][1] = x[0]; dp[7][2] = 0.;
		dp[8][0] = x[2]; dp[8][1] = 0.; dp[8][2] = x[0];
		dp[9][0] = 0.; dp[9][1] = x[2]; dp[9][2] = x[1];
	}
	else if (deg == 3)
	{
		ply.resize(20);
		dp.resize(20);
		ply[0] = 1.; dp[0][0] = 0.; dp[0][1] = 0.; dp[0][2] = 0.;

		ply[1] = x[0]; ply[2] = x[1]; ply[3] = x[2];
		ply[4] = x[0] * x[0]; ply[5] = x[1] * x[1]; ply[6] = x[2] * x[2];
		ply[7] = x[0] * x[1]; ply[8] = x[0] * x[2]; ply[9] = x[1] * x[2];
		ply[10] = x[0] * x[0] * x[0]; ply[11] = x[1] * x[1] * x[1]; ply[12] = x[2] * x[2] * x[2];
		ply[13] = x[0] * x[0] * x[1]; ply[14] = x[0] * x[0] * x[2]; 
		ply[15] = x[1] * x[1] * x[0]; ply[16] = x[1] * x[1] * x[2];
		ply[17] = x[2] * x[2] * x[0]; ply[18] = x[2] * x[2] * x[1];
		ply[19] = x[0] * x[1] * x[2];

		dp[1][0] = 1.; dp[1][1] = 0.; dp[1][2] = 0.;
		dp[2][0] = 0.; dp[2][1] = 1.; dp[2][2] = 0.;
		dp[3][0] = 0.; dp[3][1] = 0.; dp[3][2] = 1.;
		dp[4][0] = 2.*x[0]; dp[4][1] = 0.; dp[4][2] = 0.;
		dp[5][0] = 0.; dp[5][1] = 2.*x[1]; dp[5][2] = 0.;
		dp[6][0] = 0.; dp[6][1] = 0.; dp[6][2] = 2.*x[2];
		dp[7][0] = x[1]; dp[7][1] = x[0]; dp[7][2] = 0.;
		dp[8][0] = x[2]; dp[8][1] = 0.; dp[8][2] = x[0];
		dp[9][0] = 0.; dp[9][1] = x[2]; dp[9][2] = x[1];
		dp[10][0] = 3.*x[0] * x[0]; dp[10][1] = 0.; dp[10][2] = 0.;
		dp[11][0] = 0.; dp[11][1] = 3.*x[1] * x[1]; dp[11][2] = 0.;
		dp[12][0] = 0.; dp[12][1] = 0.; dp[12][2] = 3.*x[2] * x[2];
		dp[13][0] = 2.*x[0] * x[1]; dp[13][1] = x[0] * x[0]; dp[13][2] = 0.;
		dp[14][0] = 2.*x[0] * x[2]; dp[14][1] = 0.; dp[14][2] = x[0] * x[0];
		dp[15][0] = x[1] * x[1]; dp[15][1] = 2.*x[0] * x[1]; dp[15][2] = 0.;
		dp[16][0] = 0.; dp[16][1] = 2.*x[1] * x[2]; dp[16][2] = x[1] * x[1];
		dp[17][0] = x[2] * x[2]; dp[17][1] = 0.; dp[17][2] = 2.*x[0] * x[2];
		dp[18][0] = 0.; dp[18][1] = x[2] * x[2]; dp[18][2] = 2.*x[1] * x[2];
		dp[19][0] = x[1] * x[2]; dp[19][1] = x[0] * x[2]; dp[19][2] = x[0] * x[1];
	}
	/*else if (deg == 4)
	{
		ply[1] = x[0]; ply[2] = x[1];
		ply[3] = x[0] * x[0]; ply[4] = x[0] * x[1]; ply[5] = x[1] * x[1];
		ply[6] = x[0] * x[0] * x[0]; ply[7] = x[0] * x[0] * x[1]; ply[8] = x[0] * x[1] * x[1]; ply[9] = x[1] * x[1] * x[1];
		ply[10] = x[0] * x[0] * x[0] * x[0]; ply[11] = x[0] * x[0] * x[0] * x[1]; ply[12] = x[0] * x[0] * x[1] * x[1]; ply[13] = x[0] * x[1] * x[1] * x[1]; ply[14] = x[1] * x[1] * x[1] * x[1];
		dp[1][0] = 1.; dp[2][0] = 0.;
		dp[1][1] = 0.; dp[2][1] = 1.;
		dp[3][0] = 2.*x[0]; dp[4][0] = x[1]; dp[5][0] = 0.;
		dp[3][1] = 0.; dp[4][1] = x[0]; dp[5][1] = 2.*x[1];
		dp[6][0] = 3.*x[0] * x[0]; dp[7][0] = 2.*x[0] * x[1]; dp[8][0] = x[1] * x[1]; dp[9][0] = 0.;
		dp[6][1] = 0.; dp[7][1] = x[0] * x[0]; dp[8][1] = 2.*x[0] * x[1]; dp[9][1] = 3.*x[1] * x[1];
		dp[10][0] = 4.*x[0] * x[0] * x[0]; dp[11][0] = 3.*x[0] * x[0] * x[1]; dp[12][0] = 2.*x[0] * x[1] * x[1]; dp[13][0] = x[1] * x[1] * x[1]; dp[14][0] = 0.;
		dp[10][1] = 0.; dp[11][1] = x[0] * x[0] * x[0]; dp[12][1] = 2.*x[0] * x[0] * x[1]; dp[13][1] = 3.*x[0] * x[1] * x[1]; dp[14][1] = 4.*x[1] * x[1] * x[1];
	}
	else if (deg == 11)
	{
		ply.resize(4, 0.);
		dp.resize(4);
		ply[1] = x[0]; ply[2] = x[1]; ply[3] = x[0] * x[1];
		dp[1][0] = 1.; dp[2][0] = 0.; dp[3][0] = x[1];
		dp[1][1] = 0.; dp[2][1] = 1.; dp[3][1] = x[0];
	}*/
	else
	{
		cerr << "Other degree not considered yet!\n";
		getchar();
	}
}

void Laplace::KernelFunction(double x[3], double a, double& kf, double dkf[3])
{
	KernelFunction_2(x, a, kf, dkf);
	//KernelFunction_3(x, a, kf, dkf);
}

void Laplace::KernelFunction_2(double x[3], double a, double& kf, double dkf[3])//quadratic
{
	double s = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]) / a;
	if (s > 1.)
	{
		kf = 0.; dkf[0] = 0.; dkf[1] = 0.;
		return;
	}
	kf = 1. - 6.*s*s + 8.*s*s*s - 3.*s*s*s*s;
	double dfds = -12.*s + 24.*s*s - 12.*s*s*s;
	double dsdx[3] = { x[0] / (a*a*s), x[1] / (a*a*s), x[2] / (a*a*s) };
	if (fabs(s) < 1.e-8)
	{
		dsdx[0] = 0.; dsdx[1] = 0.; dsdx[2] = 0.;
	}
	dkf[0] = dfds*dsdx[0];
	dkf[1] = dfds*dsdx[1];
	dkf[2] = dfds*dsdx[2];
}

void Laplace::MMatrix(const vector<double>& ply, const vector<array<double, 3>>& dp, double kf, double dkf[3], MatrixXd& m, MatrixXd& dmdx, MatrixXd& dmdy, MatrixXd& dmdz)
{
	for (uint i = 0; i < ply.size(); i++)
	{
		for (uint j = 0; j < ply.size(); j++)
		{
			m(i, j) += (ply[i] * ply[j] * kf);
		}
	}
	for (uint i = 0; i < ply.size(); i++)
	{
		for (uint j = 0; j < ply.size(); j++)
		{
			dmdx(i, j) += (dp[i][0] * ply[j] * kf + ply[i] * dp[j][0] * kf + ply[i] * ply[j] * dkf[0]);
		}
	}
	for (uint i = 0; i < ply.size(); i++)
	{
		for (uint j = 0; j < ply.size(); j++)
		{
			dmdy(i, j) += (dp[i][1] * ply[j] * kf + ply[i] * dp[j][1] * kf + ply[i] * ply[j] * dkf[1]);
		}
	}
	for (uint i = 0; i < ply.size(); i++)
	{
		for (uint j = 0; j < ply.size(); j++)
		{
			dmdz(i, j) += (dp[i][2] * ply[j] * kf + ply[i] * dp[j][2] * kf + ply[i] * ply[j] * dkf[2]);
		}
	}
}

void Laplace::ElementMatrix_RKPM(int start, const vector<array<double, 3>>& dNdx, double detJ, vector<vector<double>>& EK)
{
	uint i, j;
	for (i = start; i<dNdx.size(); i++)
	{
		for (j = start; j<dNdx.size(); j++)
		{
			EK[i][j] += (dNdx[i][0] * dNdx[j][0] + dNdx[i][1] * dNdx[j][1] + dNdx[i][2] * dNdx[j][2])*detJ;
		}
	}
}

void Laplace::ElementForce_RKPM(int start, const vector<double>& Nx, double detJ, double Fb, vector<double>& EF)
{
	for (uint i = start; i<Nx.size(); i++)
	{
		EF[i] += Nx[i] * detJ * Fb;
	}
}

void Laplace::NitscheInterface(const BezierElement3D& cell, double beta, vector<vector<double>>& EK, vector<double>& EF)
{
	vector<double> Nx(cell.IEN.size());
	vector<array<double, 3>> dNdx(cell.IEN.size());
	double bzpts[16][3];
	double detJ, Fb;
	double gps[3];
	array<double, 3> x;
	//int nglev(3);
	int nic(pow(2., double(nglev)));//# of integration subcells in each direction
	double clen(pow(2., double(-nglev)));
	double clen2(clen*clen);
	int i, j, ic, jc;
	double nm[3] = { 0., 0., 0. };
	int ptfloc[6][16] = { { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 }, { 0, 1, 2, 3, 16, 17, 18, 19, 32, 33, 34, 35, 48, 49, 50, 51 },
	{ 3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63 }, { 12, 13, 14, 15, 28, 29, 30, 31, 44, 45, 46, 47, 60, 61, 62, 63 },
	{ 0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60 }, { 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63 } };
	if (cell.bc[0] == 2)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[0][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[0][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[0][ip]][2];
		}
		for (ic = 0; ic < nic; ic++)
		{
			for (jc = 0; jc < nic; jc++)
			{
				for (i = 0; i < Gpt.size(); i++)
				{
					for (j = 0; j < Gpt.size(); j++)
					{
						gps[0] = (Gpt[j] + double(jc))*clen;
						gps[1] = (Gpt[i] + double(ic))*clen;
						gps[2] = 0.;
						BasisFunction_RKPM(gps[0], gps[1], gps[2], cell, Nx, dNdx, detJ, x);
						GetFaceInfo(gps[0], gps[1], bzpts, detJ, nm);
						detJ = wght[i] * wght[j] * detJ * clen2;
						nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
						ElementMatrix_RKPM_Interface(cell.cmat.size(), Nx, dNdx, nm, beta, detJ, EK);
					}
				}
			}
		}
	}
	if (cell.bc[1] == 2)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[1][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[1][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[1][ip]][2];
		}
		for (ic = 0; ic < nic; ic++)
		{
			for (jc = 0; jc < nic; jc++)
			{
				for (i = 0; i < Gpt.size(); i++)
				{
					for (j = 0; j < Gpt.size(); j++)
					{
						gps[0] = (Gpt[j] + double(jc))*clen;
						gps[1] = 0.;
						gps[2] = (Gpt[i] + double(ic))*clen;
						BasisFunction_RKPM(gps[0], gps[1], gps[2], cell, Nx, dNdx, detJ, x);
						GetFaceInfo(gps[0], gps[2], bzpts, detJ, nm);
						detJ = wght[i] * wght[j] * detJ * clen2;
						//nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
						ElementMatrix_RKPM_Interface(cell.cmat.size(), Nx, dNdx, nm, beta, detJ, EK);
					}
				}
			}
		}
	}
	if (cell.bc[2] == 2)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[2][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[2][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[2][ip]][2];
		}
		for (ic = 0; ic < nic; ic++)
		{
			for (jc = 0; jc < nic; jc++)
			{
				for (i = 0; i < Gpt.size(); i++)
				{
					for (j = 0; j < Gpt.size(); j++)
					{
						gps[0] = 1.;
						gps[1] = (Gpt[j] + double(jc))*clen;
						gps[2] = (Gpt[i] + double(ic))*clen;
						BasisFunction_RKPM(gps[0], gps[1], gps[2], cell, Nx, dNdx, detJ, x);
						GetFaceInfo(gps[1], gps[2], bzpts, detJ, nm);
						detJ = wght[i] * wght[j] * detJ * clen2;
						//nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
						ElementMatrix_RKPM_Interface(cell.cmat.size(), Nx, dNdx, nm, beta, detJ, EK);
					}
				}
			}
		}
	}
	if (cell.bc[3] == 2)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[3][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[3][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[3][ip]][2];
		}
		for (ic = 0; ic < nic; ic++)
		{
			for (jc = 0; jc < nic; jc++)
			{
				for (i = 0; i < Gpt.size(); i++)
				{
					for (j = 0; j < Gpt.size(); j++)
					{
						gps[0] = (Gpt[j] + double(jc))*clen;
						gps[1] = 1.;
						gps[2] = (Gpt[i] + double(ic))*clen;
						BasisFunction_RKPM(gps[0], gps[1], gps[2], cell, Nx, dNdx, detJ, x);
						GetFaceInfo(gps[0], gps[2], bzpts, detJ, nm);
						detJ = wght[i] * wght[j] * detJ * clen2;
						nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
						ElementMatrix_RKPM_Interface(cell.cmat.size(), Nx, dNdx, nm, beta, detJ, EK);
					}
				}
			}
		}
	}
	if (cell.bc[4] == 2)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[4][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[4][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[4][ip]][2];
		}
		for (ic = 0; ic < nic; ic++)
		{
			for (jc = 0; jc < nic; jc++)
			{
				for (i = 0; i < Gpt.size(); i++)
				{
					for (j = 0; j < Gpt.size(); j++)
					{
						gps[0] = 0.;
						gps[1] = (Gpt[j] + double(jc))*clen;
						gps[2] = (Gpt[i] + double(ic))*clen;
						BasisFunction_RKPM(gps[0], gps[1], gps[2], cell, Nx, dNdx, detJ, x);
						GetFaceInfo(gps[1], gps[2], bzpts, detJ, nm);
						detJ = wght[i] * wght[j] * detJ * clen2;
						nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
						ElementMatrix_RKPM_Interface(cell.cmat.size(), Nx, dNdx, nm, beta, detJ, EK);
					}
				}
			}
		}
	}
	if (cell.bc[5] == 2)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[5][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[5][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[5][ip]][2];
		}
		for (ic = 0; ic < nic; ic++)
		{
			for (jc = 0; jc < nic; jc++)
			{
				for (i = 0; i < Gpt.size(); i++)
				{
					for (j = 0; j < Gpt.size(); j++)
					{
						gps[0] = (Gpt[j] + double(jc))*clen;
						gps[1] = (Gpt[i] + double(ic))*clen;
						gps[2] = 1.;
						BasisFunction_RKPM(gps[0], gps[1], gps[2], cell, Nx, dNdx, detJ, x);
						GetFaceInfo(gps[0], gps[1], bzpts, detJ, nm);
						detJ = wght[i] * wght[j] * detJ * clen2;
						//nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
						ElementMatrix_RKPM_Interface(cell.cmat.size(), Nx, dNdx, nm, beta, detJ, EK);
					}
				}
			}
		}
	}
}

void Laplace::GetFaceInfo(double u, double v, const double pts[16][3], double& detJ, double nm[3])
{
	double Nu[4] = { (1. - u)*(1. - u)*(1. - u), 3.*(1. - u)*(1. - u)*u, 3.*(1. - u)*u*u, u*u*u };
	double Nv[4] = { (1. - v)*(1. - v)*(1. - v), 3.*(1. - v)*(1. - v)*v, 3.*(1. - v)*v*v, v*v*v };
	double dNdu[4] = { -3.*(1. - u)*(1. - u), 3. - 12.*u + 9.*u*u, 3.*(2. - 3.*u)*u, 3.*u*u };
	double dNdv[4] = { -3.*(1. - v)*(1. - v), 3. - 12.*v + 9.*v*v, 3.*(2. - 3.*v)*v, 3.*v*v };
	double Nt[16], dNdt[16][2];
	int i, j, a, b, loc(0);
	for (i = 0; i<4; i++)
	{
		for (j = 0; j<4; j++)
		{
			Nt[loc] = Nu[j] * Nv[i];
			dNdt[loc][0] = dNdu[j] * Nv[i];
			dNdt[loc][1] = Nu[j] * dNdv[i];
			loc++;
		}
	}
	double dxdt[3][2] = { { 0., 0. }, { 0., 0. }, { 0., 0. } };
	loc = 0;
	for (i = 0; i<4; i++)
	{
		for (j = 0; j<4; j++)
		{
			for (a = 0; a<3; a++)
			{
				for (b = 0; b<2; b++)
				{
					dxdt[a][b] += pts[loc][a] * dNdt[loc][b];
				}
			}
			loc++;
		}
	}
	//double g[2][2] = { { 0., 0. }, { 0., 0. } };
	//for (i = 0; i<2; i++)
	//{
	//	for (j = 0; j<2; j++)
	//	{
	//		g[i][j] = dxdt[0][i] * dxdt[0][j] + dxdt[1][i] * dxdt[1][j] + dxdt[2][i] * dxdt[2][j];
	//	}
	//}
	////double gdet = g[0][0] * g[0][0] - g[0][1] * g[1][0];
	////double gi[2][2] = { { g[1][1] / gdet, -g[0][1] / gdet }, { -g[1][0] / gdet, g[0][0] / gdet } };
	//detJ = 0.25*sqrt(g[0][0] * g[1][1] - g[0][1] * g[1][0]);

	nm[0] = dxdt[1][0] * dxdt[2][1] - dxdt[2][0] * dxdt[1][1];
	nm[1] = dxdt[2][0] * dxdt[0][1] - dxdt[0][0] * dxdt[2][1];
	nm[2] = dxdt[0][0] * dxdt[1][1] - dxdt[1][0] * dxdt[0][1];
	double len = sqrt(nm[0] * nm[0] + nm[1] * nm[1] + nm[2] * nm[2]);
	detJ = 0.25*len;
	nm[0] /= len;
	nm[1] /= len;
	nm[2] /= len;
}

void Laplace::ElementMatrix_RKPM_Interface(int mfid, const vector<double>& Nx, const vector<array<double, 3>>& dNdx, double nm[3], double beta, double detJ, vector<vector<double>>& EK)
{
	unsigned int i, j;
	double tmp1(0.), tmp2(0.), tmp3(0.);
	for (i = mfid; i<dNdx.size(); i++)//symmetry term
	{
		for (j = mfid; j<dNdx.size(); j++)
		{
			tmp1 = (dNdx[i][0] * nm[0] + dNdx[i][1] * nm[1] + dNdx[i][2] * nm[2])*Nx[j];
			EK[i][j] += (-tmp1)*detJ;
		}
		for (j = 0; j < mfid; j++)
		{
			tmp1 = (dNdx[i][0] * nm[0] + dNdx[i][1] * nm[1] + dNdx[i][2] * nm[2])*Nx[j];
			EK[i][j] += tmp1*detJ;
		}
	}
	for (i = mfid; i<dNdx.size(); i++)//consistency term
	{
		for (j = mfid; j<dNdx.size(); j++)
		{
			tmp2 = (dNdx[j][0] * nm[0] + dNdx[j][1] * nm[1] + dNdx[j][2] * nm[2])*Nx[i];
			EK[i][j] += (-tmp2)*detJ;
		}
	}
	for (i = 0; i<mfid; i++)//consistency term
	{
		for (j = mfid; j<dNdx.size(); j++)
		{
			tmp2 = (dNdx[j][0] * nm[0] + dNdx[j][1] * nm[1] + dNdx[j][2] * nm[2])*Nx[i];
			EK[i][j] += tmp2*detJ;
		}
	}
	for (i = mfid; i<dNdx.size(); i++)//penalty term
	{
		for (j = mfid; j<dNdx.size(); j++)
		{
			tmp3 = beta * Nx[i] * Nx[j];
			EK[i][j] += tmp3*detJ;
		}
	}
	for (i = mfid; i<dNdx.size(); i++)//penalty term
	{
		for (j = 0; j<mfid; j++)
		{
			tmp3 = beta * Nx[i] * Nx[j];
			EK[i][j] -= tmp3*detJ;
		}
	}
	for (i = 0; i<mfid; i++)//penalty term
	{
		for (j = mfid; j<dNdx.size(); j++)
		{
			tmp3 = beta * Nx[i] * Nx[j];
			EK[i][j] -= tmp3*detJ;
		}
	}
	for (i = 0; i<mfid; i++)//penalty term
	{
		for (j = 0; j<mfid; j++)
		{
			tmp3 = beta * Nx[i] * Nx[j];
			EK[i][j] += tmp3*detJ;
		}
	}
}

void Laplace::NitscheBoundary_RKPM(const BezierElement3D& cell, double beta, vector<vector<double>>& EK, vector<double>& EF)
{
	vector<double> Nx(cell.IEN.size());
	vector<array<double, 3>> dNdx(cell.IEN.size());
	double bzpts[16][3];
	double detJ, Fb;
	double gps[3];
	array<double, 3> x;
	//int nglev(3);
	int nic(pow(2., double(nglev)));//# of integration subcells in each direction
	double clen(pow(2., double(-nglev)));
	double clen2(clen*clen);
	int i, j, ic, jc;
	double nm[3] = { 0., 0., 0. };
	int ptfloc[6][16] = { { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 }, { 0, 1, 2, 3, 16, 17, 18, 19, 32, 33, 34, 35, 48, 49, 50, 51 },
	{ 3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63 }, { 12, 13, 14, 15, 28, 29, 30, 31, 44, 45, 46, 47, 60, 61, 62, 63 },
	{ 0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60 }, { 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63 } };
	if (cell.bc[0] == 1)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[0][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[0][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[0][ip]][2];
		}
		for (ic = 0; ic < nic; ic++)
		{
			for (jc = 0; jc < nic; jc++)
			{
				for (i = 0; i < Gpt.size(); i++)
				{
					for (j = 0; j < Gpt.size(); j++)
					{
						gps[0] = (Gpt[j] + double(jc))*clen;
						gps[1] = (Gpt[i] + double(ic))*clen;
						gps[2] = 0.;
						BasisFunction_RKPM(gps[0], gps[1], gps[2], cell, Nx, dNdx, detJ, x);
						GetFaceInfo(gps[0], gps[1], bzpts, detJ, nm);
						detJ = wght[i] * wght[j] * detJ * clen2;
						nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
						ElementMatrix_RKPM_Boundary(cell.cmat.size(), Nx, dNdx, nm, beta, detJ, EK);
						Fb = exact_sol(x);
						ElementForce_RKPM_Boundary(cell.cmat.size(), Nx, dNdx, nm, beta, detJ, Fb, EF);
					}
				}
			}
		}
	}
	if (cell.bc[1] == 1)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[1][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[1][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[1][ip]][2];
		}
		for (ic = 0; ic < nic; ic++)
		{
			for (jc = 0; jc < nic; jc++)
			{
				for (i = 0; i < Gpt.size(); i++)
				{
					for (j = 0; j < Gpt.size(); j++)
					{
						gps[0] = (Gpt[j] + double(jc))*clen;
						gps[1] = 0.;
						gps[2] = (Gpt[i] + double(ic))*clen;
						BasisFunction_RKPM(gps[0], gps[1], gps[2], cell, Nx, dNdx, detJ, x);
						GetFaceInfo(gps[0], gps[2], bzpts, detJ, nm);
						detJ = wght[i] * wght[j] * detJ * clen2;
						//nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
						ElementMatrix_RKPM_Boundary(cell.cmat.size(), Nx, dNdx, nm, beta, detJ, EK);
						Fb = exact_sol(x);
						ElementForce_RKPM_Boundary(cell.cmat.size(), Nx, dNdx, nm, beta, detJ, Fb, EF);
					}
				}
			}
		}
	}
	if (cell.bc[2] == 1)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[2][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[2][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[2][ip]][2];
		}
		for (ic = 0; ic < nic; ic++)
		{
			for (jc = 0; jc < nic; jc++)
			{
				for (i = 0; i < Gpt.size(); i++)
				{
					for (j = 0; j < Gpt.size(); j++)
					{
						gps[0] = 1.;
						gps[1] = (Gpt[j] + double(jc))*clen;
						gps[2] = (Gpt[i] + double(ic))*clen;
						BasisFunction_RKPM(gps[0], gps[1], gps[2], cell, Nx, dNdx, detJ, x);
						GetFaceInfo(gps[1], gps[2], bzpts, detJ, nm);
						detJ = wght[i] * wght[j] * detJ * clen2;
						//nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
						ElementMatrix_RKPM_Boundary(cell.cmat.size(), Nx, dNdx, nm, beta, detJ, EK);
						Fb = exact_sol(x);
						ElementForce_RKPM_Boundary(cell.cmat.size(), Nx, dNdx, nm, beta, detJ, Fb, EF);
					}
				}
			}
		}
	}
	if (cell.bc[3] == 1)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[3][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[3][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[3][ip]][2];
		}
		for (ic = 0; ic < nic; ic++)
		{
			for (jc = 0; jc < nic; jc++)
			{
				for (i = 0; i < Gpt.size(); i++)
				{
					for (j = 0; j < Gpt.size(); j++)
					{
						gps[0] = (Gpt[j] + double(jc))*clen;
						gps[1] = 1.;
						gps[2] = (Gpt[i] + double(ic))*clen;
						BasisFunction_RKPM(gps[0], gps[1], gps[2], cell, Nx, dNdx, detJ, x);
						GetFaceInfo(gps[0], gps[2], bzpts, detJ, nm);
						detJ = wght[i] * wght[j] * detJ * clen2;
						nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
						ElementMatrix_RKPM_Boundary(cell.cmat.size(), Nx, dNdx, nm, beta, detJ, EK);
						Fb = exact_sol(x);
						ElementForce_RKPM_Boundary(cell.cmat.size(), Nx, dNdx, nm, beta, detJ, Fb, EF);
					}
				}
			}
		}
	}
	if (cell.bc[4] == 1)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[4][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[4][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[4][ip]][2];
		}
		for (ic = 0; ic < nic; ic++)
		{
			for (jc = 0; jc < nic; jc++)
			{
				for (i = 0; i < Gpt.size(); i++)
				{
					for (j = 0; j < Gpt.size(); j++)
					{
						gps[0] = 0.;
						gps[1] = (Gpt[j] + double(jc))*clen;
						gps[2] = (Gpt[i] + double(ic))*clen;
						BasisFunction_RKPM(gps[0], gps[1], gps[2], cell, Nx, dNdx, detJ, x);
						GetFaceInfo(gps[1], gps[2], bzpts, detJ, nm);
						detJ = wght[i] * wght[j] * detJ * clen2;
						nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
						ElementMatrix_RKPM_Boundary(cell.cmat.size(), Nx, dNdx, nm, beta, detJ, EK);
						Fb = exact_sol(x);
						ElementForce_RKPM_Boundary(cell.cmat.size(), Nx, dNdx, nm, beta, detJ, Fb, EF);
					}
				}
			}
		}
	}
	if (cell.bc[5] == 1)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[5][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[5][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[5][ip]][2];
		}
		for (ic = 0; ic < nic; ic++)
		{
			for (jc = 0; jc < nic; jc++)
			{
				for (i = 0; i < Gpt.size(); i++)
				{
					for (j = 0; j < Gpt.size(); j++)
					{
						gps[0] = (Gpt[j] + double(jc))*clen;
						gps[1] = (Gpt[i] + double(ic))*clen;
						gps[2] = 1.;
						BasisFunction_RKPM(gps[0], gps[1], gps[2], cell, Nx, dNdx, detJ, x);
						GetFaceInfo(gps[0], gps[1], bzpts, detJ, nm);
						detJ = wght[i] * wght[j] * detJ * clen2;
						//nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
						ElementMatrix_RKPM_Boundary(cell.cmat.size(), Nx, dNdx, nm, beta, detJ, EK);
						Fb = exact_sol(x);
						ElementForce_RKPM_Boundary(cell.cmat.size(), Nx, dNdx, nm, beta, detJ, Fb, EF);
					}
				}
			}
		}
	}
}

void Laplace::ElementMatrix_RKPM_Boundary(int mfid, const vector<double>& Nx, const vector<array<double, 3>>& dNdx, double nm[3], double beta, double detJ, vector<vector<double>>& EK)
{
	//only meshfree nodes have contribution on boundary since domain is split
	unsigned int i, j;
	double tmp1(0.), tmp2(0.), tmp3(0.);
	for (i = mfid; i<dNdx.size(); i++)
	{
		for (j = mfid; j<dNdx.size(); j++)
		{
			tmp1 = (dNdx[i][0] * nm[0] + dNdx[i][1] * nm[1] + dNdx[i][2] * nm[2])*Nx[j];
			tmp2 = (dNdx[j][0] * nm[0] + dNdx[j][1] * nm[1] + dNdx[j][2] * nm[2])*Nx[i];
			tmp3 = beta * Nx[i] * Nx[j];
			EK[i][j] += (-tmp1 - tmp2 + tmp3)*detJ;
		}
	}
}

void Laplace::ElementForce_RKPM_Boundary(int mfid, const vector<double>& Nx, const vector<array<double, 3>>& dNdx, double nm[3], double beta, double detJ, double Fb, vector<double>& EF)
{
	double tmp1(0.), tmp2(0.);
	//double sum(0.);
	for (unsigned int i = mfid; i<Nx.size(); i++)
	{
		tmp1 = Fb*(dNdx[i][0] * nm[0] + dNdx[i][1] * nm[1] + dNdx[i][2] * nm[2]);
		tmp2 = beta*Fb*Nx[i];
		EF[i] += (-tmp1 + tmp2) * detJ;
		//sum += EF[i];
	}
}

void Laplace::Set_Problem(const vector<int>& IDBC_in)
{
	npt = IDBC_in.size();
	IDBC = IDBC_in;
	neq = 0;
	for (unsigned int i = 0; i < IDBC.size(); i++)
	{
		if (IDBC[i] != -1) neq++;
	}

	//nglev = 1;//level of subcells for integration
}

void Laplace::Run_Coupling(const vector<BezierElement3D>& bzmesh, string fn)
{
	//CheckReproduce(bzmesh, fn);

	GaussInfo(4);
	SparseMatrix<double> GK(neq, neq);
	VectorXd GF(neq);
	GF.setZero();
	cout << "Building linear system...\n";
	BuildLinearSystem_Coupling(bzmesh, GK, GF);
	Solver_Coupling(GK, GF);
	VisualizeDisp_Coupling(bzmesh, fn);
	vector<double> err;
	ErrorEstimate_Coupling(bzmesh, err);//with exact solution
	VisualizeError(bzmesh, err, fn);
	OutputError(err, fn);
	
	cout << "\nDone run Laplace...\n";
}

void Laplace::Solver_Coupling(SparseMatrix<double>& GK, VectorXd& GF)
{
	cout << "Solving linear system...\n";

	clock_t begin = clock();

	//SimplicialLDLT<SparseMatrix<double>> solver;
	//VectorXd sol = solver.compute(GK).solve(GF);

	//SimplicialLLT<SparseMatrix<double>> solver;
	//VectorXd sol = solver.compute(GK).solve(GF);

	//ConjugateGradient<SparseMatrix<double>, Lower | Upper, IdentityPreconditioner> solver;
	ConjugateGradient<SparseMatrix<double>, Lower | Upper> solver;
	//ConjugateGradient<SparseMatrix<double>, Lower | Upper, SimplicialCholesky<SparseMatrix<double>>> solver;//very slow
	//double tol(double(neq)*10.e-16);
	double tol(1.e-8);
	//cout << "tol before: " << solver.tolerance() << "\n";
	solver.setTolerance(tol);
	//cout << "tol after: " << solver.tolerance() << "\n";
	solver.setMaxIterations(10*npt);
	solver.compute(GK);
	cout << "done compute!\n";
	//VectorXd sol = solver.solve(GF);
	//VectorXd x0 = VectorXd::Ones(neq);
	VectorXd sol = solver.solveWithGuess(GF, u0);
	cout << "\n# iterations: " << solver.iterations() << '\n';
	cout << "estimated error: " << solver.error() << '\n';
	//getchar();

	//BiCGSTAB<SparseMatrix<double>, IncompleteLUT<double>> solver;
	//solver.compute(GK);
	//cout << "done compute!\n";
	//VectorXd sol = solver.solve(GF);
	//cout << "\n# iterations: " << solver.iterations() << '\n';
	//cout << "estimated error: " << solver.error() << '\n';

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << "solver time: " << elapsed_secs << "\n";

	uh.resize(npt);
	for (int i = 0; i<npt; i++)
	{
		if (IDBC[i] != -1)
			uh[i] = sol(IDBC[i]);
		else
			uh[i] = 0.;
	}

	//cout << "\nGF: ";
	//for (uint i = 0; i < GF.size(); i++)
	//{
	//	cout << GF[i] << " ";
	//}
	//getchar();
	//cout << "\nuh: ";
	//for (uint i = 0; i < uh.size(); i++)
	//{
	//	cout << uh[i] << " ";
	//}
	//getchar();


	cout << "Done solving!\n";
}

void Laplace::VisualizeDisp_Coupling(const vector<BezierElement3D>& bzmesh, string fn)
{
	vector<array<double, 3>> spt;//sample points
	vector<double> sdisp;
	vector<array<int, 8>> sele;
	vector<double> errL2;
	vector<array<double, 3>> lpt;//visulize parameter lines
	vector<array<int, 2>> led;//line connectivity
	double detJ;

	for (unsigned int e = 0; e<bzmesh.size(); e++)
	{
		int ns(2);
		//if (bzmesh[e].type == 1) ns = 5;
		vector<double> su(ns);
		for (int i = 0; i < ns; i++)
		{
			su[i] = double(i) / (double(ns) - 1.);
		}

		int loc(0);
		int pstart = spt.size();
		for (int a = 0; a<ns; a++)
		{
			for (int b = 0; b<ns; b++)
			{
				for (int c = 0; c < ns; c++)
				{
					double pt1[3];
					double disp;
					//bzmesh[e].Para2Phys(su[c], su[b], su[a], pt1);
					DispCal_Coupling(su[c], su[b], su[a], bzmesh[e], pt1, disp, detJ);
					array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
					spt.push_back(pt);
					sdisp.push_back(disp);
				}
			}
		}
		int nns[2] = { ns*ns*ns, ns*ns };
		for (int a = 0; a<ns - 1; a++)
		{
			for (int b = 0; b<ns - 1; b++)
			{
				for (int c = 0; c < ns - 1; c++)
				{
					array<int, 8> el;
					el[0] = pstart + a*nns[1] + b*ns + c;
					el[1] = pstart + a*nns[1] + b*ns + c + 1;
					el[2] = pstart + a*nns[1] + (b + 1)*ns + c + 1;
					el[3] = pstart + a*nns[1] + (b + 1)*ns + c;
					el[4] = pstart + (a + 1)*nns[1] + b*ns + c;
					el[5] = pstart + (a + 1)*nns[1] + b*ns + c + 1;
					el[6] = pstart + (a + 1)*nns[1] + (b + 1)*ns + c + 1;
					el[7] = pstart + (a + 1)*nns[1] + (b + 1)*ns + c;
					sele.push_back(el);
				}
			}
		}
		//edges
		int lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(su[a], 0., 0., pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(su[a], 1., 0., pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(su[a], 0., 1., pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(su[a], 1., 1., pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(0., su[a], 0., pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(1., su[a], 0., pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(0., su[a], 1., pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(1., su[a], 1., pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(0., 0., su[a], pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(1., 0., su[a], pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(0., 1., su[a], pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(1., 1., su[a], pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
	}

	string fname = fn + "_disp.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i<spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 9 * sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "8 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3]
				<< " " << sele[i][4] << " " << sele[i][5] << " " << sele[i][6] << " " << sele[i][7] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "12\n";
		}
		//fout<<"\nPOINT_DATA "<<sdisp.size()<<"\nVECTORS disp float\n";
		//for(i=0;i<sdisp.size();i++)
		//{
		//	fout << sdisp[i][0] << " " << sdisp[i][1] << " " << sdisp[i][2] << "\n";
		//}
		//fout << "\nPOINT_DATA " << sse.size() << "\nVECTORS strain float\n";
		//for (i = 0; i<sse.size(); i++)
		//{
		//	fout << sse[i][0] << " " << sse[i][1] << " " << sse[i][2] << "\n";
		//}
		//fout<<"POINT_DATA "<<sss.size()<<"\nVECTORS stress float\n";
		//for(i=0;i<sss.size();i++)
		//{
		//	fout<<sss[i][0]<<" "<<sss[i][1]<<" "<<sss[i][2]<<"\n";
		//}
		//fout<<"POINT_DATA "<<sdisp_err.size()<<"\nVECTORS disp_err float\n";
		//for(i=0;i<sdisp_err.size();i++)
		//{
		//	fout<<sdisp_err[i][0]<<" "<<sdisp_err[i][1]<<" 0\n";
		//}

		fout << "POINT_DATA " << sdisp.size() << "\nSCALARS err float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<sdisp.size(); i++)
		{
			fout << sdisp[i] << "\n";
		}

		//fout<<"\nCELL_DATA "<<errL2.size()<<"\nSCALARS Error float 1\nLOOKUP_TABLE default\n";
		//for(i=0;i<errL2.size();i++)
		//{
		//	fout<<errL2[i]<<"\n";
		//	//fout<<eles[i].type<<"\n";
		//}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	string fname1(fn + "-lines.vtk");
	ofstream fout1;
	fout1.open(fname1.c_str());
	if (fout1.is_open())
	{
		fout1 << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout1 << "POINTS " << lpt.size() << " float\n";
		for (i = 0; i<lpt.size(); i++)
		{
			fout1 << lpt[i][0] << " " << lpt[i][1] << " " << lpt[i][2] << "\n";
		}
		fout1 << "\nCELLS " << led.size() << " " << 3 * led.size() << '\n';
		for (i = 0; i<led.size(); i++)
		{
			fout1 << "2 " << led[i][0] << " " << led[i][1] << '\n';
		}
		fout1 << "\nCELL_TYPES " << led.size() << '\n';
		for (i = 0; i<led.size(); i++)
		{
			fout1 << "3\n";
		}
		fout1.close();
	}
	else
	{
		cout << "Cannot open " << fname1 << "!\n";
	}
}

void Laplace::DispCal_Coupling(double u, double v, double w, const BezierElement3D& bzel, double pt[3], double& disp, double& detJ)
{
	//double detJ;
	vector<double> Nx(bzel.IEN.size());
	vector<array<double, 3>> dNdx(bzel.IEN.size());

	if (bzel.type == 0)
	{
		double pt0[3];
		bzel.Para2Phys(u, v, w, pt0);
		pt[0] = pt0[0]; pt[1] = pt0[1]; pt[2] = pt0[2];
		BasisFunction_IGA(u, v, w, bzel, Nx, dNdx, detJ);
		disp = 0.;
		for (uint i = 0; i < bzel.IEN.size(); i++)
		{
			disp += Nx[i] * uh[bzel.IEN[i]];
		}
	}
	else
	{
		array<double, 3> x;
		BasisFunction_RKPM(u, v, w, bzel, Nx, dNdx, detJ, x);
		pt[0] = x[0]; pt[1] = x[1]; pt[2] = x[2];
		disp = 0.;
		for (uint i = bzel.cmat.size(); i < bzel.IEN.size(); i++)
		{
			disp += Nx[i] * uh[bzel.IEN[i]];
		}
	}
}

void Laplace::ErrorEstimate_Coupling(const vector<BezierElement3D>& bzmesh, vector<double>& err)
{
	cout << "Error estimating...\n";
	err.clear();
	err.resize(bzmesh.size(), 0.);
#pragma omp parallel for
	for (int e = 0; e < bzmesh.size(); e++)
	{
		double L2, H1;
		ElementError_Coupling(bzmesh[e], L2, H1);
		err[e] = L2;
	}
}

void Laplace::ElementError_Coupling(const BezierElement3D& bzel, double& L2, double& H1)
{
	L2 = 0.; H1 = 0.;
	uint i, j, k;
	double disp, detJ, sol_e;
	array<double, 3> pt;
	int nic(pow(2., double(nglev)));//# of integration cells in each direction
	double clen(pow(2., double(-nglev)));
	double clen3(clen*clen*clen);
	if (bzel.type == 0)
	{
		for (i = 0; i < Gpt.size(); i++)
		{
			for (j = 0; j < Gpt.size(); j++)
			{
				for (k = 0; k < Gpt.size(); k++)
				{
					DispCal_Coupling(Gpt[i], Gpt[j], Gpt[k], bzel, pt.data(), disp, detJ);
					sol_e = exact_sol(pt);
					L2 += (wght[i] * wght[j] * wght[k] * detJ*(disp - sol_e)*(disp - sol_e));
				}
			}
		}
	}
	else
	{
		double gps[3] = { 0., 0., 0. };
		int icu, icv, icw;
		for (icw = 0; icw < nic; icw++)
		{
			for (icv = 0; icv < nic; icv++)
			{
				for (icu = 0; icu < nic; icu++)
				{
					for (k = 0; k<Gpt.size(); k++)
					{
						for (j = 0; j<Gpt.size(); j++)
						{
							for (i = 0; i < Gpt.size(); i++)
							{
								gps[0] = (Gpt[i] + double(icu))*clen;
								gps[1] = (Gpt[j] + double(icv))*clen;
								gps[2] = (Gpt[k] + double(icw))*clen;
								DispCal_Coupling(gps[0], gps[1], gps[2], bzel, pt.data(), disp, detJ);
								sol_e = exact_sol(pt);
								L2 += (wght[i] * wght[j] * wght[k] * detJ * clen3 * (disp - sol_e)*(disp - sol_e));
							}
						}
					}
				}
			}
		}
	}
}

void Laplace::OutputError(const vector<double>& err, string fn)
{
	double L2all(0.);
	for (uint i = 0; i < err.size(); i++)
	{
		L2all += err[i];
	}
	L2all = sqrt(L2all);

	cout << hemax << " " << L2all << "\n";

	ofstream fout;
	string fname(fn + "_err.txt");
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << hemax << " " << L2all << "\n";
		fout.close();
	}
	else
	{
		cerr << "Can't open " << fname << "\n";
	}
}

void Laplace::NitscheBoundary_IGA(const BezierElement3D& cell, double beta, vector<vector<double>>& EK, vector<double>& EF)
{
	vector<double> Nx(cell.IEN.size());
	vector<array<double, 3>> dNdx(cell.IEN.size());
	double bzpts[16][3];
	double detJ, Fb;
	double gps[3];
	array<double, 3> x;
	//int nglev(3);
	int nic(pow(2., double(nglev)));//# of integration subcells in each direction
	double clen(pow(2., double(-nglev)));
	double clen2(clen*clen);
	int i, j, ic, jc;
	double nm[3] = { 0., 0., 0. };
	int ptfloc[6][16] = { { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 }, { 0, 1, 2, 3, 16, 17, 18, 19, 32, 33, 34, 35, 48, 49, 50, 51 },
	{ 3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63 }, { 12, 13, 14, 15, 28, 29, 30, 31, 44, 45, 46, 47, 60, 61, 62, 63 },
	{ 0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60 }, { 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63 } };
	if (cell.bc[0] == 1)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[0][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[0][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[0][ip]][2];
		}
		for (ic = 0; ic < nic; ic++)
		{
			for (jc = 0; jc < nic; jc++)
			{
				for (i = 0; i < Gpt.size(); i++)
				{
					for (j = 0; j < Gpt.size(); j++)
					{
						gps[0] = (Gpt[j] + double(jc))*clen;
						gps[1] = (Gpt[i] + double(ic))*clen;
						gps[2] = 0.;
						BasisFunction_IGA(gps[0], gps[1], gps[2], cell, Nx, dNdx, detJ);
						GetFaceInfo(gps[0], gps[1], bzpts, detJ, nm);
						detJ = wght[i] * wght[j] * detJ * clen2;
						nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
						ElementMatrix_IGA_Boundary(0, Nx, dNdx, nm, beta, detJ, EK);
						cell.Para2Phys(gps[0], gps[1], gps[2], x.data());
						Fb = exact_sol(x);
						ElementForce_IGA_Boundary(0, Nx, dNdx, nm, beta, detJ, Fb, EF);
					}
				}
			}
		}
	}
	if (cell.bc[1] == 1)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[1][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[1][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[1][ip]][2];
		}
		for (ic = 0; ic < nic; ic++)
		{
			for (jc = 0; jc < nic; jc++)
			{
				for (i = 0; i < Gpt.size(); i++)
				{
					for (j = 0; j < Gpt.size(); j++)
					{
						gps[0] = (Gpt[j] + double(jc))*clen;
						gps[1] = 0.;
						gps[2] = (Gpt[i] + double(ic))*clen;
						BasisFunction_IGA(gps[0], gps[1], gps[2], cell, Nx, dNdx, detJ);
						GetFaceInfo(gps[0], gps[2], bzpts, detJ, nm);
						detJ = wght[i] * wght[j] * detJ * clen2;
						//nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
						ElementMatrix_IGA_Boundary(0, Nx, dNdx, nm, beta, detJ, EK);
						cell.Para2Phys(gps[0], gps[1], gps[2], x.data());
						Fb = exact_sol(x);
						ElementForce_IGA_Boundary(0, Nx, dNdx, nm, beta, detJ, Fb, EF);
					}
				}
			}
		}
	}
	if (cell.bc[2] == 1)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[2][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[2][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[2][ip]][2];
		}
		for (ic = 0; ic < nic; ic++)
		{
			for (jc = 0; jc < nic; jc++)
			{
				for (i = 0; i < Gpt.size(); i++)
				{
					for (j = 0; j < Gpt.size(); j++)
					{
						gps[0] = 1.;
						gps[1] = (Gpt[j] + double(jc))*clen;
						gps[2] = (Gpt[i] + double(ic))*clen;
						BasisFunction_IGA(gps[0], gps[1], gps[2], cell, Nx, dNdx, detJ);
						GetFaceInfo(gps[1], gps[2], bzpts, detJ, nm);
						detJ = wght[i] * wght[j] * detJ * clen2;
						//nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
						ElementMatrix_IGA_Boundary(0, Nx, dNdx, nm, beta, detJ, EK);
						cell.Para2Phys(gps[0], gps[1], gps[2], x.data());
						Fb = exact_sol(x);
						ElementForce_IGA_Boundary(0, Nx, dNdx, nm, beta, detJ, Fb, EF);
					}
				}
			}
		}
	}
	if (cell.bc[3] == 1)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[3][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[3][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[3][ip]][2];
		}
		for (ic = 0; ic < nic; ic++)
		{
			for (jc = 0; jc < nic; jc++)
			{
				for (i = 0; i < Gpt.size(); i++)
				{
					for (j = 0; j < Gpt.size(); j++)
					{
						gps[0] = (Gpt[j] + double(jc))*clen;
						gps[1] = 1.;
						gps[2] = (Gpt[i] + double(ic))*clen;
						BasisFunction_IGA(gps[0], gps[1], gps[2], cell, Nx, dNdx, detJ);
						GetFaceInfo(gps[0], gps[2], bzpts, detJ, nm);
						detJ = wght[i] * wght[j] * detJ * clen2;
						nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
						ElementMatrix_IGA_Boundary(0, Nx, dNdx, nm, beta, detJ, EK);
						cell.Para2Phys(gps[0], gps[1], gps[2], x.data());
						Fb = exact_sol(x);
						ElementForce_IGA_Boundary(0, Nx, dNdx, nm, beta, detJ, Fb, EF);
					}
				}
			}
		}
	}
	if (cell.bc[4] == 1)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[4][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[4][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[4][ip]][2];
		}
		for (ic = 0; ic < nic; ic++)
		{
			for (jc = 0; jc < nic; jc++)
			{
				for (i = 0; i < Gpt.size(); i++)
				{
					for (j = 0; j < Gpt.size(); j++)
					{
						gps[0] = 0.;
						gps[1] = (Gpt[j] + double(jc))*clen;
						gps[2] = (Gpt[i] + double(ic))*clen;
						BasisFunction_IGA(gps[0], gps[1], gps[2], cell, Nx, dNdx, detJ);
						GetFaceInfo(gps[1], gps[2], bzpts, detJ, nm);
						detJ = wght[i] * wght[j] * detJ * clen2;
						nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
						ElementMatrix_IGA_Boundary(0, Nx, dNdx, nm, beta, detJ, EK);
						cell.Para2Phys(gps[0], gps[1], gps[2], x.data());
						Fb = exact_sol(x);
						ElementForce_IGA_Boundary(0, Nx, dNdx, nm, beta, detJ, Fb, EF);
					}
				}
			}
		}
	}
	if (cell.bc[5] == 1)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[5][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[5][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[5][ip]][2];
		}
		for (ic = 0; ic < nic; ic++)
		{
			for (jc = 0; jc < nic; jc++)
			{
				for (i = 0; i < Gpt.size(); i++)
				{
					for (j = 0; j < Gpt.size(); j++)
					{
						gps[0] = (Gpt[j] + double(jc))*clen;
						gps[1] = (Gpt[i] + double(ic))*clen;
						gps[2] = 1.;
						BasisFunction_IGA(gps[0], gps[1], gps[2], cell, Nx, dNdx, detJ);
						GetFaceInfo(gps[0], gps[1], bzpts, detJ, nm);
						detJ = wght[i] * wght[j] * detJ * clen2;
						//nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
						ElementMatrix_IGA_Boundary(0, Nx, dNdx, nm, beta, detJ, EK);
						cell.Para2Phys(gps[0], gps[1], gps[2], x.data());
						Fb = exact_sol(x);
						ElementForce_IGA_Boundary(0, Nx, dNdx, nm, beta, detJ, Fb, EF);
					}
				}
			}
		}
	}
}

void Laplace::ElementMatrix_IGA_Boundary(int mfid, const vector<double>& Nx, const vector<array<double, 3>>& dNdx, double nm[3], double beta, double detJ, vector<vector<double>>& EK)
{
	//only meshfree nodes have contribution on boundary since domain is split
	unsigned int i, j;
	double tmp1(0.), tmp2(0.), tmp3(0.);
	for (i = mfid; i<dNdx.size(); i++)
	{
		for (j = mfid; j<dNdx.size(); j++)
		{
			tmp1 = (dNdx[i][0] * nm[0] + dNdx[i][1] * nm[1] + dNdx[i][2] * nm[2])*Nx[j];
			tmp2 = (dNdx[j][0] * nm[0] + dNdx[j][1] * nm[1] + dNdx[j][2] * nm[2])*Nx[i];
			tmp3 = beta * Nx[i] * Nx[j];
			EK[i][j] += (-tmp1 - tmp2 + tmp3)*detJ;
		}
	}
}

void Laplace::ElementForce_IGA_Boundary(int mfid, const vector<double>& Nx, const vector<array<double, 3>>& dNdx, double nm[3], double beta, double detJ, double Fb, vector<double>& EF)
{
	double tmp1(0.), tmp2(0.);
	//double sum(0.);
	for (unsigned int i = mfid; i<Nx.size(); i++)
	{
		tmp1 = Fb*(dNdx[i][0] * nm[0] + dNdx[i][1] * nm[1] + dNdx[i][2] * nm[2]);
		tmp2 = beta*Fb*Nx[i];
		EF[i] += (-tmp1 + tmp2) * detJ;
		//sum += EF[i];
	}
}

void Laplace::CheckReproduce(const vector<BezierElement3D>& bzmesh, string fn)
{
	vector<array<double, 3>> spt;//sample points
	vector<array<double,3>> sdisp;
	vector<array<int, 8>> sele;
	vector<double> errL2;
	vector<array<double, 3>> lpt;//visulize parameter lines
	vector<array<int, 2>> led;//line connectivity
	double detJ;

	for (unsigned int e = 0; e<bzmesh.size(); e++)
	{
		int ns(2);
		//if (bzmesh[e].type == 1) ns = 5;
		vector<double> su(ns);
		for (int i = 0; i < ns; i++)
		{
			su[i] = double(i) / (double(ns) - 1.);
		}

		int loc(0);
		int pstart = spt.size();
		if (bzmesh[e].type != 0)
		{
			vector<double> Nx(bzmesh[e].IEN.size());
			vector<array<double, 3>> dNdx(bzmesh[e].IEN.size());
			array<double, 3> x;
			double detJ;
			for (int a = 0; a<ns; a++)
			{
				for (int b = 0; b<ns; b++)
				{
					for (int c = 0; c < ns; c++)
					{
						array<double, 3> disp = { 0., 0., 0. };
						BasisFunction_RKPM(su[c], su[b], su[a], bzmesh[e], Nx, dNdx, detJ, x);
						//const
						for (uint i = bzmesh[e].cmat.size(); i < Nx.size(); i++)
						{
							disp[0] += dNdx[i][0];
						}
						//disp[0] -= 1.;
						//linear
						for (uint i = bzmesh[e].cmat.size(); i < Nx.size(); i++)
						{
							disp[1] += (dNdx[i][0] * (x[0] - bzmesh[e].mp[i - bzmesh[e].cmat.size()][0]) + Nx[i]);
						}
						//p2
						for (uint i = bzmesh[e].cmat.size(); i < Nx.size(); i++)
						{
							disp[2] += (dNdx[i][0] * (x[0] - bzmesh[e].mp[i - bzmesh[e].cmat.size()][0])* (x[0] - bzmesh[e].mp[i - bzmesh[e].cmat.size()][0]) + 
								Nx[i] * 2. * (x[0] - bzmesh[e].mp[i - bzmesh[e].cmat.size()][0]));
						}

						spt.push_back(x);
						sdisp.push_back(disp);
					}
				}
			}
		}
		else
		{
			for (int a = 0; a<ns; a++)
			{
				for (int b = 0; b<ns; b++)
				{
					for (int c = 0; c < ns; c++)
					{
						array<double, 3> pt;
						array<double, 3> disp = { 0., 0., 0. };
						bzmesh[e].Para2Phys(su[c], su[b], su[a], pt.data());
						spt.push_back(pt);
						sdisp.push_back(disp);
					}
				}
			}
		}
		
		int nns[2] = { ns*ns*ns, ns*ns };
		for (int a = 0; a<ns - 1; a++)
		{
			for (int b = 0; b<ns - 1; b++)
			{
				for (int c = 0; c < ns - 1; c++)
				{
					array<int, 8> el;
					el[0] = pstart + a*nns[1] + b*ns + c;
					el[1] = pstart + a*nns[1] + b*ns + c + 1;
					el[2] = pstart + a*nns[1] + (b + 1)*ns + c + 1;
					el[3] = pstart + a*nns[1] + (b + 1)*ns + c;
					el[4] = pstart + (a + 1)*nns[1] + b*ns + c;
					el[5] = pstart + (a + 1)*nns[1] + b*ns + c + 1;
					el[6] = pstart + (a + 1)*nns[1] + (b + 1)*ns + c + 1;
					el[7] = pstart + (a + 1)*nns[1] + (b + 1)*ns + c;
					sele.push_back(el);
				}
			}
		}
	}

	string fname = fn + "_rep.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i<spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 9 * sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "8 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3]
				<< " " << sele[i][4] << " " << sele[i][5] << " " << sele[i][6] << " " << sele[i][7] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "12\n";
		}
		fout<<"\nPOINT_DATA "<<sdisp.size()<<"\nVECTORS disp float\n";
		for(i=0;i<sdisp.size();i++)
		{
			fout << sdisp[i][0] << " " << sdisp[i][1] << " " << sdisp[i][2] << "\n";
		}
		//fout << "POINT_DATA " << sdisp.size() << "\nSCALARS err float 1\nLOOKUP_TABLE default\n";
		//for (uint i = 0; i<sdisp.size(); i++)
		//{
		//	fout << sdisp[i] << "\n";
		//}
		//fout<<"\nCELL_DATA "<<errL2.size()<<"\nSCALARS Error float 1\nLOOKUP_TABLE default\n";
		//for(i=0;i<errL2.size();i++)
		//{
		//	fout<<errL2[i]<<"\n";
		//	//fout<<eles[i].type<<"\n";
		//}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}








//coupling with Bezier

void Laplace::Set_Problem_Bezier(const vector<int>& IDBC_in, const vector<double>& gh_in)
{
	npt = IDBC_in.size();
	IDBC = IDBC_in;
	gh = gh_in;
	neq = 0;
	for (unsigned int i = 0; i < IDBC.size(); i++)
	{
		if (IDBC[i] != -1) neq++;
	}

	u0.resize(gh_in.size());
	for (unsigned int i = 0; i < gh_in.size(); i++)
	{
		u0(i) = gh_in[i];
	}

	//nglev = 1;//level of subcells for integration
}

void Laplace::Run_Coupling_Bezier(const vector<BezierElement3D>& bzmesh, string fn)
{
	GaussInfo(4);

	SparseMatrix<double> GK(neq, neq);
	VectorXd GF(neq);
	GF.setZero();
	cout << "Building linear system...\n";
	BuildLinearSystem_Coupling_Bezier(bzmesh, GK, GF);

	//Solver_Coupling_Bezier(GK, GF);
	//Solver_MKL(GK, GF);                 //Angran Comments
	//Solver_Matlab(GK, GF);

	//cout << "Write matrix...\n";
	////clock_t begin = clock();
	//Write2MKL_SparseMat(GK, GF, fn);
	////Write2Matlab_SparseMat(GK, GF, fn);
	////clock_t end = clock();
	////double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	////cout << "Writing time: " << elapsed_secs << "\n";

	//ReadSolution(fn);
	//for (int i = 0; i < uh.size(); i++)
	//	cout << i<<" "<<uh[i] << endl;
	VisualizeDisp_Coupling_Bezier(bzmesh, fn);
	vector<vector<double>> err(2);//0 for L2 and 1 for H1
	//ErrorEstimate_Coupling_Bezier(bzmesh, err);//with exact solution
	//VisualizeError(bzmesh, err[0], fn);
	////VisualizeError_BezierCouple(bzmesh, err[0], fn);
	//OutputError_All(err, fn);

	cout << "\nDone run Laplace...\n";
}

void Laplace::BasisFunction_Bezier(double u, double v, double w, const vector<array<double, 3>>& pt, vector<double>& Nx, vector<array<double, 3>>& dNdx, double& detJ)
{
	double Nu[4] = { (1. - u)*(1. - u)*(1. - u), 3.*(1. - u)*(1. - u)*u, 3.*(1. - u)*u*u, u*u*u };
	double Nv[4] = { (1. - v)*(1. - v)*(1. - v), 3.*(1. - v)*(1. - v)*v, 3.*(1. - v)*v*v, v*v*v };
	double Nw[4] = { (1. - w)*(1. - w)*(1. - w), 3.*(1. - w)*(1. - w)*w, 3.*(1. - w)*w*w, w*w*w };
	double dNdu[4] = { -3.*(1. - u)*(1. - u), 3. - 12.*u + 9.*u*u, 3.*(2. - 3.*u)*u, 3.*u*u };
	double dNdv[4] = { -3.*(1. - v)*(1. - v), 3. - 12.*v + 9.*v*v, 3.*(2. - 3.*v)*v, 3.*v*v };
	double dNdw[4] = { -3.*(1. - w)*(1. - w), 3. - 12.*w + 9.*w*w, 3.*(2. - 3.*w)*w, 3.*w*w };
	double dNdt[64][3];
	int i, j, k, a, b, loc(0);
	for (i = 0; i<4; i++)
	{
		for (j = 0; j<4; j++)
		{
			for (k = 0; k < 4; k++)
			{
				Nx[loc] = Nu[k] * Nv[j] * Nw[i];
				dNdt[loc][0] = dNdu[k] * Nv[j] * Nw[i];
				dNdt[loc][1] = Nu[k] * dNdv[j] * Nw[i];
				dNdt[loc][2] = Nu[k] * Nv[j] * dNdw[i];
				loc++;
			}
		}
	}
	Matrix3d dxdt = Matrix3d::Zero();
	loc = 0;
	for (i = 0; i<4; i++)
	{
		for (j = 0; j<4; j++)
		{
			for (k = 0; k < 4; k++)
			{
				for (a = 0; a<3; a++)
				{
					for (b = 0; b<3; b++)
					{
						dxdt(a, b) += pt[loc][a] * dNdt[loc][b];
					}
				}
				loc++;
			}
		}
	}
	//cout << dxdt << "\n\n";
	Matrix3d dtdx = dxdt.inverse();
	for (i = 0; i<64; i++)
	{
		dNdx[i][0] = dNdt[i][0] * dtdx(0, 0) + dNdt[i][1] * dtdx(1, 0) + dNdt[i][2] * dtdx(2, 0);
		dNdx[i][1] = dNdt[i][0] * dtdx(0, 1) + dNdt[i][1] * dtdx(1, 1) + dNdt[i][2] * dtdx(2, 1);
		dNdx[i][2] = dNdt[i][0] * dtdx(0, 2) + dNdt[i][1] * dtdx(1, 2) + dNdt[i][2] * dtdx(2, 2);
	}
	detJ = dxdt.determinant();
	detJ = 0.125*detJ;
}

void Laplace::BasisFunction_BZAll(double u, double v, double w, const BezierElement3D& bzel, vector<double>& Nx, vector<array<double, 3>>& dNdx, double& detJ)
{
	double Nx0[64];
	double dNdx0[64][3];
	BasisFunction(u, v, w, bzel.pts, Nx0, dNdx0, detJ);
	for (int i = 0; i < bzel.cmat.size(); i++)
	{
		Nx[i] = 0.;
		dNdx[i][0] = 0.; dNdx[i][1] = 0.; dNdx[i][2] = 0.;
		for (int j = 0; j < 64; j++)
		{
			Nx[i] += bzel.cmat[i][j] * Nx0[j];
			dNdx[i][0] += bzel.cmat[i][j] * dNdx0[j][0];
			dNdx[i][1] += bzel.cmat[i][j] * dNdx0[j][1];
			dNdx[i][2] += bzel.cmat[i][j] * dNdx0[j][2];
		}
	}
	int loc(0);
	for (int i = 0; i < 64; i++)
	{
		loc = bzel.IEN.size() + i;
		Nx[loc] = Nx0[i];
		dNdx[loc][0] = dNdx0[i][0];
		dNdx[loc][1] = dNdx0[i][1];
		dNdx[loc][2] = dNdx0[i][2];
	}

	//cout << "Nx: ";
	//for (uint i = 0; i < Nx.size(); i++)
	//{
	//	cout << Nx[i] << " ";
	//}
	//cout << "\n";
	//getchar();
}

void Laplace::InitializeSparseMatrix_Bezier(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK)
{
	cout << "Initialize sparse matrix...\n";
	cout << "# of eles: " << bzmesh.size() << "\n";
	vector<Triplet<double>> trilist;
	for (int e = 0; e < bzmesh.size(); e++)
	{
		if (e != 0 && e % 500 == 0)
		{
			cout << e << " ";
		}
		if (bzmesh[e].bzflag == 0)
		{
			for (int i = 0; i < bzmesh[e].IEN.size(); i++)
			{
				for (int j = 0; j < bzmesh[e].IEN.size(); j++)
				{
					int A(bzmesh[e].IEN[i]), B(bzmesh[e].IEN[j]);
					//if (A != -1 && B != -1)
					if (A != -1 && B != -1 && A >= B)
					{
						if (IDBC[A] != -1 && IDBC[B] != -1)
						{
							//GK.coeffRef(IDBC[A], IDBC[B]) = 0.;
							//cout << IDBC[A] << " " << IDBC[B] << " " << trilist.size() << "\n";
							trilist.push_back(Triplet<double>(IDBC[A], IDBC[B], 0.));
						}
					}
				}
			}
		}
		else if (bzmesh[e].bzflag == 1 && bzmesh[e].bzcouple == 0)
		{
			for (int i = 0; i < bzmesh[e].IENb.size(); i++)
			{
				for (int j = 0; j < bzmesh[e].IENb.size(); j++)
				{
					int A(bzmesh[e].IENb[i]), B(bzmesh[e].IENb[j]);
					//if (A != -1 && B != -1)
					if (A != -1 && B != -1 && A >= B)
					{
						if (IDBC[A] != -1 && IDBC[B] != -1)
						{
							//GK.coeffRef(IDBC[A], IDBC[B]) = 0.;
							trilist.push_back(Triplet<double>(IDBC[A], IDBC[B], 0.));
						}
					}
				}
			}
		}
		else if (bzmesh[e].bzflag == 1 && bzmesh[e].bzcouple == 1)
		{
			vector<int> IEN_all(bzmesh[e].IEN.size() + bzmesh[e].IENb.size());
			for (uint i = 0; i < bzmesh[e].IEN.size(); i++) IEN_all[i] = bzmesh[e].IEN[i];
			for (uint i = 0; i < bzmesh[e].IENb.size(); i++) IEN_all[i + bzmesh[e].IEN.size()] = bzmesh[e].IENb[i];
			for (int i = 0; i < IEN_all.size(); i++)
			{
				for (int j = 0; j < IEN_all.size(); j++)
				{
					int A(IEN_all[i]), B(IEN_all[j]);
					//if (A != -1 && B != -1)
					if (A != -1 && B != -1 && A >= B)
					{
						if (IDBC[A] != -1 && IDBC[B] != -1)
						{
							//GK.coeffRef(IDBC[A], IDBC[B]) = 0.;
							trilist.push_back(Triplet<double>(IDBC[A], IDBC[B], 0.));
						}
					}
				}
			}
		}
		//cout << trilist.size() << "\n";
	}
	GK.setFromTriplets(trilist.begin(), trilist.end());
	GK.makeCompressed();
	cout << "done initializing\n";
	cout << "npt: " << npt << "\n";
	cout << "neq: " << neq << "\n";
}

void Laplace::BuildLinearSystem_Coupling_Bezier(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF)
{
	InitializeSparseMatrix_Bezier(bzmesh, GK);//matrix need to be large enough, otherwise parallel part got problem
	//find max h
	double hmax(0.), dst;
	int bzop[4][2] = { { 0, 63 }, { 3, 60 }, { 15, 48 }, { 12, 51 } };
	for (int e = 0; e < bzmesh.size(); e++)
	{
		for (int i = 0; i<4; i++)
		{
			dst = (bzmesh[e].pts[bzop[i][0]][0] - bzmesh[e].pts[bzop[i][1]][0])*(bzmesh[e].pts[bzop[i][0]][0] - bzmesh[e].pts[bzop[i][1]][0]) +
				(bzmesh[e].pts[bzop[i][0]][1] - bzmesh[e].pts[bzop[i][1]][1])*(bzmesh[e].pts[bzop[i][0]][1] - bzmesh[e].pts[bzop[i][1]][1]) +
				(bzmesh[e].pts[bzop[i][0]][2] - bzmesh[e].pts[bzop[i][1]][2])*(bzmesh[e].pts[bzop[i][0]][2] - bzmesh[e].pts[bzop[i][1]][2]);
			if (dst>hmax) hmax = dst;
		}

	}
	hmax = sqrt(hmax);
	double beta(1.e6 / hmax);
	cout << "hmax: " << hmax << "\n";
	hemax = hmax;

	//assembly
	cout << "# of eles: " << bzmesh.size() << "\n";
#pragma omp parallel for
	for (int e = 0; e<bzmesh.size(); e++)
	{
		if (e != 0 && e % 500 == 0) cout << e << " ";
		//vector<vector<double>> EK;
		//vector<double> Nx, EF;
		//vector<array<double, 3>> dNdx;
		double detJ, Fb;
		array<double, 3> pt;
		if (bzmesh[e].bzflag == 0)//spline element
		{
			vector<vector<double>> EK(bzmesh[e].IEN.size(), vector<double>(bzmesh[e].IEN.size(), 0.));
			vector<double> EF(bzmesh[e].IEN.size(), 0.);
			vector<double> Nx(bzmesh[e].IEN.size());
			vector<array<double, 3>> dNdx(bzmesh[e].IEN.size());
			int i, j, k;
			for (i = 0; i<Gpt.size(); i++)
			{
				for (j = 0; j<Gpt.size(); j++)
				{
					for (k = 0; k < Gpt.size(); k++)
					{
						BasisFunction_IGA(Gpt[i], Gpt[j], Gpt[k], bzmesh[e], Nx, dNdx, detJ);
						detJ = wght[i] * wght[j] * wght[k] * detJ;
						ElementMatrix_1(dNdx, detJ, EK);
						bzmesh[e].Para2Phys(Gpt[i], Gpt[j], Gpt[k], pt.data());
						//Fb = f_source(pt);
						//ElementForce_1(Nx, detJ, Fb, EF);
					}
				}
			}

			//boundary condition, newly added for comparison
			if (bzmesh[e].type == 1)//boundary
			{
				NitscheBoundary_Spline(bzmesh[e], beta, EK, EF);
			}

			//Assembly_Coupling(EK, EF, bzmesh[e].IEN, GK, GF);
			Assembly_Coupling_Bezier(EK, EF, bzmesh[e].IEN, GK, GF);
		}
		else if (bzmesh[e].bzflag == 1 && bzmesh[e].bzcouple == 0)//Bezier element
		{
			vector<vector<double>> EK(bzmesh[e].IENb.size(), vector<double>(bzmesh[e].IENb.size(), 0.));
			vector<double> EF(bzmesh[e].IENb.size(), 0.);
			vector<double> Nx(bzmesh[e].IENb.size());
			vector<array<double, 3>> dNdx(bzmesh[e].IENb.size());
			int i, j, k;
			for (k = 0; k<Gpt.size(); k++)
			{
				for (j = 0; j<Gpt.size(); j++)
				{
					for (i = 0; i < Gpt.size(); i++)
					{
						BasisFunction_Bezier(Gpt[i], Gpt[j], Gpt[k], bzmesh[e].pts, Nx, dNdx, detJ);
						detJ = wght[i] * wght[j] * wght[k] * detJ;
						ElementMatrix_1(dNdx, detJ, EK);
						bzmesh[e].Para2Phys(Gpt[i], Gpt[j], Gpt[k], pt.data());
						Fb = f_source(pt);
						ElementForce_1(Nx, detJ, Fb, EF);
					}
				}
			}
			//boundary condition
			if (bzmesh[e].bcflag == 1)
			{
				NitscheBoundary_Bezier(bzmesh[e], beta, EK, EF);
			}

			//Assembly_Coupling(EK, EF, bzmesh[e].IENb, GK, GF);
			Assembly_Coupling_Bezier(EK, EF, bzmesh[e].IENb, GK, GF);
		}

		if (bzmesh[e].bzcouple == 1)//coupling element, bzflag must be 1
		{
			vector<int> IEN_all(bzmesh[e].IEN.size() + bzmesh[e].IENb.size());
			for (uint i = 0; i < bzmesh[e].IEN.size(); i++) IEN_all[i] = bzmesh[e].IEN[i];
			for (uint i = 0; i < bzmesh[e].IENb.size(); i++) IEN_all[i + bzmesh[e].IEN.size()] = bzmesh[e].IENb[i];
			vector<vector<double>> EK(IEN_all.size(), vector<double>(IEN_all.size(), 0.));
			vector<double> EF(IEN_all.size(), 0.);
			//Nitsche's method coupling interface
			NitscheInterface_Bezier(bzmesh[e], beta, EK, EF);

			//Assembly_Coupling(EK, EF, IEN_all, GK, GF);
			Assembly_Coupling_Bezier(EK, EF, IEN_all, GK, GF);
		}
	}
	cout << "Done assembly!\n";

	//cout << "GK:\n";
	//cout << GK << "\n";
	//cout << "GF:\n";
	//cout << GF << "\n";
	//getchar();

	//VectorXd sol0 = VectorXd::Zero(neq);
	//int loc(0);
	//for (int i = 0; i < npt; i++)
	//{
	//	if (IDBC[i] != -1)
	//	{
	//		sol0[loc] = u0[i];
	//		loc++;
	//	}
	//}
	//VectorXd rmn = GK*sol0 - GF;
	//cout << rmn << "\n"; getchar();

	//ofstream fout;
	////fout.open("../io/ptest3/mat.txt");
	//fout.open("../io/ptest3/gf.txt");
	//if (fout.is_open())
	//{
	//	fout << GF << "\n";
	//	fout.close();
	//}
	//cout << "done output matrix!\n";
	//getchar();
}

void Laplace::Assembly_Coupling_Bezier(const vector<vector<double>>& EK, const vector<double>& EF, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF)
{
	int i, j, A, B;
	for (i = 0; i<IEN.size(); i++)
	{
		for (j = 0; j<IEN.size(); j++)
		{
			A = IEN[i]; B = IEN[j];
			//if (A != -1 && B != -1)
			if (A != -1 && B != -1 && A >= B)//lower triangle part of matrix
			{
				if (IDBC[A] != -1 && IDBC[B] != -1)
				{
					GK.coeffRef(IDBC[A], IDBC[B]) += EK[i][j];
				}
				//else if (IDBC[A] != -1 && IDBC[B] == -1)//strongly imposed
				//{
				//	GF(IDBC[A]) -= EK[i][j] * gh[B];
				//}
			}
			//else if (IDBC[A] != -1 && IDBC[B] == -1)
			//{
			//	GF(IDBC[A]) -= EK[i][j] * gh[B];
			//}
		}
	}

	for (i = 0; i < IEN.size(); i++)//strongly impose BC
	{
		A = IEN[i];
		if (A != -1 && IDBC[A] != -1)
		{
			GF(IDBC[A]) += EF[i];
			for (j = 0; j < IEN.size(); j++)
			{
				B = IEN[j];
				if (B != -1 && IDBC[B] == -1)
				{
					GF(IDBC[A]) -= EK[i][j] * gh[B];
				}
			}
		}
	}
	//for (i = 0; i < IEN.size(); i++)
	//{
	//	A = IEN[i];
	//	if (A != -1 && IDBC[A] != -1)
	//	{
	//		GF(IDBC[A]) += EF[i];
	//	}
	//}
}

void Laplace::NitscheBoundary_Bezier(const BezierElement3D& cell, double beta, vector<vector<double>>& EK, vector<double>& EF)
{
	vector<double> Nx(cell.IENb.size());
	vector<array<double, 3>> dNdx(cell.IENb.size());
	double bzpts[16][3];
	double detJ, Fb;
	array<double, 3> x;
	int i, j;
	double nm[3] = { 0., 0., 0. };
	int ptfloc[6][16] = { { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 }, { 0, 1, 2, 3, 16, 17, 18, 19, 32, 33, 34, 35, 48, 49, 50, 51 },
	{ 3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63 }, { 12, 13, 14, 15, 28, 29, 30, 31, 44, 45, 46, 47, 60, 61, 62, 63 },
	{ 0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60 }, { 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63 } };

	if (cell.bc[0] == 1)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[0][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[0][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[0][ip]][2];
		}
		for (i = 0; i < Gpt.size(); i++)
		{
			for (j = 0; j < Gpt.size(); j++)
			{
				BasisFunction_Bezier(Gpt[j], Gpt[i], 0., cell.pts, Nx, dNdx, detJ);
				GetFaceInfo(Gpt[j], Gpt[i], bzpts, detJ, nm);
				detJ = wght[i] * wght[j] * detJ;
				nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
				ElementMatrix_Bezier_Boundary(Nx, dNdx, nm, beta, detJ, EK);
				cell.Para2Phys(Gpt[j], Gpt[i], 0., x.data());
				Fb = exact_sol(x);
				ElementForce_Bezier_Boundary(Nx, dNdx, nm, beta, detJ, Fb, EF);
			}
		}
	}
	if (cell.bc[1] == 1)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[1][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[1][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[1][ip]][2];
		}
		for (i = 0; i < Gpt.size(); i++)
		{
			for (j = 0; j < Gpt.size(); j++)
			{
				BasisFunction_Bezier(Gpt[j], 0., Gpt[i], cell.pts, Nx, dNdx, detJ);
				GetFaceInfo(Gpt[j], Gpt[i], bzpts, detJ, nm);
				detJ = wght[i] * wght[j] * detJ;
				//nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
				ElementMatrix_Bezier_Boundary(Nx, dNdx, nm, beta, detJ, EK);
				cell.Para2Phys(Gpt[j], 0., Gpt[i], x.data());
				Fb = exact_sol(x);
				ElementForce_Bezier_Boundary(Nx, dNdx, nm, beta, detJ, Fb, EF);
			}
		}
	}
	if (cell.bc[2] == 1)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[2][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[2][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[2][ip]][2];
		}
		for (i = 0; i < Gpt.size(); i++)
		{
			for (j = 0; j < Gpt.size(); j++)
			{
				BasisFunction_Bezier(1., Gpt[j], Gpt[i], cell.pts, Nx, dNdx, detJ);
				GetFaceInfo(Gpt[j], Gpt[i], bzpts, detJ, nm);
				detJ = wght[i] * wght[j] * detJ;
				//nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
				ElementMatrix_Bezier_Boundary(Nx, dNdx, nm, beta, detJ, EK);
				cell.Para2Phys(1., Gpt[j], Gpt[i], x.data());
				Fb = exact_sol(x);
				ElementForce_Bezier_Boundary(Nx, dNdx, nm, beta, detJ, Fb, EF);
			}
		}
	}
	if (cell.bc[3] == 1)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[3][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[3][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[3][ip]][2];
		}
		for (i = 0; i < Gpt.size(); i++)
		{
			for (j = 0; j < Gpt.size(); j++)
			{
				BasisFunction_Bezier(Gpt[j], 1., Gpt[i], cell.pts, Nx, dNdx, detJ);
				GetFaceInfo(Gpt[j], Gpt[i], bzpts, detJ, nm);
				detJ = wght[i] * wght[j] * detJ;
				nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
				ElementMatrix_Bezier_Boundary(Nx, dNdx, nm, beta, detJ, EK);
				cell.Para2Phys(Gpt[j], 1., Gpt[i], x.data());
				Fb = exact_sol(x);
				ElementForce_Bezier_Boundary(Nx, dNdx, nm, beta, detJ, Fb, EF);
			}
		}
	}
	if (cell.bc[4] == 1)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[4][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[4][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[4][ip]][2];
		}
		for (i = 0; i < Gpt.size(); i++)
		{
			for (j = 0; j < Gpt.size(); j++)
			{
				BasisFunction_Bezier(0., Gpt[j], Gpt[i], cell.pts, Nx, dNdx, detJ);
				GetFaceInfo(Gpt[j], Gpt[i], bzpts, detJ, nm);
				detJ = wght[i] * wght[j] * detJ;
				nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
				ElementMatrix_Bezier_Boundary(Nx, dNdx, nm, beta, detJ, EK);
				cell.Para2Phys(0., Gpt[j], Gpt[i], x.data());
				Fb = exact_sol(x);
				ElementForce_Bezier_Boundary(Nx, dNdx, nm, beta, detJ, Fb, EF);
			}
		}
	}
	if (cell.bc[5] == 1)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[5][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[5][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[5][ip]][2];
		}
		for (i = 0; i < Gpt.size(); i++)
		{
			for (j = 0; j < Gpt.size(); j++)
			{
				BasisFunction_Bezier(Gpt[j], Gpt[i], 1., cell.pts, Nx, dNdx, detJ);
				GetFaceInfo(Gpt[j], Gpt[i], bzpts, detJ, nm);
				detJ = wght[i] * wght[j] * detJ;
				//nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
				ElementMatrix_Bezier_Boundary(Nx, dNdx, nm, beta, detJ, EK);
				cell.Para2Phys(Gpt[j], Gpt[i], 1., x.data());
				Fb = exact_sol(x);
				ElementForce_Bezier_Boundary(Nx, dNdx, nm, beta, detJ, Fb, EF);
			}
		}
	}
}

void Laplace::ElementMatrix_Bezier_Boundary(const vector<double>& Nx, const vector<array<double, 3>>& dNdx, double nm[3], double beta, double detJ, vector<vector<double>>& EK)
{
	//only meshfree nodes have contribution on boundary since domain is split
	unsigned int i, j;
	double tmp1(0.), tmp2(0.), tmp3(0.);
	for (i = 0; i<dNdx.size(); i++)
	{
		for (j = 0; j<dNdx.size(); j++)
		{
			tmp1 = (dNdx[i][0] * nm[0] + dNdx[i][1] * nm[1] + dNdx[i][2] * nm[2])*Nx[j];
			tmp2 = (dNdx[j][0] * nm[0] + dNdx[j][1] * nm[1] + dNdx[j][2] * nm[2])*Nx[i];
			tmp3 = beta * Nx[i] * Nx[j];
			EK[i][j] += (-tmp1 - tmp2 + tmp3)*detJ;
		}
	}
}

void Laplace::ElementForce_Bezier_Boundary(const vector<double>& Nx, const vector<array<double, 3>>& dNdx, double nm[3], double beta, double detJ, double Fb, vector<double>& EF)
{
	double tmp1(0.), tmp2(0.);
	//double sum(0.);
	for (unsigned int i = 0; i<Nx.size(); i++)
	{
		tmp1 = Fb*(dNdx[i][0] * nm[0] + dNdx[i][1] * nm[1] + dNdx[i][2] * nm[2]);
		tmp2 = beta*Fb*Nx[i];
		EF[i] += (-tmp1 + tmp2) * detJ;
		//sum += EF[i];
	}
}

void Laplace::NitscheInterface_Bezier(const BezierElement3D& cell, double beta, vector<vector<double>>& EK, vector<double>& EF)
{
	vector<double> Nx(cell.IEN.size() + cell.IENb.size());
	vector<array<double, 3>> dNdx(cell.IEN.size() + cell.IENb.size());
	double bzpts[16][3];
	double detJ, Fb;
	array<double, 3> x;
	int i, j;
	double nm[3] = { 0., 0., 0. };
	int ptfloc[6][16] = { { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 }, { 0, 1, 2, 3, 16, 17, 18, 19, 32, 33, 34, 35, 48, 49, 50, 51 },
	{ 3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63 }, { 12, 13, 14, 15, 28, 29, 30, 31, 44, 45, 46, 47, 60, 61, 62, 63 },
	{ 0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60 }, { 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63 } };
	if (cell.bc[0] == 2)
	{
		//cout << "bc face: 0\n";
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[0][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[0][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[0][ip]][2];
		}
		for (i = 0; i < Gpt.size(); i++)
		{
			for (j = 0; j < Gpt.size(); j++)
			{
				BasisFunction_BZAll(Gpt[j], Gpt[i], 0., cell, Nx, dNdx, detJ);
				GetFaceInfo(Gpt[j], Gpt[i], bzpts, detJ, nm);
				detJ = wght[i] * wght[j] * detJ;
				nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
				ElementMatrix_Bezier_Interface(cell.IEN.size(), Nx, dNdx, nm, beta, detJ, EK);
			}
		}
	}
	if (cell.bc[1] == 2)
	{
		//cout << "bc face: 1\n";
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[1][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[1][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[1][ip]][2];
		}
		for (i = 0; i < Gpt.size(); i++)
		{
			for (j = 0; j < Gpt.size(); j++)
			{
				BasisFunction_BZAll(Gpt[j], 0., Gpt[i], cell, Nx, dNdx, detJ);
				GetFaceInfo(Gpt[j], Gpt[i], bzpts, detJ, nm);
				detJ = wght[i] * wght[j] * detJ;
				ElementMatrix_Bezier_Interface(cell.IEN.size(), Nx, dNdx, nm, beta, detJ, EK);
			}
		}
	}
	if (cell.bc[2] == 2)
	{
		//cout << "bc face: 2\n";
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[2][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[2][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[2][ip]][2];
		}
		for (i = 0; i < Gpt.size(); i++)
		{
			for (j = 0; j < Gpt.size(); j++)
			{
				BasisFunction_BZAll(1., Gpt[j], Gpt[i], cell, Nx, dNdx, detJ);
				GetFaceInfo(Gpt[j], Gpt[i], bzpts, detJ, nm);
				detJ = wght[i] * wght[j] * detJ;
				ElementMatrix_Bezier_Interface(cell.IEN.size(), Nx, dNdx, nm, beta, detJ, EK);
			}
		}
	}
	if (cell.bc[3] == 2)
	{
		//cout << "bc face: 3\n";
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[3][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[3][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[3][ip]][2];
		}
		for (i = 0; i < Gpt.size(); i++)
		{
			for (j = 0; j < Gpt.size(); j++)
			{
				BasisFunction_BZAll(Gpt[j], 1., Gpt[i], cell, Nx, dNdx, detJ);
				GetFaceInfo(Gpt[j], Gpt[i], bzpts, detJ, nm);
				detJ = wght[i] * wght[j] * detJ;
				nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
				ElementMatrix_Bezier_Interface(cell.IEN.size(), Nx, dNdx, nm, beta, detJ, EK);
			}
		}
	}
	if (cell.bc[4] == 2)
	{
		//cout << "bc face: 4\n";
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[4][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[4][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[4][ip]][2];
		}
		for (i = 0; i < Gpt.size(); i++)
		{
			for (j = 0; j < Gpt.size(); j++)
			{
				BasisFunction_BZAll(0., Gpt[j], Gpt[i], cell, Nx, dNdx, detJ);
				GetFaceInfo(Gpt[j], Gpt[i], bzpts, detJ, nm);
				detJ = wght[i] * wght[j] * detJ;
				nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
				ElementMatrix_Bezier_Interface(cell.IEN.size(), Nx, dNdx, nm, beta, detJ, EK);
			}
		}
	}
	if (cell.bc[5] == 2)
	{
		//cout << "bc face: 5\n";
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[5][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[5][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[5][ip]][2];
		}
		for (i = 0; i < Gpt.size(); i++)
		{
			for (j = 0; j < Gpt.size(); j++)
			{
				BasisFunction_BZAll(Gpt[j], Gpt[i], 1., cell, Nx, dNdx, detJ);
				GetFaceInfo(Gpt[j], Gpt[i], bzpts, detJ, nm);
				detJ = wght[i] * wght[j] * detJ;
				ElementMatrix_Bezier_Interface(cell.IEN.size(), Nx, dNdx, nm, beta, detJ, EK);
			}
		}
	}
}

void Laplace::ElementMatrix_Bezier_Interface(int mfid, const vector<double>& Nx, const vector<array<double, 3>>& dNdx, double nm[3], double beta, double detJ, vector<vector<double>>& EK)
{
	unsigned int i, j;
	double tmp1(0.), tmp2(0.), tmp3(0.);
	//cout << mfid << " " << Nx.size() << " " << dNdx.size() << "\n";
	for (i = mfid; i<dNdx.size(); i++)//symmetry term
	{
		for (j = mfid; j<dNdx.size(); j++)
		{
			tmp1 = (dNdx[i][0] * nm[0] + dNdx[i][1] * nm[1] + dNdx[i][2] * nm[2])*Nx[j];
			EK[i][j] += (-tmp1)*detJ;
		}
		for (j = 0; j < mfid; j++)
		{
			tmp1 = (dNdx[i][0] * nm[0] + dNdx[i][1] * nm[1] + dNdx[i][2] * nm[2])*Nx[j];
			EK[i][j] += tmp1*detJ;
		}
	}
	for (i = mfid; i<dNdx.size(); i++)//consistency term
	{
		for (j = mfid; j<dNdx.size(); j++)
		{
			tmp2 = (dNdx[j][0] * nm[0] + dNdx[j][1] * nm[1] + dNdx[j][2] * nm[2])*Nx[i];
			EK[i][j] += (-tmp2)*detJ;
		}
	}
	for (i = 0; i<mfid; i++)//consistency term
	{
		for (j = mfid; j<dNdx.size(); j++)
		{
			tmp2 = (dNdx[j][0] * nm[0] + dNdx[j][1] * nm[1] + dNdx[j][2] * nm[2])*Nx[i];
			EK[i][j] += tmp2*detJ;
		}
	}
	for (i = mfid; i<dNdx.size(); i++)//penalty term
	{
		for (j = mfid; j<dNdx.size(); j++)
		{
			tmp3 = beta * Nx[i] * Nx[j];
			EK[i][j] += tmp3*detJ;
		}
	}
	for (i = mfid; i<dNdx.size(); i++)//penalty term
	{
		for (j = 0; j<mfid; j++)
		{
			tmp3 = beta * Nx[i] * Nx[j];
			EK[i][j] -= tmp3*detJ;
		}
	}
	for (i = 0; i<mfid; i++)//penalty term
	{
		for (j = mfid; j<dNdx.size(); j++)
		{
			tmp3 = beta * Nx[i] * Nx[j];
			EK[i][j] -= tmp3*detJ;
		}
	}
	for (i = 0; i<mfid; i++)//penalty term
	{
		for (j = 0; j<mfid; j++)
		{
			tmp3 = beta * Nx[i] * Nx[j];
			EK[i][j] += tmp3*detJ;
		}
	}
}

void Laplace::VisualizeDisp_Coupling_Bezier(const vector<BezierElement3D>& bzmesh, string fn)
{
	vector<array<double, 3>> spt;//sample points
	vector<double> sdisp;
	vector<array<int, 8>> sele;
	vector<double> errL2;
	vector<array<double, 3>> lpt;//visulize parameter lines
	vector<array<int, 2>> led;//line connectivity
	double detJ;

	for (unsigned int e = 0; e<bzmesh.size(); e++)
	{
		int ns(4);
		if (bzmesh[e].type == 1) ns = 5;
		vector<double> su(ns);
		for (int i = 0; i < ns; i++)
		{
			su[i] = double(i) / (double(ns) - 1.);
		}

		int loc(0);
		int pstart = spt.size();
		for (int a = 0; a<ns; a++)
		{
			for (int b = 0; b<ns; b++)
			{
				for (int c = 0; c < ns; c++)
				{
					double pt1[3], dudx[3];
					double disp;
					//bzmesh[e].Para2Phys(su[c], su[b], su[a], pt1);
					DispCal_Coupling_Bezier(su[c], su[b], su[a], bzmesh[e], pt1, disp, dudx, detJ);
					array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
					spt.push_back(pt);
					sdisp.push_back(disp);
				}
			}
		}
		int nns[2] = { ns*ns*ns, ns*ns };
		for (int a = 0; a<ns - 1; a++)
		{
			for (int b = 0; b<ns - 1; b++)
			{
				for (int c = 0; c < ns - 1; c++)
				{
					array<int, 8> el;
					el[0] = pstart + a*nns[1] + b*ns + c;
					el[1] = pstart + a*nns[1] + b*ns + c + 1;
					el[2] = pstart + a*nns[1] + (b + 1)*ns + c + 1;
					el[3] = pstart + a*nns[1] + (b + 1)*ns + c;
					el[4] = pstart + (a + 1)*nns[1] + b*ns + c;
					el[5] = pstart + (a + 1)*nns[1] + b*ns + c + 1;
					el[6] = pstart + (a + 1)*nns[1] + (b + 1)*ns + c + 1;
					el[7] = pstart + (a + 1)*nns[1] + (b + 1)*ns + c;
					sele.push_back(el);
				}
			}
		}
		//edges
		int lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(su[a], 0., 0., pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(su[a], 1., 0., pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(su[a], 0., 1., pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(su[a], 1., 1., pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(0., su[a], 0., pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(1., su[a], 0., pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(0., su[a], 1., pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(1., su[a], 1., pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(0., 0., su[a], pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(1., 0., su[a], pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(0., 1., su[a], pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
		lstart = lpt.size();
		for (int a = 0; a < ns; a++)
		{
			double pt1[3];
			bzmesh[e].Para2Phys(1., 1., su[a], pt1);
			array<double, 3> pt = { pt1[0], pt1[1], pt1[2] };
			lpt.push_back(pt);
		}
		for (int a = 0; a < ns - 1; a++)
		{
			array<int, 2> ed = { lstart + a, lstart + a + 1 };
			led.push_back(ed);
		}
	}

	string fname = fn + "_disp.vtk";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << spt.size() << " float\n";
		for (i = 0; i<spt.size(); i++)
		{
			fout << spt[i][0] << " " << spt[i][1] << " " << spt[i][2] << "\n";
		}
		fout << "\nCELLS " << sele.size() << " " << 9 * sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "8 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3]
				<< " " << sele[i][4] << " " << sele[i][5] << " " << sele[i][6] << " " << sele[i][7] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "12\n";
		}
		//fout<<"\nPOINT_DATA "<<sdisp.size()<<"\nVECTORS disp float\n";
		//for(i=0;i<sdisp.size();i++)
		//{
		//	fout << sdisp[i][0] << " " << sdisp[i][1] << " " << sdisp[i][2] << "\n";
		//}
		//fout << "\nPOINT_DATA " << sse.size() << "\nVECTORS strain float\n";
		//for (i = 0; i<sse.size(); i++)
		//{
		//	fout << sse[i][0] << " " << sse[i][1] << " " << sse[i][2] << "\n";
		//}
		//fout<<"POINT_DATA "<<sss.size()<<"\nVECTORS stress float\n";
		//for(i=0;i<sss.size();i++)
		//{
		//	fout<<sss[i][0]<<" "<<sss[i][1]<<" "<<sss[i][2]<<"\n";
		//}
		//fout<<"POINT_DATA "<<sdisp_err.size()<<"\nVECTORS disp_err float\n";
		//for(i=0;i<sdisp_err.size();i++)
		//{
		//	fout<<sdisp_err[i][0]<<" "<<sdisp_err[i][1]<<" 0\n";
		//}

		fout << "POINT_DATA " << sdisp.size() << "\nSCALARS err float 1\nLOOKUP_TABLE default\n";
		for (uint i = 0; i<sdisp.size(); i++)
		{
			fout << sdisp[i] << "\n";
		}

		//fout<<"\nCELL_DATA "<<errL2.size()<<"\nSCALARS Error float 1\nLOOKUP_TABLE default\n";
		//for(i=0;i<errL2.size();i++)
		//{
		//	fout<<errL2[i]<<"\n";
		//	//fout<<eles[i].type<<"\n";
		//}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}

	string fname1(fn + "-lines.vtk");
	ofstream fout1;
	fout1.open(fname1.c_str());
	if (fout1.is_open())
	{
		fout1 << "# vtk DataFile Version 2.0\nSquare plate test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout1 << "POINTS " << lpt.size() << " float\n";
		for (i = 0; i<lpt.size(); i++)
		{
			fout1 << lpt[i][0] << " " << lpt[i][1] << " " << lpt[i][2] << "\n";
		}
		fout1 << "\nCELLS " << led.size() << " " << 3 * led.size() << '\n';
		for (i = 0; i<led.size(); i++)
		{
			fout1 << "2 " << led[i][0] << " " << led[i][1] << '\n';
		}
		fout1 << "\nCELL_TYPES " << led.size() << '\n';
		for (i = 0; i<led.size(); i++)
		{
			fout1 << "3\n";
		}
		fout1.close();
	}
	else
	{
		cout << "Cannot open " << fname1 << "!\n";
	}
}

void Laplace::DispCal_Coupling_Bezier(double u, double v, double w, const BezierElement3D& bzel, double pt[3], double& disp, double dudx[3], double& detJ)
{
	if (bzel.bzflag == 0)
	{
		vector<double> Nx(bzel.IEN.size());
		vector<array<double, 3>> dNdx(bzel.IEN.size());
		bzel.Para2Phys(u, v, w, pt);
		BasisFunction_IGA(u, v, w, bzel, Nx, dNdx, detJ);
		disp = 0.;
		dudx[0] = 0.; dudx[1] = 0.; dudx[2] = 0.;
		for (uint i = 0; i < bzel.IEN.size(); i++)
		{
			disp += Nx[i] * uh[bzel.IEN[i]];
			dudx[0] += dNdx[i][0] * uh[bzel.IEN[i]];
			dudx[1] += dNdx[i][1] * uh[bzel.IEN[i]];
			dudx[2] += dNdx[i][2] * uh[bzel.IEN[i]];
		}
	}
	else
	{
		vector<double> Nx(bzel.IENb.size());
		vector<array<double, 3>> dNdx(bzel.IENb.size());
		bzel.Para2Phys(u, v, w, pt);
		BasisFunction_Bezier(u, v, w, bzel.pts, Nx, dNdx, detJ);
		disp = 0.;
		dudx[0] = 0.; dudx[1] = 0.; dudx[2] = 0.;
		for (uint i = 0; i < bzel.IENb.size(); i++)
		{
			disp += Nx[i] * uh[bzel.IENb[i]];
			dudx[0] += dNdx[i][0] * uh[bzel.IENb[i]];
			dudx[1] += dNdx[i][1] * uh[bzel.IENb[i]];
			dudx[2] += dNdx[i][2] * uh[bzel.IENb[i]];
		}
	}
}

void Laplace::ErrorEstimate_Coupling_Bezier(const vector<BezierElement3D>& bzmesh, vector<vector<double>>& err)
{
	double hmax(0.), dst;
	int bzop[4][2] = { { 0, 63 }, { 3, 60 }, { 15, 48 }, { 12, 51 } };
	for (int e = 0; e < bzmesh.size(); e++)
	{
		for (int i = 0; i<4; i++)
		{
			dst = (bzmesh[e].pts[bzop[i][0]][0] - bzmesh[e].pts[bzop[i][1]][0])*(bzmesh[e].pts[bzop[i][0]][0] - bzmesh[e].pts[bzop[i][1]][0]) +
				(bzmesh[e].pts[bzop[i][0]][1] - bzmesh[e].pts[bzop[i][1]][1])*(bzmesh[e].pts[bzop[i][0]][1] - bzmesh[e].pts[bzop[i][1]][1]) +
				(bzmesh[e].pts[bzop[i][0]][2] - bzmesh[e].pts[bzop[i][1]][2])*(bzmesh[e].pts[bzop[i][0]][2] - bzmesh[e].pts[bzop[i][1]][2]);
			if (dst>hmax) hmax = dst;
		}

	}
	hmax = sqrt(hmax);
	cout << "hmax: " << hmax << "\n";
	hemax = hmax;

	cout << "Error estimating...\n";
	//err.clear();
	err[0].resize(bzmesh.size(), 0.);
	err[1].resize(bzmesh.size(), 0.);
#pragma omp parallel for
	for (int e = 0; e < bzmesh.size(); e++)
	{
		double L2, H1;
		ElementError_Coupling_Bezier(bzmesh[e], L2, H1);
		err[0][e] = L2;//squared
		err[1][e] = H1;//squared
		//cout << err[e] << "\n"; getchar();
	}
}

void Laplace::ElementError_Coupling_Bezier(const BezierElement3D& bzel, double& L2, double& H1)
{
	L2 = 0.; H1 = 0.;
	uint i, j, k;
	double disp, detJ, sol_e;
	array<double, 3> pt, dudx, dudx_e;
	for (i = 0; i < Gpt.size(); i++)
	{
		for (j = 0; j < Gpt.size(); j++)
		{
			for (k = 0; k < Gpt.size(); k++)
			{
				DispCal_Coupling_Bezier(Gpt[i], Gpt[j], Gpt[k], bzel, pt.data(), disp, dudx.data(), detJ);
				sol_e = exact_sol(pt);
				grad_sol(pt, dudx_e);
				detJ = wght[i] * wght[j] * wght[k] * detJ;
				L2 += (detJ*(disp - sol_e)*(disp - sol_e));
				H1 += detJ*((dudx[0] - dudx_e[0])*(dudx[0] - dudx_e[0]) + (dudx[1] - dudx_e[1])*(dudx[1] - dudx_e[1]) + (dudx[2] - dudx_e[2])*(dudx[2] - dudx_e[2]));
			}
		}
	}
}

void Laplace::OutputError_All(const vector<vector<double>>& err, string fn)
{
	double L2all(0.), H1all(0.);
	for (uint i = 0; i < err[0].size(); i++)
	{
		L2all += err[0][i];
		H1all += err[1][i];
	}
	L2all = sqrt(L2all);
	H1all = sqrt(H1all);

	cout << "DOF: " << npt << "\n";
	cout << "L2: " << hemax << " " << L2all << "\n";
	cout << "H1: " << hemax << " " << H1all << "\n";

	ofstream fout;
	string fname(fn + "_err.txt");
	fout.open(fname.c_str());
	if (fout.is_open())
	{
		fout << npt << " " << L2all << "\n";
		fout << npt << " " << H1all << "\n";
		fout << hemax << " " << L2all << "\n";
		fout << hemax << " " << H1all << "\n";
		fout << neq << "\n";
		fout.close();
	}
	else
	{
		cerr << "Can't open " << fname << "\n";
	}
}

void Laplace::Solver_Coupling_Bezier(SparseMatrix<double>& GK, VectorXd& GF)
{
	cout << "Solving linear system...\n";

	clock_t begin = clock();

	SimplicialLDLT<SparseMatrix<double>> solver;
	VectorXd sol = solver.compute(GK).solve(GF);

	//SimplicialLLT<SparseMatrix<double>> solver;
	//VectorXd sol = solver.compute(GK).solve(GF);

	////ConjugateGradient<SparseMatrix<double>, Lower | Upper, IdentityPreconditioner> solver;
	//ConjugateGradient<SparseMatrix<double>, Lower | Upper> solver;
	////ConjugateGradient<SparseMatrix<double>, Lower | Upper, SimplicialCholesky<SparseMatrix<double>>> solver;//very slow
	////double tol(double(neq)*10.e-16);
	//double tol(1.e-6);
	////cout << "tol before: " << solver.tolerance() << "\n";
	//solver.setTolerance(tol);
	////cout << "tol after: " << solver.tolerance() << "\n";
	//solver.setMaxIterations(10 * npt);
	//solver.compute(GK);
	//cout << "done compute!\n";
	////VectorXd sol = solver.solve(GF);
	////VectorXd x0 = VectorXd::Ones(neq);
	//VectorXd sol = solver.solveWithGuess(GF, u0);
	//cout << "\n# iterations: " << solver.iterations() << '\n';
	//cout << "estimated error: " << solver.error() << '\n';
	////getchar();

	//BiCGSTAB<SparseMatrix<double>, IncompleteLUT<double>> solver;
	//solver.compute(GK);
	//cout << "done compute!\n";
	//VectorXd sol = solver.solve(GF);
	//cout << "\n# iterations: " << solver.iterations() << '\n';
	//cout << "estimated error: " << solver.error() << '\n';

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << "solver time: " << elapsed_secs << "\n";

	uh.resize(npt);
	for (int i = 0; i<npt; i++)
	{
		if (IDBC[i] != -1)
		{
			uh[i] = sol(IDBC[i]);
		}
		else
		{
			//uh[i] = 0.;
			uh[i] = gh[i];
		}
		//cout << IDBC[i] << " ";
		//cout << uh[i] << " ";
	}
	//getchar();
	//cout << sol << "\n"; getchar();

	cout << "Done solving!\n";
}

void Laplace::Solver_MKL(SparseMatrix<double>& GK, VectorXd& GF)
{
	cout << "\n\nSolving linear system...\n";

	MKL_INT n(GF.size());//neq
	vector<MKL_INT> ia(GK.outerSize() + 1), ja(GK.nonZeros());//ia.size()==GK.outerSize()+1, ja.size()==GK.nonZeros()
	vector<double>	a(GK.nonZeros()), b(GF.size()), x(GF.size());//a.size()==GK.nonZeros(), b.size()==GF.size(), sol.size()==GF.size()
																 //assign matrix and free memory of GK and GF
	llint count(0);
	for (llint k = 0; k < GK.outerSize(); k++)//k is the column id
	{
		ia[k] = count + 1;
		for (SparseMatrix<double>::InnerIterator it(GK, k); it; ++it)
		{
			if (it.row() >= k)//symmetric
			{
				ja[count] = it.row() + 1;
				a[count] = it.value();
				count++;
			}
		}
	}
	ia[GK.outerSize()] = GK.nonZeros() + 1;
	for (uint k = 0; k < GF.size(); k++)
	{
		b[k] = GF[k];
	}

	clock_t begin = clock();
	//MKL PARDISO
	MKL_INT mtype = -2; /* Real symmetric matrix */
						//MKL_INT mtype = 2; /* Real symmetric pos def matrix */
	MKL_INT nrhs = 1; /* Number of right hand sides. */
					  /* Internal solver memory pointer pt, */
					  /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
					  /* or void *pt[64] should be OK on both architectures */
	void *pt[64];
	/* Pardiso control parameters. */
	MKL_INT iparm[64];
	MKL_INT maxfct, mnum, phase, error, msglvl;
	/* Auxiliary variables. */
	MKL_INT i;
	double ddum; /* Double dummy */
	MKL_INT idum; /* Integer dummy. */
				  /* -------------------------------------------------------------------- */
				  /* .. Setup Pardiso control parameters. */
				  /* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++) {
		iparm[i] = 0;
	}
	iparm[0] = 1; /* No solver default */
	iparm[1] = 2; /* Fill-in reordering from METIS */
				  /* Numbers of processors, value of OMP_NUM_THREADS */
	iparm[2] = 1;
	iparm[3] = 0; /* No iterative-direct algorithm */
	iparm[4] = 0; /* No user fill-in reducing permutation */
	iparm[5] = 0; /* Write solution into x */
	iparm[6] = 0; /* Not in use */
	iparm[7] = 2; /* Max numbers of iterative refinement steps */
	iparm[8] = 0; /* Not in use */
	iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
	iparm[11] = 0; /* Not in use */
	iparm[12] = 0; /* Not in use */
	iparm[13] = 0; /* Output: Number of perturbed pivots */
	iparm[14] = 0; /* Not in use */
	iparm[15] = 0; /* Not in use */
	iparm[16] = 0; /* Not in use */
	iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1; /* Output: Mflops for LU factorization */
	iparm[19] = 0; /* Output: Numbers of CG Iterations */
	iparm[59] = 0; /* Out-Of-Core mode, default is In-Core model 0 */
	maxfct = 1; /* Maximum number of numerical factorizations. */
	mnum = 1; /* Which factorization to use. */
	msglvl = 1; /* Print statistical information in file */
	error = 0; /* Initialize error flag */
			   /* -------------------------------------------------------------------- */
			   /* .. Initialize the internal solver memory pointer. This is only */
			   /* necessary for the FIRST call of the PARDISO solver. */
			   /* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++) {
		pt[i] = 0;
	}
	/* -------------------------------------------------------------------- */
	/* .. Reordering and Symbolic Factorization. This step also allocates */
	/* all memory that is necessary for the factorization. */
	/* -------------------------------------------------------------------- */
	phase = 11;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, a.data(), ia.data(), ja.data(), &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0) {
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	}
	printf("\nReordering completed ... ");
	printf("\nNumber of nonzeros in factors = %d", iparm[17]);
	printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
	/* -------------------------------------------------------------------- */
	/* .. Numerical factorization. */
	/* -------------------------------------------------------------------- */
	phase = 22;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, a.data(), ia.data(), ja.data(), &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);
	if (error != 0) {
		printf("\nERROR during numerical factorization: %d", error);
		exit(2);
	}
	printf("\nFactorization completed ... ");
	/* -------------------------------------------------------------------- */
	/* .. Back substitution and iterative refinement. */
	/* -------------------------------------------------------------------- */
	phase = 33;
	iparm[7] = 2; /* Max numbers of iterative refinement steps. */
				  /* Set right hand side to one. */
				  //for (i = 0; i < n; i++) {
				  //	b[i] = 1;
				  //}
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, a.data(), ia.data(), ja.data(), &idum, &nrhs,
		iparm, &msglvl, b.data(), x.data(), &error);
	if (error != 0) {
		printf("\nERROR during solution: %d", error);
		exit(3);
	}
	printf("\nSolve completed ... ");
	//printf("\nThe solution of the system is: ");
	//for (i = 0; i < n; i++) {
	//	printf("\n x [%d] = % f", i, x[i]);
	//}
	//printf("\n");
	/* -------------------------------------------------------------------- */
	/* .. Termination and release of memory. */
	/* -------------------------------------------------------------------- */
	phase = -1; /* Release internal memory. */
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
		&n, &ddum, ia.data(), ja.data(), &idum, &nrhs,
		iparm, &msglvl, &ddum, &ddum, &error);

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout << "solver time: " << elapsed_secs << "\n";

	uh.resize(npt);
	for (uint k = 0; k<npt; k++)
	{
		if (IDBC[k] != -1)
		{
			uh[k] = x[IDBC[k]];
		}
		else
		{
			uh[k] = gh[k];
		}
	}
	cout << "Done solving!\n";
}

//void Laplace::Solver_Matlab(SparseMatrix<double>& GK, VectorXd& GF)
//{
//	vector<double> irow(GK.nonZeros()), jcol(GK.nonZeros()), val(GK.nonZeros()), rhs(GF.size());
//	int count(0);
//	for (int k = 0; k < GK.outerSize(); k++)
//	{
//		for (SparseMatrix<double>::InnerIterator it(GK, k); it; ++it)
//		{
//			irow[count] = double(it.row() + 1);
//			jcol[count] = double(it.col() + 1);
//			val[count] = it.value();
//			count++;
//		}
//	}
//	for (int i = 0; i < GF.size(); i++)
//	{
//		rhs[i] = GF[i];
//	}
//
//	cout << "Solving linear system...\n";
//	clock_t begin = clock();
//
//	MatlabSolver solver;
//	vector<double> sol(neq);
//	solver.Initilize_SparseMatrix(irow, jcol, val, rhs);
//	//solver.CG_Solver(sol.data());
//	solver.LU_Solver(sol.data());
//
//	clock_t end = clock();
//	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//	cout << "solver time: " << elapsed_secs << "\n";
//
//	uh.resize(npt);
//	for (int i = 0; i<npt; i++)
//	{
//		if (IDBC[i] != -1)
//			uh[i] = sol[IDBC[i]];
//		else
//			uh[i] = 0.;
//	}
//
//	cout << "Done solving!\n";
//}

//void Laplace::Solver_MKL(SparseMatrix<double, Eigen::RowMajor>& GK, VectorXd& GF)
//{
//	MKL_INT n = neq;
//	vector<MKL_INT> ia(GK.outerSize() + 1), ja(GK.nonZeros());
//	vector<double>	a(GK.nonZeros()), b(GF.size()), sol(GF.size());
//	int count(0);
//	for (int k = 0; k < GK.outerSize(); k++)
//	{
//		ia[k] = count + 1;
//		for (SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(GK, k); it; ++it)
//		{
//			ja[count] = double(it.col() + 1);
//			a[count] = it.value();
//			count++;
//		}
//	}
//	ia[GK.outerSize()] = GK.nonZeros();
//	for (int i = 0; i < GF.size(); i++)
//	{
//		b[i] = GF[i];
//	}
//
//	cout << "Solving linear system...\n";
//
//	clock_t begin = clock();
//
//	MKL_INT mtype = 2; /* Real symmetric positive matrix */
//	MKL_INT nrhs = 1; /* Number of right hand sides. */
//	/* Internal solver memory pointer pt, */
//	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
//	/* or void *pt[64] should be OK on both architectures */
//	void *pt[64];
//	/* Pardiso control parameters. */
//	MKL_INT iparm[64];
//	MKL_INT maxfct, mnum, phase, error, msglvl;
//	/* Auxiliary variables. */
//	MKL_INT i;
//	double ddum; /* Double dummy */
//	MKL_INT idum; /* Integer dummy. */
//	/* -------------------------------------------------------------------- */
//	/* .. Setup Pardiso control parameters. */
//	/* -------------------------------------------------------------------- */
//	for (i = 0; i < 64; i++) {
//		iparm[i] = 0;
//	}
//	iparm[0] = 1; /* No solver default */
//	iparm[1] = 2; /* Fill-in reordering from METIS */
//	/* Numbers of processors, value of OMP_NUM_THREADS */
//	iparm[2] = 1;
//	iparm[3] = 0; /* No iterative-direct algorithm */
//	iparm[4] = 0; /* No user fill-in reducing permutation */
//	iparm[5] = 0; /* Write solution into x */
//	iparm[6] = 0; /* Not in use */
//	iparm[7] = 2; /* Max numbers of iterative refinement steps */
//	iparm[8] = 0; /* Not in use */
//	iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
//	iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
//	iparm[11] = 0; /* Not in use */
//	iparm[12] = 0; /* Not in use */
//	iparm[13] = 0; /* Output: Number of perturbed pivots */
//	iparm[14] = 0; /* Not in use */
//	iparm[15] = 0; /* Not in use */
//	iparm[16] = 0; /* Not in use */
//	iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
//	iparm[18] = -1; /* Output: Mflops for LU factorization */
//	iparm[19] = 0; /* Output: Numbers of CG Iterations */
//	maxfct = 1; /* Maximum number of numerical factorizations. */
//	mnum = 1; /* Which factorization to use. */
//	msglvl = 1; /* Print statistical information in file */
//	error = 0; /* Initialize error flag */
//	/* -------------------------------------------------------------------- */
//	/* .. Initialize the internal solver memory pointer. This is only */
//	/* necessary for the FIRST call of the PARDISO solver. */
//	/* -------------------------------------------------------------------- */
//	for (i = 0; i < 64; i++) {
//		pt[i] = 0;
//	}
//	/* -------------------------------------------------------------------- */
//	/* .. Reordering and Symbolic Factorization. This step also allocates */
//	/* all memory that is necessary for the factorization. */
//	/* -------------------------------------------------------------------- */
//	phase = 11;
//	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
//		&n, a.data(), ia.data(), ja.data(), &idum, &nrhs,
//		iparm, &msglvl, &ddum, &ddum, &error);
//	if (error != 0) {
//		printf("\nERROR during symbolic factorization: %d", error);
//		exit(1);
//	}
//	printf("\nReordering completed ... ");
//	printf("\nNumber of nonzeros in factors = %d", iparm[17]);
//	printf("\nNumber of factorization MFLOPS = %d", iparm[18]);
//	/* -------------------------------------------------------------------- */
//	/* .. Numerical factorization. */
//	/* -------------------------------------------------------------------- */
//	phase = 22;
//	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
//		&n, a.data(), ia.data(), ja.data(), &idum, &nrhs,
//		iparm, &msglvl, &ddum, &ddum, &error);
//	if (error != 0) {
//		printf("\nERROR during numerical factorization: %d", error);
//		exit(2);
//	}
//	printf("\nFactorization completed ... ");
//	/* -------------------------------------------------------------------- */
//	/* .. Back substitution and iterative refinement. */
//	/* -------------------------------------------------------------------- */
//	phase = 33;
//	iparm[7] = 2; /* Max numbers of iterative refinement steps. */
//	/* Set right hand side to one. */
//	for (i = 0; i < n; i++) {
//		b[i] = 1;
//	}
//	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
//		&n, a.data(), ia.data(), ja.data(), &idum, &nrhs,
//		iparm, &msglvl, b.data(), sol.data(), &error);
//	if (error != 0) {
//		printf("\nERROR during solution: %d", error);
//		exit(3);
//	}
//	printf("\nSolve completed ... ");
//	//printf("\nThe solution of the system is: ");
//	//for (i = 0; i < n; i++) {
//	//	printf("\n x [%d] = % f", i, x[i]);
//	//}
//	//printf("\n");
//	/* -------------------------------------------------------------------- */
//	/* .. Termination and release of memory. */
//	/* -------------------------------------------------------------------- */
//	phase = -1; /* Release internal memory. */
//	PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
//		&n, &ddum, ia.data(), ja.data(), &idum, &nrhs,
//		iparm, &msglvl, &ddum, &ddum, &error);
//
//
//
//
//
//	clock_t end = clock();
//	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//	cout << "solver time: " << elapsed_secs << "\n";
//
//	uh.resize(npt);
//	for (int i = 0; i<npt; i++)
//	{
//		if (IDBC[i] != -1)
//			uh[i] = sol[IDBC[i]];
//		else
//			uh[i] = 0.;
//	}
//
//	cout << "Done solving!\n";
//}

void Laplace::Write2MKL_SparseMat(SparseMatrix<double>& GK, VectorXd& GF, string fn)
{
	//MKL sparse matrix is row major but eigen is colomn major. However, doesn't matter here because GK is symmetric
	//vector<double> irow(GK.nonZeros()), jcol(GK.outerSize() + 1), val(GK.nonZeros());
	//llint nnz((GK.nonZeros() + GK.outerSize()) / 2);
	llint nnz(GK.nonZeros());
	vector<llint> irow(nnz), jcol(GK.outerSize() + 1);
	vector<double> val(nnz);
	llint count(0);
	for (llint k = 0; k < GK.outerSize(); k++)//k is the column id
	{
		jcol[k] = count + 1;
		for (SparseMatrix<double>::InnerIterator it(GK, k); it; ++it)
		{
			if (it.row() >= k)//symmetric
			{
				irow[count] = it.row() + 1;
				val[count] = it.value();
				count++;
			}
		}
	}
	jcol[GK.outerSize()] = nnz + 1;

	int ndgt(16);
	string fname(fn + "_spmat.txt");
	ofstream fout;
	//fout.open(fname.c_str());
	fout.open(fname.c_str(), ios::out | ios::binary);
	if (fout.is_open())
	{
		fout << jcol.size() << " " << irow.size() << " " << val.size() << " " << GF.size() << "\n";
		for (llint i = 0; i < jcol.size(); i++)
		{
			fout << jcol[i] << " ";
			if (i != 0 && i % 100 == 0) fout << '\n';
			//cout << jcol[i] << " ";
			//if (i != 0 && i % 1000 == 0) getchar();
		}
		//cout << "col..."; getchar();
		for (llint i = 0; i < irow.size(); i++)
		{
			fout << irow[i] << " ";
			if (i != 0 && i % 100 == 0) fout << '\n';
		}
		for (llint i = 0; i < val.size(); i++)
		{
			fout << setprecision(ndgt) << val[i] << " ";
			if (i != 0 && i % 100 == 0) fout << '\n';
		}
		for (llint i = 0; i < GF.size(); i++)
		{
			fout << setprecision(ndgt) << GF[i] << " ";
			if (i != 0 && i % 100 == 0) fout << '\n';
		}
		fout.close();
	}
	else
	{
		cerr << "Can't open " << fname << "\n";
	}
}

void Laplace::Write2Matlab_SparseMat(SparseMatrix<double>& GK, VectorXd& GF, string fn)
{
	vector<llint> irow(GK.nonZeros()), jcol(GK.nonZeros());
	vector<double> val(GK.nonZeros());
	llint count(0);
	for (llint k = 0; k < GK.outerSize(); k++)//k is the column id
	{
		//jcol[k] = count + 1;
		for (SparseMatrix<double>::InnerIterator it(GK, k); it; ++it)
		{
			//if (it.row() >= k)//symmetric
			{
				irow[count] = it.row() + 1;
				jcol[count] = it.col() + 1;
				val[count] = it.value();
				count++;
			}
			//cout << "(col,row,val) = " << it.col() << " " << it.row() << " " << it.value() << "\n";
			//getchar();
		}
	}

	int ndgt(16);
	string fname1(fn + "_spmat_i.txt");
	string fname2(fn + "_spmat_j.txt");
	string fname3(fn + "_spmat_v.txt");
	string fname4(fn + "_rhs.txt");
	ofstream fout;
	fout.open(fname1.c_str());
	if (fout.is_open())
	{
		for (llint i = 0; i < irow.size(); i++)
		{
			fout << irow[i] << "\n";
		}
		fout.close();
	}
	else
	{
		cerr << "Can't open " << fname1 << "\n";
		return;
	}

	fout.open(fname2.c_str());
	if (fout.is_open())
	{
		for (llint i = 0; i < jcol.size(); i++)
		{
			fout << jcol[i] << "\n";
		}
		fout.close();
	}
	else
	{
		cerr << "Can't open " << fname2 << "\n";
		return;
	}

	fout.open(fname3.c_str());
	if (fout.is_open())
	{
		for (llint i = 0; i < val.size(); i++)
		{
			fout << setprecision(ndgt) << val[i] << "\n";
		}
		fout.close();
	}
	else
	{
		cerr << "Can't open " << fname3 << "\n";
		return;
	}

	fout.open(fname4.c_str());
	if (fout.is_open())
	{
		for (llint i = 0; i < GF.size(); i++)
		{
			fout << setprecision(ndgt) << GF[i] << "\n";
		}
		fout.close();
	}
	else
	{
		cerr << "Can't open " << fname4 << "\n";
		return;
	}
}

void Laplace::ReadSolution(string fn)
{
	vector<double> sol(neq);
	llint dim;
	string fname(fn + "_sol.txt");
	ifstream fin;
	//fin.open(fname.c_str());
	fin.open(fname.c_str(), ios::in | ios::binary);
	if (fin.is_open())
	{
		//fin >> dim;
		//sol.resize(dim);
		for (llint i = 0; i < neq; i++)
		{
			fin >> sol[i];
		}
		fin.close();
	}
	else
	{
		cerr << "Can't open " << fname << "!\n";
	}

	//if (sol.size() != neq)
	//{
	//	cerr << "Solution dimension wrong!\n";
	//	getchar();
	//}

	uh.resize(npt);
	for (int i = 0; i<npt; i++)
	{
		if (IDBC[i] != -1)
			uh[i] = sol[IDBC[i]];
		else
		{
			uh[i] = gh[i];
			//uh[i] = 0.;
		}
	}

	cout << "Done solving!\n";
}





//for comparison purely B-spline-like
void Laplace::NitscheBoundary_Spline(const BezierElement3D& cell, double beta, vector<vector<double>>& EK, vector<double>& EF)
{
	vector<double> Nx(cell.IEN.size());
	vector<array<double, 3>> dNdx(cell.IEN.size());
	double bzpts[16][3];
	double detJ, Fb;
	array<double, 3> x;
	int i, j;
	double nm[3] = { 0., 0., 0. };
	int ptfloc[6][16] = { { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 }, { 0, 1, 2, 3, 16, 17, 18, 19, 32, 33, 34, 35, 48, 49, 50, 51 },
	{ 3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63 }, { 12, 13, 14, 15, 28, 29, 30, 31, 44, 45, 46, 47, 60, 61, 62, 63 },
	{ 0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60 }, { 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63 } };

	if (cell.bc[0] == 1)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[0][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[0][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[0][ip]][2];
		}
		for (i = 0; i < Gpt.size(); i++)
		{
			for (j = 0; j < Gpt.size(); j++)
			{
				BasisFunction_IGA(Gpt[j], Gpt[i], 0., cell, Nx, dNdx, detJ);
				//BasisFunction_Bezier(Gpt[j], Gpt[i], 0., cell.pts, Nx, dNdx, detJ);
				GetFaceInfo(Gpt[j], Gpt[i], bzpts, detJ, nm);
				detJ = wght[i] * wght[j] * detJ;
				nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
				ElementMatrix_Bezier_Boundary(Nx, dNdx, nm, beta, detJ, EK);
				cell.Para2Phys(Gpt[j], Gpt[i], 0., x.data());
				Fb = exact_sol(x);
				//Fb = cell.bcval[0];
				ElementForce_Bezier_Boundary(Nx, dNdx, nm, beta, detJ, Fb, EF);
			}
		}
	}
	if (cell.bc[1] == 1)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[1][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[1][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[1][ip]][2];
		}
		for (i = 0; i < Gpt.size(); i++)
		{
			for (j = 0; j < Gpt.size(); j++)
			{
				BasisFunction_IGA(Gpt[j], 0., Gpt[i], cell, Nx, dNdx, detJ);
				//BasisFunction_Bezier(Gpt[j], 0., Gpt[i], cell.pts, Nx, dNdx, detJ);
				GetFaceInfo(Gpt[j], Gpt[i], bzpts, detJ, nm);
				detJ = wght[i] * wght[j] * detJ;
				//nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
				ElementMatrix_Bezier_Boundary(Nx, dNdx, nm, beta, detJ, EK);
				cell.Para2Phys(Gpt[j], 0., Gpt[i], x.data());
				Fb = exact_sol(x);
				//Fb = cell.bcval[1];
				ElementForce_Bezier_Boundary(Nx, dNdx, nm, beta, detJ, Fb, EF);
			}
		}
	}
	if (cell.bc[2] == 1)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[2][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[2][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[2][ip]][2];
		}
		for (i = 0; i < Gpt.size(); i++)
		{
			for (j = 0; j < Gpt.size(); j++)
			{
				BasisFunction_IGA(1., Gpt[j], Gpt[i], cell, Nx, dNdx, detJ);
				//BasisFunction_Bezier(1., Gpt[j], Gpt[i], cell.pts, Nx, dNdx, detJ);
				GetFaceInfo(Gpt[j], Gpt[i], bzpts, detJ, nm);
				detJ = wght[i] * wght[j] * detJ;
				//nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
				ElementMatrix_Bezier_Boundary(Nx, dNdx, nm, beta, detJ, EK);
				cell.Para2Phys(1., Gpt[j], Gpt[i], x.data());
				Fb = exact_sol(x);
				//Fb = cell.bcval[2];
				ElementForce_Bezier_Boundary(Nx, dNdx, nm, beta, detJ, Fb, EF);
			}
		}
	}
	if (cell.bc[3] == 1)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[3][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[3][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[3][ip]][2];
		}
		for (i = 0; i < Gpt.size(); i++)
		{
			for (j = 0; j < Gpt.size(); j++)
			{
				BasisFunction_IGA(Gpt[j], 1., Gpt[i], cell, Nx, dNdx, detJ);
				//BasisFunction_Bezier(Gpt[j], 1., Gpt[i], cell.pts, Nx, dNdx, detJ);
				GetFaceInfo(Gpt[j], Gpt[i], bzpts, detJ, nm);
				detJ = wght[i] * wght[j] * detJ;
				nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
				ElementMatrix_Bezier_Boundary(Nx, dNdx, nm, beta, detJ, EK);
				cell.Para2Phys(Gpt[j], 1., Gpt[i], x.data());
				Fb = exact_sol(x);
				//Fb = cell.bcval[3];
				ElementForce_Bezier_Boundary(Nx, dNdx, nm, beta, detJ, Fb, EF);
			}
		}
	}
	if (cell.bc[4] == 1)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[4][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[4][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[4][ip]][2];
		}
		for (i = 0; i < Gpt.size(); i++)
		{
			for (j = 0; j < Gpt.size(); j++)
			{
				BasisFunction_IGA(0., Gpt[j], Gpt[i], cell, Nx, dNdx, detJ);
				//BasisFunction_Bezier(0., Gpt[j], Gpt[i], cell.pts, Nx, dNdx, detJ);
				GetFaceInfo(Gpt[j], Gpt[i], bzpts, detJ, nm);
				detJ = wght[i] * wght[j] * detJ;
				nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
				ElementMatrix_Bezier_Boundary(Nx, dNdx, nm, beta, detJ, EK);
				cell.Para2Phys(0., Gpt[j], Gpt[i], x.data());
				Fb = exact_sol(x);
				//Fb = cell.bcval[4];
				ElementForce_Bezier_Boundary(Nx, dNdx, nm, beta, detJ, Fb, EF);
			}
		}
	}
	if (cell.bc[5] == 1)
	{
		for (int ip = 0; ip < 16; ip++)
		{
			bzpts[ip][0] = cell.pts[ptfloc[5][ip]][0];
			bzpts[ip][1] = cell.pts[ptfloc[5][ip]][1];
			bzpts[ip][2] = cell.pts[ptfloc[5][ip]][2];
		}
		for (i = 0; i < Gpt.size(); i++)
		{
			for (j = 0; j < Gpt.size(); j++)
			{
				BasisFunction_IGA(Gpt[j], Gpt[i], 1., cell, Nx, dNdx, detJ);
				//BasisFunction_Bezier(Gpt[j], Gpt[i], 1., cell.pts, Nx, dNdx, detJ);
				GetFaceInfo(Gpt[j], Gpt[i], bzpts, detJ, nm);
				detJ = wght[i] * wght[j] * detJ;
				//nm[0] = -nm[0]; nm[1] = -nm[1]; nm[2] = -nm[2];
				ElementMatrix_Bezier_Boundary(Nx, dNdx, nm, beta, detJ, EK);
				cell.Para2Phys(Gpt[j], Gpt[i], 1., x.data());
				Fb = exact_sol(x);
				//Fb = cell.bcval[5];
				ElementForce_Bezier_Boundary(Nx, dNdx, nm, beta, detJ, Fb, EF);
			}
		}
	}
}












double Laplace::f_source_lap(const array<double, 3>& x)
{
	//double tmp=f_source_5(x);
	//double tmp = f_source_6(x);
	double tmp = f_source_7_lap(x);
	//double tmp = f_source_8(x);
	//double tmp = f_source_9(x);
	return tmp;
}

//double Laplace::f_source_1(const array<double, 3>& x)
//{
//	double nm = (range[0][1] - range[0][0])*(range[0][1] - range[0][0])*
//		(range[1][1] - range[1][0])*(range[1][1] - range[1][0])*
//		(range[2][1] - range[2][0])*(range[2][1] - range[2][0]);
//	double tmp = -2.*(x[1] - range[1][0])* (range[1][1] - x[1])*(x[2] - range[2][0]) * (range[2][1] - x[2])
//		- 2.*(x[0] - range[0][0])* (range[0][1] - x[0])*(x[2] - range[2][0]) * (range[2][1] - x[2])
//		- 2.*(x[1] - range[1][0])* (range[1][1] - x[1])*(x[0] - range[0][0]) * (range[0][1] - x[0]);
//	return -tmp/nm;
//}
//
//double Laplace::f_source_2(const array<double, 3>& x)//analysis domain is a unit cube
//{
//	double tmp1 = -1600.*exp(-654. + 800.*x[0] * (1. - x[0]) + 800.*x[1] * (1. - x[1]) + 800.*x[2] * (1. - x[2]))/3.;
//	double tmp2 = exp((7. - 20.*x[0])*(7. - 20.*x[0]) + (7. - 20.*x[1])*(7. - 20.*x[1]) + (7. - 20.*x[2])*(7. - 20.*x[2]));
//	double tmp3 = exp((13. - 20.*x[0])*(13. - 20.*x[0]) + (13. - 20.*x[1])*(13. - 20.*x[1]) + (13. - 20.*x[2])*(13. - 20.*x[2]));
//	double tmp4 = 1011. - 1040.*x[0] + 800.*x[0] * x[0] - 1040.*x[1] + 800.*x[1] * x[1] - 1040.*x[2] + 800.*x[2] * x[2];
//	double tmp5 = 291. - 560.*x[0] + 800.*x[0] * x[0] - 560.*x[1] + 800.*x[1] * x[1] - 560.*x[2] + 800.*x[2] * x[2];
//	return tmp1*(tmp2*tmp4+tmp3*tmp5);
//}
//
//double Laplace::f_source_3(const array<double, 3>& x)//analysis domain is a unit cube
//{
//	double tmp1 = -1600.*exp(-300. + 400.*x[0] * (1. - x[0]) + 400.*x[1] * (1. - x[1]) + 400.*x[2] * (1. - x[2])) / 3.;
//	double tmp2 = 597. - 800.*x[0] * (1. - x[0]) - 800.*x[1] * (1. - x[1]) - 800.*x[2] * (1. - x[2]);
//	return tmp1*tmp2;
//}
//
//double Laplace::f_source_4(const array<double, 3>& x)//analysis domain is a unit cube
//{
//	double tmp0 = 1. + 20.4124*x[0] + 20.4124*x[1] - 40.8248*x[2];
//	double tmp1 = 1. / cosh(tmp0);
//	return 5000.*tmp1*tmp1*tanh(tmp0);
//}
//
//double Laplace::f_source_5(const array<double, 3>& x)//analysis domain is arbitrary
//{
//	double tmp = -10.*(x[0] + x[1] + x[2]);
//	return tmp;
//}
//
//double Laplace::f_source_6(const array<double, 3>& x)//analysis domain is arbitrary
//{
//	double tmp = 2.*(x[0] * (1. - x[0])*x[1] * (1. - x[1]) + x[0] * (1. - x[0])*x[2] * (1. - x[2]) + x[2] * (1. - x[2])*x[1] * (1. - x[1]));
//	return tmp;
//}
//
double Laplace::f_source_7_lap(const array<double, 3>& x)
{
	//double tmp0 = acoef*(nmpl[0] * (x[0] - range[0][0]) + nmpl[1] * (x[1] - range[1][0]) + nmpl[2] * (x[2] - range[2][0]));
	double tmp0 = acoef*(nmpl[0] * (x[0] - xorg[0]) + nmpl[1] * (x[1] - xorg[1]) + nmpl[2] * (x[2] - xorg[2]));
	double tmp1 = 1. / cosh(tmp0);
	return bcoef*tmp1*tmp1*tanh(tmp0);
}

//double Laplace::f_source_8(const array<double, 3>& x)
//{
//	double x1[3] = { (x[0] - range[0][0]) / dmlen[0], (x[1] - range[1][0]) / dmlen[1], (x[2] - range[2][0]) / dmlen[2] };
//	double len2[3] = { dmlen[0] * dmlen[0], dmlen[1] * dmlen[1], dmlen[2] * dmlen[2] };
//	double tmp0 = -acoef / exp(acoef*((x1[0] - nmpl[0])*(x1[0] - nmpl[0]) + (x1[1] - nmpl[1])*(x1[1] - nmpl[1]) + (x1[2] - nmpl[2])*(x1[2] - nmpl[2])));
//	double tmp1 = (-2. + acoef)*(1. / len2[0] + 1. / len2[1] + 1. / len2[2]);
//	double tmp2 = 4.*acoef*(x1[0] * (x1[0] - 1.) / len2[0] + x1[1] * (x1[1] - 1.) / len2[1] + x1[2] * (x1[2] - 1.) / len2[2]);
//	return tmp0*(tmp1+tmp2);
//}
//
//double Laplace::f_source_9(const array<double, 3>& x)
//{
//	double PI(3.141592654);
//	return 3.*PI*PI*sin(PI*x[0]) * sin(PI*x[1]) * sin(PI*x[2]);
//}
//
double Laplace::exact_sol_lap(const array<double, 3>& x)
{
	//double tmp = exact_sol_5(x);
	//double tmp = exact_sol_6(x);
	double tmp = exact_sol_7_lap(x);
	//double tmp = exact_sol_8_lap(x);
	//double tmp = exact_sol_9(x);
	return tmp;
}

//double Laplace::exact_sol_1(const array<double, 3>& x)
//{
//	double nm = (range[0][1] - range[0][0])*(range[0][1] - range[0][0])*
//		(range[1][1] - range[1][0])*(range[1][1] - range[1][0])*
//		(range[2][1] - range[2][0])*(range[2][1] - range[2][0]);
//	return (x[0] - range[0][0])*(range[0][1] - x[0])*(x[1] - range[1][0])*(range[1][1] - x[1])*(x[2] - range[2][0])*(range[2][1] - x[2])/nm;
//}
//
//double Laplace::exact_sol_2(const array<double, 3>& x)
//{
//	double tmp1 = 1. / exp((20.*x[0] - 7.)*(20.*x[0] - 7.) + (20.*x[1] - 7.)*(20.*x[1] - 7.) + (20.*x[2] - 7.)*(20.*x[2] - 7.));
//	double tmp2 = 1. / exp((20.*x[0] - 13.)*(20.*x[0] - 13.) + (20.*x[1] - 13.)*(20.*x[1] - 13.) + (20.*x[2] - 13.)*(20.*x[2] - 13.));
//	return 2.*(tmp1+tmp2)/3.;
//}
//
//double Laplace::exact_sol_3(const array<double, 3>& x)
//{
//	return 2. / (3.*exp((20.*x[0] - 10.)*(20.*x[0] - 10.) + (20.*x[1] - 10.)*(20.*x[1] - 10.) + (20.*x[2] - 10.)*(20.*x[2] - 10.)));
//}
//
//double Laplace::exact_sol_4(const array<double, 3>& x)
//{
//	return tanh(1. - 50.*(-0.408248*x[0]-0.408248*x[1] +0.816497*x[2]));
//}
//
//double Laplace::exact_sol_5(const array<double, 3>& x)
//{
//	double tmp = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2])*(x[0] + x[1] + x[2]) + x[0] * x[1] * x[2];
//	return tmp;
//}
//
//double Laplace::exact_sol_6(const array<double, 3>& x)
//{
//	double tmp = x[0] * (1. - x[0])*x[1] * (1. - x[1])*x[2] * (1. - x[2]);;
//	return tmp;
//}
//
double Laplace::exact_sol_7_lap(const array<double, 3>& x)
{
	//double x1[3] = { (x[0] - range[0][0]), (x[1] - range[1][0]), (x[2] - range[2][0]) };
	double x1[3] = { (x[0] - xorg[0]), (x[1] - xorg[1]), (x[2] - xorg[2]) };
	return tanh(acoef*(nmpl[0] * x1[0] + nmpl[1] * x1[1] + nmpl[2] * x1[2]));
}

//double Laplace::exact_sol_8(const array<double, 3>& x)
//{
//	double x1[3] = { (x[0] - range[0][0]) / dmlen[0], (x[1] - range[1][0]) / dmlen[1], (x[2] - range[2][0]) / dmlen[2] };
//	double tmp = 1. / exp(acoef*((x1[0] - nmpl[0])*(x1[0] - nmpl[0]) + (x1[1] - nmpl[1])*(x1[1] - nmpl[1]) + (x1[2] - nmpl[2])*(x1[2] - nmpl[2])));
//	return 1. / exp(acoef*((x1[0] - nmpl[0])*(x1[0] - nmpl[0]) + (x1[1] - nmpl[1])*(x1[1] - nmpl[1]) + (x1[2] - nmpl[2])*(x1[2] - nmpl[2]) ));
//}
//
//double Laplace::exact_sol_9(const array<double, 3>& x)
//{
//	double PI(3.141592654);
//	return sin(PI*x[0]) * sin(PI*x[1]) * sin(PI*x[2]);
//}

void Laplace::GetRemoveRegion(double xy[3][2])
{
	rmv[0][0] = xy[0][0]; rmv[0][1] = xy[0][1];
	rmv[1][0] = xy[1][0]; rmv[1][1] = xy[1][1];
	rmv[2][0] = xy[2][0]; rmv[2][1] = xy[2][1];
}

void Laplace::GetEqParameter(double xy[3][2], double nrm[3], double a)
{
	range[0][0] = xy[0][0]; range[0][1] = xy[0][1];
	range[1][0] = xy[1][0]; range[1][1] = xy[1][1];
	range[2][0] = xy[2][0]; range[2][1] = xy[2][1];
	nmpl[0] = nrm[0]; nmpl[1] = nrm[1]; nmpl[2] = nrm[2];
	acoef = a;
	dmlen[0] = range[0][1] - range[0][0];
	dmlen[1] = range[1][1] - range[1][0];
	dmlen[2] = range[2][1] - range[2][0];
	double nm2[3] = { nmpl[0] * nmpl[0], nmpl[1] * nmpl[1], nmpl[2] * nmpl[2] };
	double len2[3] = { dmlen[0] * dmlen[0], dmlen[1] * dmlen[1], dmlen[2] * dmlen[2] };
	//bcoef = 2.*acoef*acoef*(nm2[2] * len2[0] * len2[1] + (nm2[1] * len2[0] + nm2[0] * len2[1])*len2[2]) / (len2[0] * len2[1] * len2[2]);
	bcoef = 2.*acoef*acoef*(nm2[0] + nm2[1] + nm2[2]);
	double a1(0.7);//rod
	double a2(1. - a1);
	xorg[0] = a1*range[0][0] + a2*range[0][1];
	xorg[1] = a1*range[1][0] + a2*range[1][1];
	xorg[2] = a1*range[2][0] + a2*range[2][1];
}

void Laplace::PipelineTmp(const vector<BezierElement3D>& bzmesh, string fn)
{
	uh.resize(npt);
	string fname("../io/pipeline/base_disp.vtk"), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i<4; i++) getline(fin, stmp);//skip lines
		fin >> stmp >> npts >> stmp;
		for (int i = 0; i<npts+1; i++)
		{
			getline(fin, stmp);
		}
		fin >> stmp >> neles >> itmp;
		for (int i = 0; i<neles+1; i++)
		{
			getline(fin, stmp);
		}
		for (int i = 0; i<neles + 5; i++) getline(fin, stmp);
		for (int i = 0; i < npts; i++) fin >> uh[i];
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}
	VisualizeVTK_1(bzmesh, fn);//with smoother boundary
	//VisualizeVTK_2(bzmesh, fn);//smooth boundary, remove elements
}