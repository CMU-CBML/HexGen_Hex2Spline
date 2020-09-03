#include "LeastSquare.h"
#include "PDEdata.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <omp.h>

using namespace std;

typedef unsigned int uint;

LeastSquare::LeastSquare()
{
	range[0][0] = 0.; range[0][1] = 1.;
	range[1][0] = 0.; range[1][1] = 1.;
	range[2][0] = 0.; range[2][1] = 1.;
}

void LeastSquare::GaussInfo(int ng)
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

void LeastSquare::SetProblem(const vector<int>& IDBC_in, const vector<double>& gh_in)
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
}

void LeastSquare::BasisFunction(double u, double v, double w, const vector<array<double, 3>>& pt, double Nx[64], double dNdx[64][3], double& detJ)
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
				//dNdt[loc][0] = dNdu[k] * Nv[j] * Nw[i];
				//dNdt[loc][1] = Nu[k] * dNdv[j] * Nw[i];
				//dNdt[loc][2] = Nu[k] * Nv[j] * dNdw[i];
				loc++;
			}
		}
	}
	//Matrix3d dxdt = Matrix3d::Zero();
	//loc = 0;
	//for (i = 0; i<4; i++)
	//{
	//	for (j = 0; j<4; j++)
	//	{
	//		for (k = 0; k < 4; k++)
	//		{
	//			for (a = 0; a<3; a++)
	//			{
	//				for (b = 0; b<3; b++)
	//				{
	//					dxdt(a, b) += pt[loc][a] * dNdt[loc][b];
	//				}
	//			}
	//			loc++;
	//		}
	//	}
	//}
	////cout << dxdt << "\n\n";
	//Matrix3d dtdx = dxdt.inverse();
	//for (i = 0; i<64; i++)
	//{
	//	dNdx[i][0] = dNdt[i][0] * dtdx(0, 0) + dNdt[i][1] * dtdx(1, 0) + dNdt[i][2] * dtdx(2, 0);
	//	dNdx[i][1] = dNdt[i][0] * dtdx(0, 1) + dNdt[i][1] * dtdx(1, 1) + dNdt[i][2] * dtdx(2, 1);
	//	dNdx[i][2] = dNdt[i][0] * dtdx(0, 2) + dNdt[i][1] * dtdx(1, 2) + dNdt[i][2] * dtdx(2, 2);
	//}
	//detJ = dxdt.determinant();
	//detJ = 0.25*detJ;
}

void LeastSquare::BasisFunctionMap(double Nx0[64], const vector<vector<double>>& cmat, vector<double>& Nx)
{
	//Nx.resize(cmat.size(), 0.);
	for (uint i = 0; i < cmat.size(); i++)
	{
		Nx[i] = 0.;
		for (int j = 0; j < 64; j++)
		{
			Nx[i] += cmat[i][j] * Nx0[j];
		}
	}
}

void LeastSquare::ElementMatrix(double Nx[64], double detJ, double EK[64][64])
{
	int i, j;
	for (i = 0; i<64; i++)
	{
		for (j = 0; j<64; j++)
		{
			EK[i][j] += Nx[i] * Nx[j];
		}
	}
}

void LeastSquare::ElementMatrix(const vector<double>& Nx, vector<vector<double>>& EK)
{
	int i, j;
	for (i = 0; i < Nx.size(); i++)
	{
		for (j = 0; j < Nx.size(); j++)
		{
			EK[i][j] += Nx[i] * Nx[j];
		}
	}
}

void LeastSquare::ElementMatrix_Laplace(double dNdx[64][3], double detJ, double EK[64][64])
{
	int i, j;
	for (i = 0; i<64; i++)
	{
		for (j = 0; j<64; j++)
		{
			EK[i][j] += (dNdx[i][0] * dNdx[j][0] + dNdx[i][1] * dNdx[j][1] + dNdx[i][2] * dNdx[j][2])*detJ;
		}
	}
}

void LeastSquare::ElementForce(double Nx[64], double detJ, double Fb, double EF[64])
{
	for (int i = 0; i<64; i++)
	{
		EF[i] += Nx[i] * Fb;
	}
}

void LeastSquare::ElementForce(const vector<double>& Nx, double Fb, vector<double>& EF)
{
	for (int i = 0; i < Nx.size(); i++)
	{
		EF[i] += Nx[i] * Fb;
	}
}

void LeastSquare::Assembly(double EK1[64][64], double EF1[64], const vector<vector<double>>& cmat, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF)
{
	int i, j, A, B;
	vector<vector<double>> EK(IEN.size(), vector<double>(IEN.size(), 0.));
	vector<double> EF(IEN.size(), 0.);
	//cout << "before mapping\n";

#pragma omp parallel for
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
	//cout << "after mapping\n";
	//getchar();
#pragma omp parallel for
	for (i = 0; i < IEN.size(); i++)
	{
		for (A = 0; A < 64; A++)
		{
			EF[i] += cmat[i][A] * EF1[A];
		}
	}
	//cout << "before global\n";
	//getchar();
//#pragma omp parallel for
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
			//cout << IDBC[IEN[i]] << " " << IDBC[IEN[j]] << "\n";
			//getchar();
			if (IDBC[IEN[i]] != -1 && IDBC[IEN[j]] != -1)
			{
				GK.coeffRef(IDBC[IEN[i]], IDBC[IEN[j]]) += EK[i][j];
			}
			else if (IDBC[IEN[i]] != -1 && IDBC[IEN[j]] == -1)
			{
				GF(IDBC[IEN[i]]) -= EK[i][j] * gh[IEN[j]];
			}
		}
	}
	//cout << "after global\n";
	//getchar();
#pragma omp parallel for
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

void LeastSquare::Assembly(const vector<vector<double>>& EK, const vector<double>& EF, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF)
{
	int i, j, A, B;
	for (i = 0; i<IEN.size(); i++)
	{
		if (IEN[i] != -1)
		{
			for (j = 0; j<IEN.size(); j++)
			{
				if (IEN[j] != -1)
				{
					if (IDBC[IEN[i]] != -1 && IDBC[IEN[j]] != -1)
					{
						GK.coeffRef(IDBC[IEN[i]], IDBC[IEN[j]]) += EK[i][j];
					}
					//else if (IDBC[IEN[i]] != -1 && IDBC[IEN[j]] == -1)
					//{
					//	GF(IDBC[IEN[i]]) -= EK[i][j] * gh[IEN[j]];
					//}
				}
			}
		}
	}
	for (i = 0; i < IEN.size(); i++)
	{
		if (IEN[i] != -1)
		{
			if (IDBC[IEN[i]] != -1)
			{
				GF(IDBC[IEN[i]]) += EF[i];
			}
		}
	}
}

void LeastSquare::Assembly_Fitting(double EK1[64][64], double EF1[64], const vector<vector<double>>& cmat, const vector<int>& IEN, SparseMatrix<double>& GK, VectorXd& GF)
{
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
	//cout << "after mapping\n";
	//getchar();
//#pragma omp parallel for
	for (i = 0; i < IEN.size(); i++)
	{
		for (A = 0; A < 64; A++)
		{
			EF[i] += cmat[i][A] * EF1[A];
		}
	}
	//cout << "before global\n";
	//getchar();
	//#pragma omp parallel for
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
			//cout << IDBC[IEN[i]] << " " << IDBC[IEN[j]] << "\n";
			//getchar();
			if (IDBC[IEN[i]] != -1 && IDBC[IEN[j]] != -1)
			{
				GK.coeffRef(IDBC[IEN[i]], IDBC[IEN[j]]) += EK[i][j];
			}
		}
	}
	//cout << "after global\n";
	//getchar();
//#pragma omp parallel for
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

void LeastSquare::BuildLinearSystem(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF)
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
	for (e = 0; e<bzmesh.size(); e++)
	{
		//cout << "element id: " << e << "/" << bzmesh.size() << "\n";
		if (e!=0 && e % 500 == 0)
		{
			cout << e <<" ";
		}
		for (i = 0; i<64; i++)
		{
			EF[i] = 0.;
			for (j = 0; j<64; j++)
			{
				EK[i][j] = 0.;
			}
		}
		//cout << "before loop\n";
		for (i = 0; i<Gpt.size(); i++)
		{
			for (j = 0; j<Gpt.size(); j++)
			{
				for (k = 0; k < Gpt.size(); k++)
				{
					BasisFunction(Gpt[i], Gpt[j], Gpt[k], bzmesh[e].pts, Nx, dNdx, detJ);
					detJ = wght[i] * wght[j] * wght[k] * detJ;
					//ElementMatrix(Nx, detJ, EK);
					ElementMatrix_Laplace(dNdx, detJ, EK);

					//Para2Phys(Gpt[i], Gpt[j], Gpt[k], bzmesh[e].pts, Nx, pt);
					//Fb = f_source_1(pt);
					//Fb = f_source_2(pt);
					//Fb = f_source_3(pt);
					//Fb = f_source_4(pt);
					//ElementForce(Nx, detJ, Fb, EF);
				}
			}
		}
		Assembly(EK, EF, bzmesh[e].cmat, bzmesh[e].IEN, GK, GF);
	}
}

void LeastSquare::BuildLinearSystem_Fitting(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF)
{
	//int e, i, j, k, a, b;
	//double EK[64][64];
	//double EF[64];
	//double Nx[64];
	//double dNdx[64][3];
	//double detJ(1.), Fb(0.);
	//array<double, 3> pt;

	int nresv(5 * 64);
	if (nresv < neq)
	{
		GK.reserve(VectorXi::Constant(neq, nresv));
	}
	else
	{
		GK.reserve(VectorXi::Constant(neq, neq));
	}
	int nsp(5+1);
	double rng[2] = { 0., 1. };
	vector<double> spt(nsp);
	for (int i = 0; i < nsp; i++)
	{
		spt[i] = rng[0] + double(i) * (rng[1] - rng[0]) / double(nsp - 1);
	}

	cout << "npt: " << npt << "\n";
	cout << "neq: " << neq << "\n";
	cout << "nel: " << bzmesh.size() << '\n';

//#pragma omp parallel for
	for (int e = 0; e<bzmesh.size(); e++)
	{
		if (e != 0 && e % 1000 == 0)
		{
			cout << e << " ";
		}
		//double EK[64][64];
		//double EF[64];
		//double Nx[64];
		//double dNdx[64][3];
		//double detJ(1.), Fb(0.);
		//array<double, 3> pt;
		for (int k = 0; k < bzmesh[e].bfc.size(); k++)
		{
			double Nx[64];
			double dNdx[64][3];
			vector<double> Nx1(bzmesh[e].IEN.size());
			vector<vector<double>> EK(bzmesh[e].IEN.size(), vector<double>(bzmesh[e].IEN.size(), 0.));
			vector<double> EF(bzmesh[e].IEN.size(), 0.);
			double detJ(1.), Fb(0.);
			array<double, 3> pt;
			//for (int i = 0; i<64; i++)
			//{
			//	EF[i] = 0.;
			//	for (int j = 0; j<64; j++)
			//	{
			//		EK[i][j] = 0.;
			//	}
			//}

			for (int i = 0; i<spt.size(); i++)
			{
				for (int j = 0; j<spt.size(); j++)
				{
					if (bzmesh[e].bfc[k] == 0)
					{
						BasisFunction(spt[i], spt[j], 0., bzmesh[e].pts, Nx, dNdx, detJ);
						Para2Phys(spt[i], spt[j], 0., bzmesh[e].pts, Nx, pt);
					}
					else if (bzmesh[e].bfc[k] == 1)
					{
						BasisFunction(spt[i], 0., spt[j], bzmesh[e].pts, Nx, dNdx, detJ);
						Para2Phys(spt[i], 0., spt[j], bzmesh[e].pts, Nx, pt);
					}
					else if (bzmesh[e].bfc[k] == 2)
					{
						BasisFunction(1., spt[i], spt[j], bzmesh[e].pts, Nx, dNdx, detJ);
						Para2Phys(1., spt[i], spt[j], bzmesh[e].pts, Nx, pt);
					}
					else if (bzmesh[e].bfc[k] == 3)
					{
						BasisFunction(spt[i], 1., spt[j], bzmesh[e].pts, Nx, dNdx, detJ);
						Para2Phys(spt[i], 1., spt[j], bzmesh[e].pts, Nx, pt);
					}
					else if (bzmesh[e].bfc[k] == 4)
					{
						BasisFunction(0., spt[i], spt[j], bzmesh[e].pts, Nx, dNdx, detJ);
						Para2Phys(0., spt[i], spt[j], bzmesh[e].pts, Nx, pt);
					}
					else
					{
						BasisFunction(spt[i], spt[j], 1., bzmesh[e].pts, Nx, dNdx, detJ);
						Para2Phys(spt[i], spt[j], 1., bzmesh[e].pts, Nx, pt);
					}
					BasisFunctionMap(Nx, bzmesh[e].cmat, Nx1);
					ElementMatrix(Nx1, EK);
					Fb = exact_sol_ls(pt);
					ElementForce(Nx1, Fb, EF);
				}
			}
			Assembly(EK, EF, bzmesh[e].IEN, GK, GF);
			//Assembly_Fitting(EK, EF, bzmesh[e].cmat, bzmesh[e].IEN, GK, GF);
		}
	}
}

void LeastSquare::BuildLinearSystem_Fitting_C012(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK, VectorXd& GF)
{
	int nresv(5 * 64);
	if (nresv < neq)
	{
		GK.reserve(VectorXi::Constant(neq, nresv));
	}
	else
	{
		GK.reserve(VectorXi::Constant(neq, neq));
	}
	int nsp(5 + 1);
	double rng[2] = { 0., 1. };
	vector<double> spt(nsp);
	for (int i = 0; i < nsp; i++)
	{
		spt[i] = rng[0] + double(i) * (rng[1] - rng[0]) / double(nsp - 1);
	}

	cout << "npt: " << npt << "\n";
	cout << "neq: " << neq << "\n";
	cout << "nel: " << bzmesh.size() << '\n';

	//#pragma omp parallel for
	for (int e = 0; e<bzmesh.size(); e++)
	{
		if (e != 0 && e % 1000 == 0)
		{
			cout << e << " ";
		}
		if (bzmesh[e].type == 1)
		{
			double Nx[64];
			double dNdx[64][3];
			vector<double> Nx1(bzmesh[e].IEN.size());
			vector<vector<double>> EK(bzmesh[e].IEN.size(), vector<double>(bzmesh[e].IEN.size(), 0.));
			vector<double> EF(bzmesh[e].IEN.size(), 0.);
			double detJ(1.), Fb(0.);
			array<double, 3> pt;
			uint i, j, k;
			for (k = 0; k < 6; k++)
			{
				if (bzmesh[e].bc[k] == 1)
				{
					for (i = 0; i<spt.size(); i++)
					{
						for (j = 0; j<spt.size(); j++)
						{
							if (k == 0)
							{
								BasisFunction(spt[i], spt[j], 0., bzmesh[e].pts, Nx, dNdx, detJ);
								Para2Phys(spt[i], spt[j], 0., bzmesh[e].pts, Nx, pt);
							}
							else if (k == 1)
							{
								BasisFunction(spt[i], 0., spt[j], bzmesh[e].pts, Nx, dNdx, detJ);
								Para2Phys(spt[i], 0., spt[j], bzmesh[e].pts, Nx, pt);
							}
							else if (k == 2)
							{
								BasisFunction(1., spt[i], spt[j], bzmesh[e].pts, Nx, dNdx, detJ);
								Para2Phys(1., spt[i], spt[j], bzmesh[e].pts, Nx, pt);
							}
							else if (k == 3)
							{
								BasisFunction(spt[i], 1., spt[j], bzmesh[e].pts, Nx, dNdx, detJ);
								Para2Phys(spt[i], 1., spt[j], bzmesh[e].pts, Nx, pt);
							}
							else if (k == 4)
							{
								BasisFunction(0., spt[i], spt[j], bzmesh[e].pts, Nx, dNdx, detJ);
								Para2Phys(0., spt[i], spt[j], bzmesh[e].pts, Nx, pt);
							}
							else
							{
								BasisFunction(spt[i], spt[j], 1., bzmesh[e].pts, Nx, dNdx, detJ);
								Para2Phys(spt[i], spt[j], 1., bzmesh[e].pts, Nx, pt);
							}
							BasisFunctionMap(Nx, bzmesh[e].cmat, Nx1);
							ElementMatrix(Nx1, EK);
							Fb = exact_sol(pt);
							ElementForce(Nx1, Fb, EF);
						}
					}
				}
			}
			Assembly(EK, EF, bzmesh[e].IEN, GK, GF);
		}
	}
}

void LeastSquare::Solver(SparseMatrix<double>& GK, VectorXd& GF)
{
	SimplicialLDLT<SparseMatrix<double>> solver;
	VectorXd sol = solver.compute(GK).solve(GF);
	//cout<<sol;
	//getchar();
	uh.resize(npt);
	for(int i=0; i<npt; i++)
	{
		if(IDBC[i]!=-1)
			uh[i]=sol(IDBC[i]);
		else
			uh[i]=gh[i];
	}

	//for (uint i = 0; i < sol.size(); i++)
	//{
	//	if (fabs(sol(i)-1.) > .1)
	//	{
	//		cout << sol(i) << "\n";
	//		getchar();
	//	}
	//}
	//cout << "done output sol!\n";
}

void LeastSquare::Para2Phys(double u, double v, double w, const vector<array<double, 3>>& bzpt, double Nx[64], array<double, 3>& pt)
{
	pt[0] = 0.; pt[1] = 0.; pt[2] = 0.;
	for (int i = 0; i < 64; i++)
	{
		pt[0] += bzpt[i][0] * Nx[i];
		pt[1] += bzpt[i][1] * Nx[i];
		pt[2] += bzpt[i][2] * Nx[i];
	}
}

void LeastSquare::DispCal(double u, double v, double w, const BezierElement3D& bzel, double& disp, double& detJ)
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

void LeastSquare::ElementError(const BezierElement3D& bzel, double& L2, double& H1)
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
				//sol_e = exact_sol_1(pt);
				//sol_e = exact_sol_2(pt);
				//sol_e = exact_sol_3(pt);
				sol_e = exact_sol_4(pt);
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

void LeastSquare::ErrorEstimate(const vector<BezierElement3D>& bzmesh, vector<double>& err)
{
	double L2, H1;
	err.clear();
	err.resize(bzmesh.size(),0.);
	for (uint e = 0; e < bzmesh.size(); e++)
	{
		ElementError(bzmesh[e],L2,H1);
		err[e] = L2;
	}
}

double LeastSquare::ElementErrorEstimate(const BezierElement3D& bzel)
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

void LeastSquare::BasisFunctionMore(double u, double v, double w, const BezierElement3D& bzel, double pt[3], vector<double>& Nx, vector<array<double, 3>>& dNdx, double& detJ, double dtdx1[3][3])
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
	detJ = 0.25*detJ;
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

void LeastSquare::BubbleFunction(double u, double v, double w, double dtdx[3][3], double& f, double dfdx[3])
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

void LeastSquare::DispFunction(const vector<int>& IEN, const vector<array<double, 3>>& dNdx, double dfdx[3])
{
	dfdx[0] = 0.; dfdx[1] = 0.; dfdx[2] = 0.;
	for (uint i = 0; i < IEN.size(); i++)
	{
		dfdx[0] += uh[IEN[i]] * dNdx[i][0];
		dfdx[1] += uh[IEN[i]] * dNdx[i][1];
		dfdx[2] += uh[IEN[i]] * dNdx[i][2];
	}
}

void LeastSquare::TotalErrorEstimate(const vector<BezierElement3D>& bzmesh, vector<double>& err)
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

void LeastSquare::VisualizeVTK(const vector<BezierElement3D>& bzmesh, string fn)
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
	int ns(2), ecount(0);
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
					if (a == 0 || a == ns - 1 || b == 0 || b == ns - 1)
					{
						lpt.push_back(pt);
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
					el[0] = ecount*nns[0] + a*nns[1] + b*ns + c;
					el[1] = ecount*nns[0] + a*nns[1] + b*ns + c + 1;
					el[2] = ecount*nns[0] + a*nns[1] + (b + 1)*ns + c + 1;
					el[3] = ecount*nns[0] + a*nns[1] + (b + 1)*ns + c;
					el[4] = ecount*nns[0] + (a + 1)*nns[1] + b*ns + c;
					el[5] = ecount*nns[0] + (a + 1)*nns[1] + b*ns + c + 1;
					el[6] = ecount*nns[0] + (a + 1)*nns[1] + (b + 1)*ns + c + 1;
					el[7] = ecount*nns[0] + (a + 1)*nns[1] + (b + 1)*ns + c;
					sele.push_back(el);
				}
			}
		}
		/*for(int a=0;a<ns-1;a++)
		{
		array<int,2> lc;
		lc[0]=ecount*4*(ns-1)+a;
		lc[1]=ecount*4*(ns-1)+a+1;
		led.push_back(lc);
		lc[0]=ecount*4*(ns-1)+3*ns-4+a;
		lc[1]=ecount*4*(ns-1)+3*ns-4+a+1;
		led.push_back(lc);
		}
		for(int a=0;a<ns-2;a++)
		{
		array<int,2> lc;
		lc[0]=ecount*4*(ns-1)+ns+2*a;
		lc[1]=ecount*4*(ns-1)+ns+2*a+2;
		led.push_back(lc);
		lc[0]=ecount*4*(ns-1)+ns+2*a-1;
		lc[1]=ecount*4*(ns-1)+ns+2*a+1;
		led.push_back(lc);
		}
		array<int,2> lc1;
		lc1[0]=ecount*4*(ns-1);
		lc1[1]=ecount*4*(ns-1)+ns;
		led.push_back(lc1);
		lc1[0]=ecount*4*(ns-1)+3*ns-5;
		lc1[1]=ecount*4*(ns-1)+4*ns-5;
		led.push_back(lc1);*/
		ecount++;
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

		fout<<"POINT_DATA "<<sdisp.size()<<"\nSCALARS err float 1\nLOOKUP_TABLE default\n";
		for(uint i=0;i<sdisp.size();i++)
		{
			fout<<sdisp[i]<<"\n";
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

void LeastSquare::VisualizeVTK_1(const vector<BezierElement3D>& bzmesh, string fn)
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

void LeastSquare::VisualizeVTK_2(const vector<BezierElement3D>& bzmesh, string fn)
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

void LeastSquare::VisualizeBoundarySurface(const vector<BezierElement3D>& bzmesh, const vector<array<double, 3>>& cpts, string fn)
{
	vector<array<double, 3>> spt;//s means sample
	vector<double> sdisp;
	vector<array<int, 4>> sele;
	int ns(2), ecount(0);
	vector<double> su(ns);
	double Nx[64], dNdx[64][3], detJ;
	for (int i = 0; i<ns; i++)
	{
		su[i] = double(i) / (double(ns) - 1.);
	}

	for (unsigned int e = 0; e<bzmesh.size(); e++)
	{
		for (uint k = 0; k < bzmesh[e].bfc.size(); k++)
		{
			for (int a = 0; a<ns; a++)
			{
				for (int b = 0; b<ns; b++)
				{
					array<double, 3> pt1 = { 0., 0., 0. };
					if (bzmesh[e].bfc[k] == 0)
					{
						BasisFunction(su[b], su[a], 0., bzmesh[e].pts, Nx, dNdx, detJ);
						Para2Phys(su[b], su[a], 0., bzmesh[e].pts, Nx, pt1);
					}
					else if (bzmesh[e].bfc[k] == 1)
					{
						BasisFunction(su[b], 0., su[a], bzmesh[e].pts, Nx, dNdx, detJ);
						Para2Phys(su[b], 0., su[a], bzmesh[e].pts, Nx, pt1);
					}
					else if (bzmesh[e].bfc[k] == 2)
					{
						BasisFunction(1., su[b], su[a], bzmesh[e].pts, Nx, dNdx, detJ);
						Para2Phys(1., su[b], su[a], bzmesh[e].pts, Nx, pt1);
					}
					else if (bzmesh[e].bfc[k] == 3)
					{
						BasisFunction(su[b], 1., su[a], bzmesh[e].pts, Nx, dNdx, detJ);
						Para2Phys(su[b], 1., su[a], bzmesh[e].pts, Nx, pt1);
					}
					else if (bzmesh[e].bfc[k] == 4)
					{
						BasisFunction(0., su[b], su[a], bzmesh[e].pts, Nx, dNdx, detJ);
						Para2Phys(0., su[b], su[a], bzmesh[e].pts, Nx, pt1);
					}
					else
					{
						BasisFunction(su[b], su[a], 1., bzmesh[e].pts, Nx, dNdx, detJ);
						Para2Phys(su[b], su[a], 1., bzmesh[e].pts, Nx, pt1);
					}
					vector<double> Nx1(bzmesh[e].IEN.size(), 0.);
					array<double, 3> pt = { 0., 0., 0. };
					double u(0.);
					for (uint i = 0; i < bzmesh[e].IEN.size(); i++)
					{
						if (IDBC[bzmesh[e].IEN[i]] != -1)
						{
							for (int j = 0; j < 64; j++)
							{
								Nx1[i] += bzmesh[e].cmat[i][j] * Nx[j];
							}
							pt[0] += Nx1[i] * cpts[bzmesh[e].IEN[i]][0];
							pt[1] += Nx1[i] * cpts[bzmesh[e].IEN[i]][1];
							pt[2] += Nx1[i] * cpts[bzmesh[e].IEN[i]][2];
							u += Nx1[i] * uh[bzmesh[e].IEN[i]];
						}
					}
					spt.push_back(pt);
					//u = exact_sol(pt);
					sdisp.push_back(u);
				}
			}
			//elements
			if (bzmesh[e].bfc[k] == 0)
			{
				array<int, 4> el = { spt.size() - 4, spt.size() - 2, spt.size() - 1, spt.size() - 3 };
				sele.push_back(el);
			}
			else if (bzmesh[e].bfc[k] == 1)
			{
				array<int, 4> el = { spt.size() - 4, spt.size() - 3, spt.size() - 1, spt.size() - 2 };
				sele.push_back(el);
			}
			else if (bzmesh[e].bfc[k] == 2)
			{
				array<int, 4> el = { spt.size() - 4, spt.size() - 3, spt.size() - 1, spt.size() - 2 };
				sele.push_back(el);
			}
			else if (bzmesh[e].bfc[k] == 3)
			{
				array<int, 4> el = { spt.size() - 4, spt.size() - 2, spt.size() - 1, spt.size() - 3 };
				sele.push_back(el);
			}
			else if (bzmesh[e].bfc[k] == 4)
			{
				array<int, 4> el = { spt.size() - 4, spt.size() - 2, spt.size() - 1, spt.size() - 3 };
				sele.push_back(el);
			}
			else
			{
				array<int, 4> el = { spt.size() - 4, spt.size() - 3, spt.size() - 1, spt.size() - 2 };
				sele.push_back(el);
			}
		}
	}

	string fname = fn + "_surf.vtk";
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
		fout << "\nCELLS " << sele.size() << " " << 5 * sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "4 " << sele[i][0] << " " << sele[i][1] << " " << sele[i][2] << " " << sele[i][3] << '\n';
		}
		fout << "\nCELL_TYPES " << sele.size() << '\n';
		for (i = 0; i<sele.size(); i++)
		{
			fout << "9\n";
		}

		//fout<<"\nPOINT_DATA "<<sdisp.size()<<"\nVECTORS disp float\n";
		//for(i=0;i<sdisp.size();i++)
		//{
		//	fout << sdisp[i][0] << " " << sdisp[i][1] << " " << sdisp[i][2] << "\n";
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
}

void LeastSquare::VisualizeError(const vector<BezierElement3D>& bzmesh, const vector<double>& err, string fn)
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
			fout<<err[i]<<"\n";
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

void LeastSquare::Run(const vector<BezierElement3D>& bzmesh, string fn, vector<double>& err)
{
	GaussInfo(4);
	FindRange(bzmesh);
	SparseMatrix<double> GK(neq,neq);
	VectorXd GF(neq);
	GK.setZero();
	GF.setZero();
	cout<<"Building linear system...\n";
	BuildLinearSystem(bzmesh,GK,GF);
	OutputMatrix(GK,fn);
	//Solver(GK,GF);

	//VisualizeVTK(bzmesh,fn);
	////VisualizeVTK_1(bzmesh, fn);//with smoother boundary
	////VisualizeVTK_2(bzmesh, fn);//smooth boundary, remove elements
	//ErrorEstimate(bzmesh,err);//with exact solution
	////TotalErrorEstimate(bzmesh, err);//without exact solution
	//VisualizeError(bzmesh, err, fn);
}

void LeastSquare::Run_Fitting(const vector<BezierElement3D>& bzmesh, string fn, vector<double>& sol)
{
	cout << "Run fitting...\n";
	GaussInfo(5);
	//FindRange(bzmesh);
	SparseMatrix<double> GK(neq, neq);
	VectorXd GF(neq);
	//GK.setZero();
	GF.setZero();
	cout << "Building linear system...\n";
	//BuildLinearSystem(bzmesh, GK, GF);
	BuildLinearSystem_Fitting(bzmesh, GK, GF);
	//BuildLinearSystem_Fitting_C012(bzmesh, GK, GF);
	//OutputMatrix(GK,fn);
	Solver(GK, GF);

	sol.clear();
	sol = uh;

	cout << "Done fitting!\n";

	//VisualizeVTK(bzmesh,fn);
	////VisualizeVTK_1(bzmesh, fn);//with smoother boundary
	////VisualizeVTK_2(bzmesh, fn);//smooth boundary, remove elements
	//ErrorEstimate(bzmesh,err);//with exact solution
	////TotalErrorEstimate(bzmesh, err);//without exact solution
	//VisualizeError(bzmesh, err, fn);
}

void LeastSquare::OutputMatrix(const SparseMatrix<double>& GK, string fn)
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
		fout << neq*neq<<" "<< ratio << "\n";
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fn1 << "!\n";
	}

	//cout << "# of nonzeros: " << GK.nonZeros() << "\n";
}

//void Laplace::setProblem_ptestHO(const TruncatedTspline_3D& tts)
//{
//	npt=0;
//	unsigned int i,j,dof;
//	BCList.clear();
//	gh.clear();
//	for(i=0; i<tts.cp.size(); i++)
//	{
//		if(tts.cp[i].act==1)
//		{
//			if((tts.cp[i].coor[0]==0.)||
//				(tts.cp[i].coor[1]==0.)||
//				(tts.cp[i].coor[0]==1.)||
//				(tts.cp[i].coor[1]==1.))
//			{
//				BCList.push_back(tts.paid[i]);
//				//double tmp=ptestHO_sol_1(tts.cp[i].coor[0],tts.cp[i].coor[1]);
//				//double tmp=ptestHO_sol_2(tts.cp[i].coor[0],tts.cp[i].coor[1]);
//				double tmp=ptestHO_sol_3(tts.cp[i].coor[0],tts.cp[i].coor[1]);
//				gh.push_back(tmp);
//			}
//			else
//			{
//				gh.push_back(0.);
//			}
//			npt++;
//		}
//	}
//}
//
//double LDomainSolution(double x,double y)
//{
//	const double PI(3.14159265359);
//	double r=sqrt(x*x+y*y);
//	double phi=atan2(y,x);
//	if(phi<=0.) phi+=2.*PI;
//	double u=pow(r,0.6666666667)*sin(2.*phi/3.-PI/3.);
//	return u;
//}
//
//void LDomainSolutionDeriv(double x, double y, double& u, double& ux, double& uy)
//{
//	const double PI(3.14159265359);
//	double r=sqrt(x*x+y*y);
//	double phi=atan2(y,x);
//	if(phi<=0.) phi+=2.*PI;
//	u=pow(r,0.6666666667)*sin(2.*phi/3.-PI/3.);
//	ux=0.6666666667*pow(r,-0.3333333333)*(sin(2.*phi/3.-PI/3.)*cos(phi)-cos(2.*phi/3.-PI/3.)*sin(phi));
//	uy=0.6666666667*pow(r,-0.3333333333)*(sin(2.*phi/3.-PI/3.)*sin(phi)+cos(2.*phi/3.-PI/3.)*cos(phi));
//}
//
//double ptestHO_sol_1(double x, double y)
//{
//	return (x+y);
//}
//
//double ptestHO_sol_2(double x, double y)
//{
//	return (x*x+y*y);
//	//return (x*x);
//}
//
//double ptestHO_sol_3(double x, double y)
//{
//	//return (x*x*x+y*y*y);
//	double u=x*(1.-x)*y*(1.-y);
//	return u;
//}

void LeastSquare::FindRange(const vector<BezierElement3D>& bzmesh)
{
	double tmp[3][2] = { { 1.e5, -1.e5 }, { 1.e5, -1.e5 }, { 1.e5, -1.e5 }};
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

double LeastSquare::f_source_ls(const array<double, 3>& x)
{
	//double tmp = f_source_1(x);
	//double tmp = f_source_2(x);
	//double tmp = f_source_3(x);
	//double tmp = f_source_4(x);
	//double tmp = f_source_5(x);
	//double tmp = f_source_6(x);
	double tmp = f_source_7_ls(x);
	//double tmp = f_source_8(x);
	return tmp;
}

//double LeastSquare::f_source_1(const array<double, 3>& x)
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
//double LeastSquare::f_source_2(const array<double, 3>& x)//analysis domain is a unit cube
//{
//	double tmp1 = -1600.*exp(-654. + 800.*x[0] * (1. - x[0]) + 800.*x[1] * (1. - x[1]) + 800.*x[2] * (1. - x[2]))/3.;
//	double tmp2 = exp((7. - 20.*x[0])*(7. - 20.*x[0]) + (7. - 20.*x[1])*(7. - 20.*x[1]) + (7. - 20.*x[2])*(7. - 20.*x[2]));
//	double tmp3 = exp((13. - 20.*x[0])*(13. - 20.*x[0]) + (13. - 20.*x[1])*(13. - 20.*x[1]) + (13. - 20.*x[2])*(13. - 20.*x[2]));
//	double tmp4 = 1011. - 1040.*x[0] + 800.*x[0] * x[0] - 1040.*x[1] + 800.*x[1] * x[1] - 1040.*x[2] + 800.*x[2] * x[2];
//	double tmp5 = 291. - 560.*x[0] + 800.*x[0] * x[0] - 560.*x[1] + 800.*x[1] * x[1] - 560.*x[2] + 800.*x[2] * x[2];
//	return tmp1*(tmp2*tmp4+tmp3*tmp5);
//}
//
//double LeastSquare::f_source_3(const array<double, 3>& x)//analysis domain is a unit cube
//{
//	double tmp1 = -1600.*exp(-300. + 400.*x[0] * (1. - x[0]) + 400.*x[1] * (1. - x[1]) + 400.*x[2] * (1. - x[2])) / 3.;
//	double tmp2 = 597. - 800.*x[0] * (1. - x[0]) - 800.*x[1] * (1. - x[1]) - 800.*x[2] * (1. - x[2]);
//	return tmp1*tmp2;
//}
//
//double LeastSquare::f_source_4(const array<double, 3>& x)//analysis domain is a unit cube
//{
//	double tmp0 = 1. + 20.4124*x[0] + 20.4124*x[1] - 40.8248*x[2];
//	double tmp1 = 1. / cosh(tmp0);
//	return 5000.*tmp1*tmp1*tanh(tmp0);
//}
//
//double LeastSquare::f_source_5(const array<double, 3>& x)//analysis domain is arbitrary
//{
//	double tmp = -10.*(x[0] + x[1] + x[2]);
//	return tmp;
//}
//
//double LeastSquare::f_source_6(const array<double, 3>& x)//analysis domain is arbitrary
//{
//	double tmp = 2.*(x[0] * (1. - x[0])*x[1] * (1. - x[1]) + x[0] * (1. - x[0])*x[2] * (1. - x[2]) + x[2] * (1. - x[2])*x[1] * (1. - x[1]));
//	return tmp;
//}
//
double LeastSquare::f_source_7_ls(const array<double, 3>& x)
{
	//double tmp0 = acoef*(nmpl[0] * (x[0] - range[0][0]) + nmpl[1] * (x[1] - range[1][0]) + nmpl[2] * (x[2] - range[2][0]));
	double tmp0 = acoef*(nmpl[0] * (x[0] - xorg[0]) + nmpl[1] * (x[1] - xorg[1]) + nmpl[2] * (x[2] - xorg[2]));
	double tmp1 = 1. / cosh(tmp0);
	return bcoef*tmp1*tmp1*tanh(tmp0);
}

//double LeastSquare::f_source_8(const array<double, 3>& x)
//{
//	double x1[3] = { (x[0] - range[0][0]) / dmlen[0], (x[1] - range[1][0]) / dmlen[1], (x[2] - range[2][0]) / dmlen[2] };
//	double len2[3] = { dmlen[0] * dmlen[0], dmlen[1] * dmlen[1], dmlen[2] * dmlen[2] };
//	double tmp0 = -acoef / exp(acoef*((x1[0] - nmpl[0])*(x1[0] - nmpl[0]) + (x1[1] - nmpl[1])*(x1[1] - nmpl[1]) + (x1[2] - nmpl[2])*(x1[2] - nmpl[2])));
//	double tmp1 = (-2. + acoef)*(1. / len2[0] + 1. / len2[1] + 1. / len2[2]);
//	double tmp2 = 4.*acoef*(x1[0] * (x1[0] - 1) / len2[0] + x1[1] * (x1[1] - 1) / len2[1] + x1[2] * (x1[2] - 1) / len2[2]);
//	return tmp0*(tmp1+tmp2);
//}
//
double LeastSquare::exact_sol_ls(const array<double, 3>& x)
{
	//double tmp = exact_sol_1(x);
	//double tmp = exact_sol_2(x);
	//double tmp = exact_sol_3(x);
	//double tmp = exact_sol_4(x);
	//double tmp = exact_sol_5(x);
	//double tmp = exact_sol_6(x);
	double tmp = exact_sol_7_ls(x);
	//double tmp = exact_sol_8(x);
	return tmp;
}

//double LeastSquare::exact_sol_1(const array<double, 3>& x)
//{
//	double nm = (range[0][1] - range[0][0])*(range[0][1] - range[0][0])*
//		(range[1][1] - range[1][0])*(range[1][1] - range[1][0])*
//		(range[2][1] - range[2][0])*(range[2][1] - range[2][0]);
//	return (x[0] - range[0][0])*(range[0][1] - x[0])*(x[1] - range[1][0])*(range[1][1] - x[1])*(x[2] - range[2][0])*(range[2][1] - x[2])/nm;
//}
//
//double LeastSquare::exact_sol_2(const array<double, 3>& x)
//{
//	double tmp1 = 1. / exp((20.*x[0] - 7.)*(20.*x[0] - 7.) + (20.*x[1] - 7.)*(20.*x[1] - 7.) + (20.*x[2] - 7.)*(20.*x[2] - 7.));
//	double tmp2 = 1. / exp((20.*x[0] - 13.)*(20.*x[0] - 13.) + (20.*x[1] - 13.)*(20.*x[1] - 13.) + (20.*x[2] - 13.)*(20.*x[2] - 13.));
//	return 2.*(tmp1+tmp2)/3.;
//}
//
//double LeastSquare::exact_sol_3(const array<double, 3>& x)
//{
//	return 2. / (3.*exp((20.*x[0] - 10.)*(20.*x[0] - 10.) + (20.*x[1] - 10.)*(20.*x[1] - 10.) + (20.*x[2] - 10.)*(20.*x[2] - 10.)));
//}
//
//double LeastSquare::exact_sol_4(const array<double, 3>& x)
//{
//	return tanh(1. - 50.*(-0.408248*x[0]-0.408248*x[1] +0.816497*x[2]));
//}
//
//double LeastSquare::exact_sol_5(const array<double, 3>& x)
//{
//	double tmp = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2])*(x[0] + x[1] + x[2]) + x[0] * x[1] * x[2];
//	return tmp;
//}
//
//double LeastSquare::exact_sol_6(const array<double, 3>& x)
//{
//	double tmp = x[0] * (1. - x[0])*x[1] * (1. - x[1])*x[2] * (1. - x[2]);;
//	return tmp;
//}
//
double LeastSquare::exact_sol_7_ls(const array<double, 3>& x)
{
	//double x1[3] = { (x[0] - range[0][0]), (x[1] - range[1][0]), (x[2] - range[2][0]) };
	double x1[3] = { (x[0] - xorg[0]), (x[1] - xorg[1]), (x[2] - xorg[2]) };
	return tanh(acoef*(nmpl[0] * x1[0] + nmpl[1] * x1[1] + nmpl[2] * x1[2]));
}

//double LeastSquare::exact_sol_8(const array<double, 3>& x)
//{
//	double x1[3] = { (x[0] - range[0][0]) / dmlen[0], (x[1] - range[1][0]) / dmlen[1], (x[2] - range[2][0]) / dmlen[2] };
//	return 1. / exp(acoef*((x1[0] - nmpl[0])*(x1[0] - nmpl[0]) + (x1[1] - nmpl[1])*(x1[1] - nmpl[1]) + (x1[2] - nmpl[2])*(x1[2] - nmpl[2])));
//}

//double f_source_2(double x, double y)
//{
//	return -4.;
//	//return -2.;
//}
//
//double f_source_3(double x, double y)
//{
//	//return -6.*(x+y);
//	double f=x+y-x*x-y*y;
//	return f;
//}

void LeastSquare::GetRemoveRegion(double xy[3][2])
{
	rmv[0][0] = xy[0][0]; rmv[0][1] = xy[0][1];
	rmv[1][0] = xy[1][0]; rmv[1][1] = xy[1][1];
	rmv[2][0] = xy[2][0]; rmv[2][1] = xy[2][1];
}

void LeastSquare::GetEqParameter(double xy[3][2], double nrm[3], double a)
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

















void LeastSquare::SetProblem_SurfaceFitting(const vector<int>& IDBC_in)
{
	npt = IDBC_in.size();
	IDBC = IDBC_in;
	neq = 0;
	//int pid(0);
	for (unsigned int i = 0; i < IDBC.size(); i++)
	{
		if (IDBC[i] != -1) neq++;
	}
}

void LeastSquare::Run_SurfaceFitting(const vector<BezierElement3D>& bzmesh, string fn, vector<array<double,3>>& sol)
{
	cout << "Run fitting...\n";
	SparseMatrix<double> GK(neq, neq);
	vector<VectorXd> GF(3);
	for (uint i = 0; i < GF.size(); i++)
	{
		GF[i] = VectorXd::Zero(neq);
	}
	cout << "Building linear system...\n";
	BuildLinearSystem_SurfaceFitting(bzmesh, GK, GF);
	//OutputMatrix(GK,fn);
	Solver_SurfaceFitting(GK, GF);

	sol.clear();
	sol = uh3;

	cout << "Done fitting!\n";

	//VisualizeVTK(bzmesh,fn);
	////VisualizeVTK_1(bzmesh, fn);//with smoother boundary
	////VisualizeVTK_2(bzmesh, fn);//smooth boundary, remove elements
	//ErrorEstimate(bzmesh,err);//with exact solution
	////TotalErrorEstimate(bzmesh, err);//without exact solution
	//VisualizeError(bzmesh, err, fn);
}

void LeastSquare::BuildLinearSystem_SurfaceFitting(const vector<BezierElement3D>& bzmesh, SparseMatrix<double>& GK, vector<VectorXd>& GF)
{
	int nresv(5 * 64);
	if (nresv < neq)
	{
		GK.reserve(VectorXi::Constant(neq, nresv));
	}
	else
	{
		GK.reserve(VectorXi::Constant(neq, neq));
	}
	int nsp(2);
	double rng[2] = { 0., 1. };
	vector<double> spt(nsp);
	for (int i = 0; i < nsp; i++)
	{
		spt[i] = rng[0] + double(i) * (rng[1] - rng[0]) / double(nsp - 1);
	}
	//int ploc[2][2] = { {0,3},{1,2} };

	cout << "npt: " << npt << "\n";
	cout << "neq: " << neq << "\n";
	cout << "nel: " << bzmesh.size() << '\n';

	//#pragma omp parallel for
	for (int e = 0; e<bzmesh.size(); e++)
	{
		if (e != 0 && e % 1000 == 0)
		{
			cout << e << " ";
		}
		if (bzmesh[e].type == 1)
		{
			for (int k = 0; k < 6; k++)
			{
				if (bzmesh[e].bc[k] == 1)
				{
					double Nx[64];
					double dNdx[64][3];
					vector<double> Nx1(bzmesh[e].IEN.size());
					vector<vector<double>> EK(bzmesh[e].IEN.size(), vector<double>(bzmesh[e].IEN.size(), 0.));
					vector<vector<double>> EF(3, vector<double>(bzmesh[e].IEN.size(), 0.));
					double detJ(1.), Fb[3] = { 0.,0.,0. };
					array<double, 3> pt;
					for (int i = 0; i<spt.size(); i++)
					{
						for (int j = 0; j<spt.size(); j++)
						{
							if (k == 0)
							{
								BasisFunction(spt[i], spt[j], 0., bzmesh[e].pts, Nx, dNdx, detJ);
								Para2Phys_SurfaceFitting(spt[i], spt[j], 0., bzmesh[e].cnpt, pt);
							}
							else if (k == 1)
							{
								BasisFunction(spt[i], 0., spt[j], bzmesh[e].pts, Nx, dNdx, detJ);
								Para2Phys_SurfaceFitting(spt[i], 0., spt[j], bzmesh[e].cnpt, pt);
							}
							else if (k == 2)
							{
								BasisFunction(1., spt[i], spt[j], bzmesh[e].pts, Nx, dNdx, detJ);
								Para2Phys_SurfaceFitting(1., spt[i], spt[j], bzmesh[e].cnpt, pt);
							}
							else if (k == 3)
							{
								BasisFunction(spt[i], 1., spt[j], bzmesh[e].pts, Nx, dNdx, detJ);
								Para2Phys_SurfaceFitting(spt[i], 1., spt[j], bzmesh[e].cnpt, pt);
							}
							else if (k == 4)
							{
								BasisFunction(0., spt[i], spt[j], bzmesh[e].pts, Nx, dNdx, detJ);
								Para2Phys_SurfaceFitting(0., spt[i], spt[j], bzmesh[e].cnpt, pt);
							}
							else
							{
								BasisFunction(spt[i], spt[j], 1., bzmesh[e].pts, Nx, dNdx, detJ);
								Para2Phys_SurfaceFitting(spt[i], spt[j], 1., bzmesh[e].cnpt, pt);
							}
							BasisFunctionMap(Nx, bzmesh[e].cmat, Nx1);
							ElementMatrix(Nx1, EK);
							ElementForce_SurfaceFitting(Nx1, pt.data(), EF);
						}
					}
					Assembly_SurfaceFitting(EK, EF, bzmesh[e].IEN, GK, GF);
				}
			}
		}
	}
}

void LeastSquare::ElementForce_SurfaceFitting(const vector<double>& Nx, double Fb[3], vector<vector<double>>& EF)
{
	for (int i = 0; i < Nx.size(); i++)
	{
		EF[0][i] += Nx[i] * Fb[0];
		EF[1][i] += Nx[i] * Fb[1];
		EF[2][i] += Nx[i] * Fb[2];
	}
}

void LeastSquare::Para2Phys_SurfaceFitting(double u, double v, double w, const vector<array<double, 3>>& cnpt, array<double, 3>& pt)
{
	pt[0] = 0.; pt[1] = 0.; pt[2] = 0.;
	double Nu[2] = { 1. - u,u };
	double Nv[2] = { 1. - v,v };
	double Nw[2] = { 1. - w,w };
	int ploc[8] = { 0,1,3,2,4,5,7,6 };
	int loc(0);
	double tmp;
	for (int i = 0; i < 2; i++)//w
	{
		for (int j = 0; j < 2; j++)//v
		{
			for (int k = 0; k < 2; k++)//u
			{
				tmp = Nu[k] * Nv[j] * Nw[i];
				pt[0] += cnpt[ploc[loc]][0] * tmp;
				pt[1] += cnpt[ploc[loc]][1] * tmp;
				pt[2] += cnpt[ploc[loc]][2] * tmp;
				loc++;
			}
		}
	}
}

void LeastSquare::Assembly_SurfaceFitting(const vector<vector<double>>& EK, const vector<vector<double>>& EF, const vector<int>& IEN, 
	SparseMatrix<double>& GK, vector<VectorXd>& GF)
{
	int i, j, A, B;
	for (i = 0; i<IEN.size(); i++)
	{
		if (IEN[i] != -1)
		{
			for (j = 0; j<IEN.size(); j++)
			{
				if (IEN[j] != -1)
				{
					if (IDBC[IEN[i]] != -1 && IDBC[IEN[j]] != -1)
					{
						GK.coeffRef(IDBC[IEN[i]], IDBC[IEN[j]]) += EK[i][j];
					}
				}
			}
		}
	}
	for (i = 0; i < IEN.size(); i++)
	{
		if (IEN[i] != -1)
		{
			if (IDBC[IEN[i]] != -1)
			{
				GF[0][IDBC[IEN[i]]] += EF[0][i];
				GF[1][IDBC[IEN[i]]] += EF[1][i];
				GF[2][IDBC[IEN[i]]] += EF[2][i];
			}
		}
	}
}

void LeastSquare::Solver_SurfaceFitting(SparseMatrix<double>& GK, vector<VectorXd>& GF)
{
	SimplicialLDLT<SparseMatrix<double>> solver;
	solver.compute(GK);
	vector<VectorXd> sol(3);
	sol[0] = solver.solve(GF[0]);
	sol[1] = solver.solve(GF[1]);
	sol[2] = solver.solve(GF[2]);
	uh3.resize(npt);
	for (int i = 0; i<npt; i++)
	{
		if (IDBC[i] != -1)
		{
			uh3[i][0] = sol[0][IDBC[i]];
			uh3[i][1] = sol[1][IDBC[i]];
			uh3[i][2] = sol[2][IDBC[i]];
		}
		else
		{
			uh3[i][0] = 0.;
			uh3[i][1] = 0.;
			uh3[i][2] = 0.;
		}
	}
}


