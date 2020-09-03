#include "BasicDataStructure.h"
#include <algorithm>
#include "KnotInsertion.h"

Edge::Edge()
{
	act=1;
	pt[0]=-1; pt[1]=-1;
	//face[0]=-1; face[1]=-1;
	lev=0;
	len=1.;
	//pned[0]=-1; pned[1]=-1;
	//pnfc[0]=-1; pnfc[1]=-1;
	pn[0][0]=3; pn[1][0]=3;//0 for edge, 1 for face (1st edge diretion), 2 for face (2nd edge direction), 3 for end
	pn[0][1]=-1; pn[1][1]=-1;//edge or face ID
	prt=-1;
	chd[0]=-1; chd[0]=-1;
	midpt=-1;
}

bool Edge::operator==(const Edge& ed)
{
	return ((pt[0]==ed.pt[0] && pt[1]==ed.pt[1])||(pt[0]==ed.pt[1] && pt[1]==ed.pt[0]));
}

/////////////////////////////////////////////////////


Edge3D::Edge3D()
{
	act=1;
	id_act = -1;
	pt[0]=-1; pt[1]=-1;
	lev=0;
	len=1.;
	type = 0;
	sharp = 0;
	pn[0][0]=4; pn[1][0]=4;//0 for edge, 1 for face (1st edge diretion), 2 for face (2nd edge direction), 3 for hex, 4 for end
	pn[0][1]=-1; pn[1][1]=-1;//edge, face or hex ID
	prt=-1;
	chd[0]=-1; chd[0]=-1;
	midpt=-1;

	mpflag = 0;

	c0flag = 0;
	c0flag_b = 0;

	ped_id = -1;
}

bool Edge3D::operator==(const Edge3D& ed)
{
	return ((pt[0]==ed.pt[0] && pt[1]==ed.pt[1])||(pt[0]==ed.pt[1] && pt[1]==ed.pt[0]));
}

/////////////////////////////////////////////////////

Face3D::Face3D()
{
	act=1;
	id_act = -1;
	lev=0.;
	type=0;
	ctpt=-1;
	prt=-1;
	cnct[0]=0; cnct[1]=0; cnct[2]=0; cnct[3]=0;
	edge[0]=0; edge[1]=0; edge[2]=0; edge[3]=0;
	pn[0][0]=4; pn[1][0]=4; pn[2][0]=4; pn[3][0]=4;
	pn[0][1]=-1; pn[1][1]=-1; pn[2][1]=-1; pn[3][1]=-1;

	mpflag = 0;

	c0flag = 0;

	pfc_id = -1;
}

bool Face3D::operator==(const Face3D& fc)
{
	array<int,4> tmp1={cnct[0],cnct[1],cnct[2],cnct[3]};
	array<int,4> tmp2={fc.cnct[0],fc.cnct[1],fc.cnct[2],fc.cnct[3]};
	sort(tmp1.begin(),tmp1.end());
	sort(tmp2.begin(),tmp2.end());
	return (tmp1==tmp2);
}

////////////////////////////////////////////////////

Vertex2D::Vertex2D()
{
	coor[0]=0.;	coor[1]=0.;	coor[2]=0.;
	coortmp[0]=0.;	coortmp[1]=0.;	coortmp[2]=0.;
	type=0;
	trun=0;
	aff=0;
	update=0;
	kitvU[0]=1.; kitvU[1]=1.; kitvU[2]=1.; kitvU[3]=1.;
	kitvV[0]=1.; kitvV[1]=1.; kitvV[2]=1.; kitvV[3]=1.;
	kitvUtmp[0]=1.; kitvUtmp[1]=1.; kitvUtmp[2]=1.; kitvUtmp[3]=1.;
	kitvVtmp[0]=1.; kitvVtmp[1]=1.; kitvVtmp[2]=1.; kitvVtmp[3]=1.;
}

//Vertex Vertex::operator+(const Vertex& v)
//{
//	Vertex a;
//	a.coor[0]=coor[0]+v.coor[0];
//	a.coor[1]=coor[1]+v.coor[1];
//	a.coor[2]=coor[2]+v.coor[2];
//	return a;
//}
//
//Vertex Vertex::operator-(const Vertex& v)
//{
//	Vertex a;
//	a.coor[0]=coor[0]-v.coor[0];
//	a.coor[1]=coor[1]-v.coor[1];
//	a.coor[2]=coor[2]-v.coor[2];
//	return a;
//}
//
//Vertex Vertex::operator*(double a)
//{
//	Vertex v;
//	v.coor[0]=coor[0]*a;
//	v.coor[1]=coor[1]*a;
//	v.coor[2]=coor[2]*a;
//	return v;
//}
//
//Vertex Vertex::operator/(double a)
//{
//	Vertex v;
//	v.coor[0]=coor[0]/a;
//	v.coor[1]=coor[1]/a;
//	v.coor[2]=coor[2]/a;
//	return v;
//}

//bool Vertex2D::operator==(const Vertex2D& v)
//{
//	double tol(1.e-8);
//	double dif=sqrt((coor[0]-v.coor[0])*(coor[0]-v.coor[0])+(coor[1]-v.coor[1])*(coor[1]-v.coor[1])+(coor[2]-v.coor[2])*(coor[2]-v.coor[2]));
//	return (dif<tol);
//}

Vertex3D::Vertex3D()
{
	coor[0]=0.;	coor[1]=0.;	coor[2]=0.;
	coortmp[0]=0.;	coortmp[1]=0.;	coortmp[2]=0.;
	wght = 1.;
	type=0;
	trun=0;
	aff = 0;
	update = 0;
	truntmp = 0;
	lev = 0;
	act = 0;
	sharp = 0;
	bcxp = 0;
	kitvU[0]=1.; kitvU[1]=1.; kitvU[2]=1.; kitvU[3]=1.;
	kitvV[0]=1.; kitvV[1]=1.; kitvV[2]=1.; kitvV[3]=1.;
	kitvW[0]=1.; kitvW[1]=1.; kitvW[2]=1.; kitvW[3]=1.;

	mpflag = 0;
	mpid = -1;

	bzid = -1;
	c0flag = 0;
	c0flag_b = 0;

	smth = 0;
	pvt_id = -1;
}

bool Vertex3D::operator==(const Vertex3D& v)
{
	double tol(1.e-6);
	double dif=sqrt((coor[0]-v.coor[0])*(coor[0]-v.coor[0])+(coor[1]-v.coor[1])*(coor[1]-v.coor[1])+(coor[2]-v.coor[2])*(coor[2]-v.coor[2]));
	return (dif<tol);
}

////////////////////////////////////////////////////////////

Element2D::Element2D()
{
	cnct[0]=0; cnct[1]=0; cnct[2]=0; cnct[3]=0;
	edge[0]=0; edge[1]=0; edge[2]=0; edge[3]=0;
	act=1;
	lev=0;
	type=0;
	prt=-1;
}

Element3D::Element3D()
{
	for(int i=0; i<8; i++) cnct[i]=0;
	for(int i=0; i<12; i++) edge[i]=-1;
	for(int i=0; i<6; i++) face[i]=-1;
	act=1;
	id_act = -1;
	lev=0;
	type=0;
	prt=-1;
	ref_flag = 0;
	trun = 0;
	ghost = 0;
	bc_lev = 0;
	bc = 0;
	smth = 0;
	for (int i = 0; i < 6; i++)
	{
		nbrot[i].Zero();
		nbrot[i](0, 0) = 1.;
		nbrot[i](1, 1) = 1.;
		nbrot[i](2, 2) = 1.;
	}
	dm[0][0] = 0.; dm[0][1] = 1.;
	dm[1][0] = 0.; dm[1][1] = 1.;
	dm[2][0] = 0.; dm[2][1] = 1.;

	bzpt.resize(64);
	for (int i = 0; i < 64; i++)
	{
		for (int j = 0; j < 3; j++) bzpt[i][j] = 0.;
	}

	kvlen[0] = 1.; kvlen[1] = 1.; kvlen[2] = 1.;
	for (int i = 0; i < 6; i++)
	{
		fcnb[i][0] = -1; fcnb[i][1] = -1;
	}
	for (int i = 0; i < 12; i++)
	{
		ednb[i][0] = -1; ednb[i][1] = -1;
	}
	for (int i = 0; i < 8; i++)
	{
		vtnb[i][0] = -1; vtnb[i][1] = -1;
	}
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 8; j++) pkv[i][j] = 1.;
	}

	bzflag = 0;

	pbd_id = -1;
	jacobFlag = 0;
}

void Element3D::Clear()
{
	unsigned int i;
	vector<int>().swap(IEN);
	for (i = 0; i < patch_ku.size(); i++) array<double, 5>().swap(patch_ku[i]);
	vector<array<double, 5>>().swap(patch_ku);
	for (i = 0; i < patch_kv.size(); i++) array<double, 5>().swap(patch_kv[i]);
	vector<array<double, 5>>().swap(patch_kv);
	for (i = 0; i < patch_kw.size(); i++) array<double, 5>().swap(patch_kw[i]);
	vector<array<double, 5>>().swap(patch_kw);

	for (i = 0; i < bemat.size(); i++)
	{
		vector<double>().swap(bemat[i]);
	}
	vector<vector<double>>().swap(bemat);
}

void Element3D::Initialize()
{
	for (int i = 0; i<8; i++) cnct[i] = 0;
	for (int i = 0; i<12; i++) edge[i] = -1;
	for (int i = 0; i<6; i++) face[i] = -1;
	act = 1;
	id_act = -1;
	lev = 0;
	type = 0;
	prt = -1;
	ref_flag = 0;
	trun = 0;
	ghost = 0;
	bc_lev = 0;
	for (int i = 0; i < 6; i++)
	{
		nbrot[i].Zero();
		nbrot[i](0, 0) = 1.;
		nbrot[i](1, 1) = 1.;
		nbrot[i](2, 2) = 1.;
	}
	dm[0][0] = 0.; dm[0][1] = 1.;
	dm[1][0] = 0.; dm[1][1] = 1.;
	dm[2][0] = 0.; dm[2][1] = 1.;

	bzpt.resize(64);
	for (int i = 0; i < 64; i++)
	{
		for (int j = 0; j < 3; j++) bzpt[i][j] = 0.;
	}

	kvlen[0] = 1.; kvlen[1] = 1.; kvlen[2] = 1.;
	for (int i = 0; i < 6; i++)
	{
		fcnb[i][0] = -1; fcnb[i][1] = -1;
	}
	for (int i = 0; i < 12; i++)
	{
		ednb[i][0] = -1; ednb[i][1] = -1;
	}
	for (int i = 0; i < 8; i++)
	{
		vtnb[i][0] = -1; vtnb[i][1] = -1;
	}
}

//bool PointIPP::operator==(const PointIPP& pt)
//{
//	if(index[0]==pt.index[0] && index[1]==pt.index[1] && pm[0]==pt.pm[0] && pm[1]==pt.pm[1])
//	{
//		return true;
//	}
//	else
//	{
//		return false;
//	}
//}

////////////////////////////////////////////////////////////////

BezierElement2D::BezierElement2D(int p)
{
	degree=p;
	order=p+1;
	nbf=order*order;
	pts.resize(nbf);
	for(int i=0; i<nbf; i++)
	{
		pts[i][0]=0.; pts[i][1]=0.; pts[i][2]=0.;
	}
	//for(int i=0;i<25;i++)
	//{
	//	pts4[i][0]=0.; pts4[i][1]=0.; pts4[i][2]=0.;
	//}
}

void BezierElement2D::BezierPolyn(double u, vector<double>& Nu, vector<double>& dNdu) const
{
	if(degree==3)
	{
		double Nu0[4]={(1.-u)*(1.-u)*(1.-u),3.*(1.-u)*(1.-u)*u,3.*(1.-u)*u*u,u*u*u};
		double dNdu0[4]={-3.*(1.-u)*(1.-u),3.-12.*u+9.*u*u,3.*(2.-3.*u)*u,3.*u*u};
		Nu.resize(order);
		dNdu.resize(order);
		for(int i=0; i<order; i++)
		{
			Nu[i]=Nu0[i];
			dNdu[i]=dNdu0[i];
		}
	}
	else if(degree==4)
	{
		double Nu0[5]={(1.-u)*(1.-u)*(1.-u)*(1.-u),4.*(1.-u)*(1.-u)*(1.-u)*u,6.*(1.-u)*(1.-u)*u*u,4.*(1.-u)*u*u*u,u*u*u*u};
		double dNdu0[5]={-4.*(1.-u)*(1.-u)*(1.-u),4.*(1.-u)*(1.-u)*(1.-4.*u),12.*u*(1.-3.*u+2.*u*u),4.*(3.-4.*u)*u*u,4.*u*u*u};
		Nu.resize(order);
		dNdu.resize(order);
		for(int i=0; i<order; i++)
		{
			Nu[i]=Nu0[i];
			dNdu[i]=dNdu0[i];
		}
	}
}

void BezierElement2D::Basis(double u, double v, vector<double>& Nt, vector<array<double,2>>& dNdt) const
{
	vector<double> Nu, Nv, dNdu, dNdv;
	BezierPolyn(u,Nu,dNdu);
	BezierPolyn(v,Nv,dNdv);
	Nt.resize(nbf);
	dNdt.resize(nbf);
	int i,j,loc(0);
	for(i=0; i<order; i++)
	{
		for(j=0; j<order; j++)
		{
			Nt[loc]=Nu[j]*Nv[i];
			dNdt[loc][0]=dNdu[j]*Nv[i];
			dNdt[loc][1]=Nu[j]*dNdv[i];
			loc++;
		}
	}
}

//void BezierElement::Basis4(double u, double v, double Nt[25], double dNdt[25][2]) const
//{
//	double Nu[5]={(1.-u)*(1.-u)*(1.-u)*(1.-u),4.*(1.-u)*(1.-u)*(1.-u)*u,6.*(1.-u)*(1.-u)*u*u,4.*(1.-u)*u*u*u,u*u*u*u};
//	double Nv[5]={(1.-v)*(1.-v)*(1.-v)*(1.-v),4.*(1.-v)*(1.-v)*(1.-v)*v,6.*(1.-v)*(1.-v)*v*v,4.*(1.-v)*v*v*v,v*v*v*v};
//	double dNdu[5]={-4.*(1.-u)*(1.-u)*(1.-u),4.*(1.-u)*(1.-u)*(1.-4.*u),12.*u*(1.-3.*u+2.*u*u),4.*(3.-4.*u)*u*u,4.*u*u*u};
//	double dNdv[5]={-4.*(1.-v)*(1.-v)*(1.-v),4.*(1.-v)*(1.-v)*(1.-4.*v),12.*v*(1.-3.*v+2.*v*v),4.*(3.-4.*v)*v*v,4.*v*v*v};
//	int i,j,loc(0);
//	for(i=0;i<5;i++)
//	{
//		for(j=0;j<5;j++)
//		{
//			Nt[loc]=Nu[j]*Nv[i];
//			dNdt[loc][0]=dNdu[j]*Nv[i];
//			dNdt[loc][1]=Nu[j]*dNdv[i];
//			loc++;
//		}
//	}
//}

void BezierElement2D::Para2Phys(double u, double v, double pt[3])
{
	vector<double> Nt;
	vector<array<double,2>> dNdt;
	Basis(u,v,Nt,dNdt);
	pt[0]=0.; pt[1]=0.; pt[2]=0.;
	for(unsigned int i=0;i<pts.size();i++)
	{
		pt[0]+=pts[i][0]*Nt[i];
		pt[1]+=pts[i][1]*Nt[i];
		pt[2]+=pts[i][2]*Nt[i];
	}
}

//void BezierElement::Para2Phys4(double u, double v, double pt[3])
//{
//	double Nx[25];
//	double dNdt[25][2];
//	Basis4(u,v,Nx,dNdt);
//	pt[0]=0.; pt[1]=0.; pt[2]=0.;
//	for(int i=0;i<25;i++)
//	{
//		pt[0]+=pts4[i][0]*Nx[i];
//		pt[1]+=pts4[i][1]*Nx[i];
//		pt[2]+=pts4[i][2]*Nx[i];
//	}
//}

void BezierElement2D::SurfPointNormal(double u, double v, array<double,3>& pt, array<double,3>& nm) const
{
	vector<double> Nt;
	vector<array<double,2>> dNdt;
	Basis(u,v,Nt,dNdt);
	pt[0]=0.; pt[1]=0.; pt[2]=0.;
	nm[0]=0.; nm[1]=0.; nm[2]=0.;
	double nmtmp[2][3]={{0.,0.,0.},{0.,0.,0.}};
	for(int i=0; i<nbf; i++)
	{
		pt[0]+=pts[i][0]*Nt[i];
		pt[1]+=pts[i][1]*Nt[i];
		pt[2]+=pts[i][2]*Nt[i];
		nmtmp[0][0]+=pts[i][0]*dNdt[i][0];
		nmtmp[0][1]+=pts[i][1]*dNdt[i][0];
		nmtmp[0][2]+=pts[i][2]*dNdt[i][0];
		nmtmp[1][0]+=pts[i][0]*dNdt[i][1];
		nmtmp[1][1]+=pts[i][1]*dNdt[i][1];
		nmtmp[1][2]+=pts[i][2]*dNdt[i][1];
	}
	nm[0]=nmtmp[0][1]*nmtmp[1][2]-nmtmp[0][2]*nmtmp[1][1];
	nm[1]=nmtmp[0][2]*nmtmp[1][0]-nmtmp[0][0]*nmtmp[1][2];
	nm[2]=nmtmp[0][0]*nmtmp[1][1]-nmtmp[0][1]*nmtmp[1][0];
	double len=sqrt(nm[0]*nm[0]+nm[1]*nm[1]+nm[2]*nm[2]);
	nm[0]/=len; nm[1]/=len; nm[2]/=len;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

BezierElement3D::BezierElement3D(int p)
{
	degree=p;
	order=p+1;
	nbf=order*order*order;
	pts.resize(nbf);
	for(int i=0; i<nbf; i++)
	{
		pts[i][0]=0.; pts[i][1]=0.; pts[i][2]=0.;
	}
	prt[0] = 0; prt[1] = 0;
	trun = 0;
	type = 0;

	rkpm = 0;
	for (int i = 0; i < 6; i++) bc[i] = 0;

	for (int i = 0; i < 6; i++) bcval[i] = 0.;

	bzcouple = 0;
	bzflag = 0;
	bcflag = 0;
}

void BezierElement3D::BezierPolyn(double u, vector<double>& Nu, vector<double>& dNdu) const
{
	if(degree==3)
	{
		double Nu0[4]={(1.-u)*(1.-u)*(1.-u),3.*(1.-u)*(1.-u)*u,3.*(1.-u)*u*u,u*u*u};
		double dNdu0[4]={-3.*(1.-u)*(1.-u),3.-12.*u+9.*u*u,3.*(2.-3.*u)*u,3.*u*u};
		Nu.resize(order);
		dNdu.resize(order);
		for(int i=0; i<order; i++)
		{
			Nu[i]=Nu0[i];
			dNdu[i]=dNdu0[i];
		}
	}
	else if(degree==4)
	{
		double Nu0[5]={(1.-u)*(1.-u)*(1.-u)*(1.-u),4.*(1.-u)*(1.-u)*(1.-u)*u,6.*(1.-u)*(1.-u)*u*u,4.*(1.-u)*u*u*u,u*u*u*u};
		double dNdu0[5]={-4.*(1.-u)*(1.-u)*(1.-u),4.*(1.-u)*(1.-u)*(1.-4.*u),12.*u*(1.-3.*u+2.*u*u),4.*(3.-4.*u)*u*u,4.*u*u*u};
		Nu.resize(order);
		dNdu.resize(order);
		for(int i=0; i<order; i++)
		{
			Nu[i]=Nu0[i];
			dNdu[i]=dNdu0[i];
		}
	}
}

void BezierElement3D::Basis(double u, double v, double w, vector<double>& Nt, vector<array<double,3>>& dNdt) const
{
	vector<double> Nu, Nv, Nw, dNdu, dNdv, dNdw;
	BezierPolyn(u,Nu,dNdu);
	BezierPolyn(v,Nv,dNdv);
	BezierPolyn(w,Nw,dNdw);
	Nt.resize(nbf);
	dNdt.resize(nbf);
	int i,j,k,loc(0);
	for(k=0; k<order; k++)
	{
		for(j=0; j<order; j++)
		{
			for(i=0; i<order; i++)
			{
				Nt[loc]=Nu[i]*Nv[j]*Nw[k];
				dNdt[loc][0]=dNdu[i]*Nv[j]*Nw[k];
				dNdt[loc][1]=Nu[i]*dNdv[j]*Nw[k];
				dNdt[loc][2]=Nu[i]*Nv[j]*dNdw[k];
				loc++;
			}
		}
	}
}

void BezierElement3D::Para2Phys(double u, double v, double w, double pt[3]) const
{
	vector<double> Nt;
	vector<array<double,3>> dNdt;
	Basis(u,v,w,Nt,dNdt);
	pt[0]=0.; pt[1]=0.; pt[2]=0.;
	for(int i=0;i<nbf;i++)
	{
		pt[0]+=pts[i][0]*Nt[i];
		pt[1]+=pts[i][1]*Nt[i];
		pt[2]+=pts[i][2]*Nt[i];
	}
}

void BezierElement3D::GeomMap(const array<double, 3>& u, array<double, 3>& pt) const
{
	vector<double> Nt;
	vector<array<double, 3>> dNdt;
	Basis(u[0], u[1], u[2], Nt, dNdt);
	pt[0] = 0.; pt[1] = 0.; pt[2] = 0.;
	for (int i = 0; i<nbf; i++)
	{
		pt[0] += pts[i][0] * Nt[i];
		pt[1] += pts[i][1] * Nt[i];
		pt[2] += pts[i][2] * Nt[i];
	}
}

//void BezierElement::SurfPointNormal4(double u, double v, array<double,3>& pt, array<double,3>& nm) const
//{
//	double Nx[25];
//	double dNdt[25][2];
//	Basis4(u,v,Nx,dNdt);
//	pt[0]=0.; pt[1]=0.; pt[2]=0.;
//	nm[0]=0.; nm[1]=0.; nm[2]=0.;
//	double nmtmp[2][3]={{0.,0.,0.},{0.,0.,0.}};
//	for(int i=0;i<25;i++)
//	{
//		pt[0]+=pts4[i][0]*Nx[i];
//		pt[1]+=pts4[i][1]*Nx[i];
//		pt[2]+=pts4[i][2]*Nx[i];
//		nmtmp[0][0]+=pts4[i][0]*dNdt[i][0];
//		nmtmp[0][1]+=pts4[i][1]*dNdt[i][0];
//		nmtmp[0][2]+=pts4[i][2]*dNdt[i][0];
//		nmtmp[1][0]+=pts4[i][0]*dNdt[i][1];
//		nmtmp[1][1]+=pts4[i][1]*dNdt[i][1];
//		nmtmp[1][2]+=pts4[i][2]*dNdt[i][1];
//	}
//	nm[0]=nmtmp[0][1]*nmtmp[1][2]-nmtmp[0][2]*nmtmp[1][1];
//	nm[1]=nmtmp[0][2]*nmtmp[1][0]-nmtmp[0][0]*nmtmp[1][2];
//	nm[2]=nmtmp[0][0]*nmtmp[1][1]-nmtmp[0][1]*nmtmp[1][0];
//	double len=sqrt(nm[0]*nm[0]+nm[1]*nm[1]+nm[2]*nm[2]);
//	nm[0]/=len; nm[1]/=len; nm[2]/=len;
//}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

BPatch3D::BPatch3D()
{
	deg[0] = 3; deg[1] = 3; deg[2] = 3;
	npt[0] = 0; npt[1] = 0; npt[2] = 0;
	kv.resize(3);
}

void BPatch3D::Initialize(const vector<vector<double>>& kv_in)
{
	unsigned int i;
	for (i = 0; i < kv.size(); i++)//==3
	{
		kv[i].clear();
		kv[i] = kv_in[i];
		npt[i] = kv[i].size() - deg[i] - 1;
	}
	cp.clear();
	cp.resize(npt[0] * npt[1] * npt[2]);
	for (i = 0; i < cp.size(); i++)
	{
		cp[i][0] = 0.; cp[i][1] = 0.; cp[i][2] = 0.;
	}
	wght.clear();
	wght.resize(cp.size(), 1.);
	pid.clear();
	pid.resize(cp.size(), -1);
}

void BPatch3D::BuildElement()
{
	ele.clear();
	unsigned int i, j, k;
	for (i = 0; i < IEN.size(); i++)
	{
		IEN[i].clear();
	}
	IEN.clear();
	for (k = 0; k < kv[2].size() - 1; k++)
	{
		for (j = 0; j < kv[1].size() - 1; j++)
		{
			for (i = 0; i < kv[0].size() - 1; i++)
			{
				if (kv[0][i] < kv[0][i + 1] && kv[1][j] < kv[1][j + 1] && kv[2][k] < kv[2][k + 1])
				{
					array<int, 3> tmp = { i, j, k };
					ele.push_back(tmp);
					vector<int> eien((deg[0] + 1)*(deg[1] + 1)*(deg[2] + 1));
					int loc(0);
					for (int k0 = k - deg[2]; k0 <= k; k0++)
					{
						for (int j0 = j - deg[1]; j0 <= j; j0++)
						{
							for (int i0 = i - deg[0]; i0 <= i; i0++)
							{
								eien[loc] = k0*npt[1] * npt[0] + j0*npt[0] + i0;
								loc++;
							}
						}
					}
					IEN.push_back(eien);
				}
			}
		}
	}
}

void BPatch3D::GlobalRefine(int dir)
{
	vector<array<double, 3>> cpnew;
	vector<double> wnew;
	vector<vector<double>> mat;
	vector<double> kvnew;
	int nptnew[3] = { npt[0], npt[1], npt[2] };
	unsigned int i, j, k, i0, j0, k0;
	for (i = 0; i < kv[dir].size() - 1; i++)
	{
		kvnew.push_back(kv[dir][i]);
		if (kv[dir][i] < kv[dir][i + 1])
		{
			kvnew.push_back((kv[dir][i] + kv[dir][i + 1]) / 2.);
		}
	}
	kvnew.push_back(kv[dir].back());
	nptnew[dir] = kvnew.size() - deg[dir] - 1;
	TMatrix(kv[dir], kvnew, deg[dir], mat);
	cpnew.resize(nptnew[0] * nptnew[1] * nptnew[2]);
	wnew.resize(cpnew.size(), 0.);
	for (i = 0; i < cpnew.size(); i++)
	{
		cpnew[i][0] = 0.; cpnew[i][1] = 0.; cpnew[i][2] = 0.;
	}
	for (i = 0; i < cp.size(); i++)
	{
		cp[i][0] = cp[i][0] * wght[i];
		cp[i][1] = cp[i][1] * wght[i];
		cp[i][2] = cp[i][2] * wght[i];
	}
	int locnew, locold;
	if (dir == 0)
	{
		//cout << "dir: " << dir << "\n";
		//cout << "nptnew: " << nptnew[dir] << "\n";
		for (i = 0; i < nptnew[dir]; i++)
		{
			//cout << i << "\n";
			for (i0 = 0; i0 < npt[dir]; i0++)
			{
				if (fabs(mat[i][i0]) > 1.e-10)
				{
					for (k = 0; k < npt[2]; k++)
					{
						for (j = 0; j < npt[1]; j++)
						{
							locnew = k*nptnew[1] * nptnew[0] + j*nptnew[0] + i;
							locold = k*npt[1] * npt[0] + j*npt[0] + i0;
							cpnew[locnew][0] += mat[i][i0] * cp[locold][0];
							cpnew[locnew][1] += mat[i][i0] * cp[locold][1];
							cpnew[locnew][2] += mat[i][i0] * cp[locold][2];
							wnew[locnew] += mat[i][i0] * wght[locold];
						}
					}
				}
			}
		}
	}
	else if (dir == 1)
	{
		for (j = 0; j < nptnew[dir]; j++)
		{
			for (j0 = 0; j0 < npt[dir]; j0++)
			{
				if (fabs(mat[j][j0]) > 1.e-10)
				{
					for (k = 0; k < npt[2]; k++)
					{
						for (i = 0; i < npt[0]; i++)
						{
							locnew = k*nptnew[1] * nptnew[0] + j*nptnew[0] + i;
							locold = k*npt[1] * npt[0] + j0*npt[0] + i;
							cpnew[locnew][0] += mat[j][j0] * cp[locold][0];
							cpnew[locnew][1] += mat[j][j0] * cp[locold][1];
							cpnew[locnew][2] += mat[j][j0] * cp[locold][2];
							wnew[locnew] += mat[j][j0] * wght[locold];
						}
					}
				}
			}
		}
	}
	else if (dir == 2)
	{
		for (k = 0; k < nptnew[dir]; k++)
		{
			for (k0 = 0; k0 < npt[dir]; k0++)
			{
				if (fabs(mat[k][k0]) > 1.e-10)
				{
					for (j = 0; j < npt[1]; j++)
					{
						for (i = 0; i < npt[0]; i++)
						{
							locnew = k*nptnew[1] * nptnew[0] + j*nptnew[0] + i;
							locold = k0*npt[1] * npt[0] + j*npt[0] + i;
							cpnew[locnew][0] += mat[k][k0] * cp[locold][0];
							cpnew[locnew][1] += mat[k][k0] * cp[locold][1];
							cpnew[locnew][2] += mat[k][k0] * cp[locold][2];
							wnew[locnew] += mat[k][k0] * wght[locold];
						}
					}
				}
			}
		}
	}

	npt[dir] = nptnew[dir];
	kv[dir].clear();
	kv[dir] = kvnew;
	wght.clear();
	wght = wnew;
	cp.clear();
	cp.resize(cpnew.size());
	for (i = 0; i < cpnew.size(); i++)
	{
		cp[i][0] = cpnew[i][0] / wnew[i];
		cp[i][1] = cpnew[i][1] / wnew[i];
		cp[i][2] = cpnew[i][2] / wnew[i];
	}
	pid.resize(cp.size(),-1);
}

void BPatch3D::BezierExtract(int eid, vector<vector<double>>& cmat, vector<int>& IENg)
{
	unsigned int i, j, k, i0, j0, k0;
	IENg.clear();
	IENg.resize(IEN[eid].size());
	for (i = 0; i < IEN[eid].size(); i++)
	{
		IENg[i] = pid[IEN[eid][i]];
	}
	for (i = 0; i < cmat.size(); i++)
	{
		cmat[i].clear();
	}
	cmat.clear();
	cmat.resize(IEN[eid].size(), vector<double>((deg[0] + 1)*(deg[1] + 1)*(deg[2] + 1), 0.));

	vector<double> ku0(kv[0].begin() + ele[eid][0] - deg[0], kv[0].begin() + (ele[eid][0] + 1) + deg[0] + 1);
	vector<double> kv0(kv[1].begin() + ele[eid][1] - deg[1], kv[1].begin() + (ele[eid][1] + 1) + deg[1] + 1);
	vector<double> kw0(kv[2].begin() + ele[eid][2] - deg[2], kv[2].begin() + (ele[eid][2] + 1) + deg[2] + 1);
	vector<double> ku1, kv1, kw1;
	array<double, 2> ktu = { kv[0][ele[eid][0]], kv[0][ele[eid][0] + 1] };
	array<double, 2> ktv = { kv[1][ele[eid][1]], kv[1][ele[eid][1] + 1] };
	array<double, 2> ktw = { kv[2][ele[eid][2]], kv[2][ele[eid][2] + 1] };
	vector<vector<double>> Tu, Tv, Tw;
	int iloc[3];
	BezierInsertKnots(ku0, ktu, ku1);
	BezierInsertKnots(kv0, ktv, kv1);
	BezierInsertKnots(kw0, ktw, kw1);
	TMatrix(ku0, ku1, 3, Tu);
	TMatrix(kv0, kv1, 3, Tv);
	TMatrix(kw0, kw1, 3, Tw);
	for (i0 = 0; i0 < ku1.size() - 1; i0++)
	{
		if (ku1[i0] == kv[0][ele[eid][0]] && ku1[i0 + 1] == kv[0][ele[eid][0] + 1])
		{
			iloc[0] = i0 - deg[0]; break;
		}
	}
	for (i0 = 0; i0 < kv1.size() - 1; i0++)
	{
		if (kv1[i0] == kv[1][ele[eid][1]] && kv1[i0 + 1] == kv[1][ele[eid][1] + 1])
		{
			iloc[1] = i0 - deg[1]; break;
		}
	}
	for (i0 = 0; i0 < kw1.size() - 1; i0++)
	{
		if (kw1[i0] == kv[2][ele[eid][2]] && kw1[i0 + 1] == kv[2][ele[eid][2] + 1])
		{
			iloc[2] = i0 - deg[2]; break;
		}
	}

	int loc0(0);//Bspline
	for (k0 = 0; k0 < deg[2] + 1; k0++)
	{
		for (j0 = 0; j0 < deg[1] + 1; j0++)
		{
			for (i0 = 0; i0 < deg[0] + 1; i0++)
			{
				int loc1(0);//Bezier
				for (int k1 = 0; k1 < deg[2] + 1; k1++)
				{
					for (int j1 = 0; j1 < deg[1] + 1; j1++)
					{
						for (int i1 = 0; i1 < deg[0] + 1; i1++)
						{
							cmat[loc0][loc1] = Tu[iloc[0] + i1][i0] * Tv[iloc[1] + j1][j0] * Tw[iloc[2] + k1][k0];
							loc1++;
						}
					}
				}
				loc0++;
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

MFNode::MFNode()
{
	coor[0] = 0.; coor[1] = 0.; coor[2] = 0.;
	a = 1.;
}

/////////////////////////

void RegularPatchBasis::Evaluate(double x,double a,double b)//input a(t-b)
{
	double t=a*(x-b);
	if(t<0.||t>1.)
	{
		for(int i=0;i<4;i++) {val[i]=0.; Dval[i]=0.;}
	}
	else
	{
		val[0]=(1.-3.*t+3.*t*t-t*t*t)/6.;
		val[1]=(4.-6.*t*t+3.*t*t*t)/6.;
		val[2]=(1.+3.*t+3.*t*t-3.*t*t*t)/6.;
		val[3]=t*t*t/6.;
		Dval[0]=a*(-3.+6.*t-3.*t*t)/6.;
		Dval[1]=a*(-12.*t+9.*t*t)/6.;
		Dval[2]=a*(3.+6.*t-9.*t*t)/6.;
		Dval[3]=a*t*t/2.;
		D2val[0]=a*(6.-6.*t)/6.;
		D2val[1]=a*(-12.+18.*t)/6.;
		D2val[2]=a*(6.-18.*t)/6.;
		D2val[3]=a*2.*t/2.;
	}
}

void RegularPatchBasis::Clear()
{
	for(int i=0;i<4;i++)
	{
		val[i]=0.;
		Dval[i]=0.;
		D2val[i]=0.;
	}
}

void Raw2Vtk_hex(string fn)
{
	unsigned int npt, nel;
	vector<array<double, 3>> pts;
	vector<array<int, 8>> cnct;
	double tmp;
	string fn1(fn + ".raw");
	ifstream fin;
	fin.open(fn1);
	if (fin.is_open())
	{
		fin >> npt >> nel;
		pts.resize(npt);
		cnct.resize(nel);
		for (unsigned int i = 0; i < npt; i++)
		{
			fin >> pts[i][0] >> pts[i][1] >> pts[i][2] >> tmp;
			//fin >> tmp;
		}
		for (unsigned int i = 0; i < nel; i++)
		{
			for (int j = 0; j < 8; j++) fin >> cnct[i][j];
			//fin >> tmp;
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fn1 << "!\n";
	}
	string fn2(fn + ".vtk");
	ofstream fout;
	fout.open(fn2);
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << pts.size() << " float\n";
		for (unsigned int i = 0; i<pts.size(); i++)
		{
			fout << pts[i][0] << " " << pts[i][1] << " " << pts[i][2] << "\n";
		}
		fout << "\nCELLS " << cnct.size() << " " << 9 * cnct.size() << '\n';
		for (unsigned int i = 0; i<cnct.size(); i++)
		{
			fout << "8 ";
			for (int j = 0; j<8; j++)
			{
				fout << cnct[i][j] << ' ';
			}
			fout << '\n';
		}
		fout << "\nCELL_TYPES " << cnct.size() << '\n';
		for (unsigned int i = 0; i<cnct.size(); i++)
		{
			fout << "12\n";
		}
		fout.close();
	}
	else
	{
		cerr << "Cannot open " << fn2 << "!\n";
	}
}

void Rawn2Vtk_hex(string fn)
{
	unsigned int npt, nel;
	vector<array<double, 3>> pts;
	vector<array<int, 8>> cnct;
	double tmp;
	string fn1(fn + ".rawn");
	ifstream fin;
	fin.open(fn1);
	if (fin.is_open())
	{
		fin >> npt >> nel;
		pts.resize(npt);
		cnct.resize(nel);
		for (unsigned int i = 0; i < npt; i++)
		{
			fin >> pts[i][0] >> pts[i][1] >> pts[i][2] >> tmp >> tmp >> tmp >> tmp;
		}
		for (unsigned int i = 0; i < nel; i++)
		{
			for (int j = 0; j < 8; j++) fin >> cnct[i][j];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fn1 << "!\n";
	}
	string fn2(fn + ".vtk");
	ofstream fout;
	fout.open(fn2);
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << pts.size() << " float\n";
		for (unsigned int i = 0; i<pts.size(); i++)
		{
			fout << pts[i][0] << " " << pts[i][1] << " " << pts[i][2] << "\n";
		}
		fout << "\nCELLS " << cnct.size() << " " << 9 * cnct.size() << '\n';
		for (unsigned int i = 0; i<cnct.size(); i++)
		{
			fout << "8 ";
			for (int j = 0; j<8; j++)
			{
				fout << cnct[i][j] << ' ';
			}
			fout << '\n';
		}
		fout << "\nCELL_TYPES " << cnct.size() << '\n';
		for (unsigned int i = 0; i<cnct.size(); i++)
		{
			fout << "12\n";
		}
		fout.close();
	}
	else
	{
		cerr << "Cannot open " << fn2 << "!\n";
	}
}

void Raw2Vtk_hex_delete(string fn, string fnd)
{
	unsigned int npt, nel;
	vector<array<double, 3>> pts;
	vector<array<int, 8>> cnct;
	double tmp;
	//string fn1(fn + ".raw");
	//ifstream fin;
	//fin.open(fn1);
	//if (fin.is_open())
	//{
	//	fin >> npt >> nel;
	//	pts.resize(npt);
	//	cnct.resize(nel);
	//	for (unsigned int i = 0; i < npt; i++)
	//	{
	//		fin >> pts[i][0] >> pts[i][1] >> pts[i][2] >> tmp;
	//		//fin >> tmp;
	//	}
	//	for (unsigned int i = 0; i < nel; i++)
	//	{
	//		for (int j = 0; j < 8; j++) fin >> cnct[i][j];
	//		//fin >> tmp;
	//	}
	//	fin.close();
	//}
	//else
	//{
	//	cerr << "Cannot open " << fn1 << "!\n";
	//}

	ReadVtk_hex(fn, pts, cnct);

	ifstream fin;
	vector<int> edel;
	fin.open(fnd);
	if (fin.is_open())
	{
		string stmp;
		fin >> stmp >> nel;
		edel.resize(nel);
		for (unsigned int i = 0; i < nel; i++)
		{
			fin >> edel[i];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fnd << "!\n";
	}

	vector<int> pflag(pts.size(), 0);
	vector<int> pid(pts.size(), -1);
	vector<int> eact(cnct.size(), 1);
	int npa(0), nea(0);
	for (unsigned int i = 0; i < edel.size(); i++)
	{
		eact[edel[i]] = 0;
	}
	for (int i = 0; i < cnct.size(); i++)
	{
		if (eact[i] == 1)
		{
			for (int j = 0; j < 8; j++)
			{
				pflag[cnct[i][j]] = 1;
			}
			nea++;
		}
	}
	for (int i = 0; i < pid.size(); i++)
	{
		if (pflag[i] == 1)
		{
			pid[i] = npa++;
		}
	}

	string fn2(fn + "_del.vtk");
	ofstream fout;
	fout.open(fn2);
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << npa << " float\n";
		for (unsigned int i = 0; i<pts.size(); i++)
		{
			if (pflag[i] == 1)
			{
				fout << pts[i][0] << " " << pts[i][1] << " " << pts[i][2] << "\n";
			}
			
		}
		fout << "\nCELLS " << nea << " " << 9 * nea << '\n';
		for (unsigned int i = 0; i<cnct.size(); i++)
		{
			if (eact[i] == 1)
			{
				fout << "8 ";
				for (int j = 0; j<8; j++)
				{
					fout << pid[cnct[i][j]] << ' ';
				}
				fout << '\n';
			}
		}
		fout << "\nCELL_TYPES " << nea << '\n';
		for (int i = 0; i < nea; i++)
		{
			fout << "12\n";
		}
		fout.close();
	}
	else
	{
		cerr << "Cannot open " << fn2 << "!\n";
	}
}

void ReadVtk_hex(string fn, vector<array<double, 3>>& pts, vector<array<int, 8>>& cnct)
{
	string fname(fn + ".vtk"), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i<4; i++) getline(fin, stmp);//skip lines
		fin >> stmp >> npts >> stmp;
		pts.resize(npts);
		for (int i = 0; i<npts; i++)
		{
			fin >> pts[i][0] >> pts[i][1] >> pts[i][2];
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		cnct.resize(neles);
		for (int i = 0; i<neles; i++)
		{
			fin >> itmp >> cnct[i][0] >> cnct[i][1] >> cnct[i][2] >> cnct[i][3] >>
				cnct[i][4] >> cnct[i][5] >> cnct[i][6] >> cnct[i][7];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}
}

void WriteVtk_hex(string fn, const vector<array<double, 3>>& pts, const vector<array<int, 8>>& cnct)
{
	string fn2(fn + ".vtk");
	ofstream fout;
	fout.open(fn2);
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << pts.size() << " float\n";
		for (unsigned int i = 0; i<pts.size(); i++)
		{
			fout << pts[i][0] << " " << pts[i][1] << " " << pts[i][2] << "\n";
		}
		fout << "\nCELLS " << cnct.size() << " " << 9 * cnct.size() << '\n';
		for (unsigned int i = 0; i<cnct.size(); i++)
		{
			fout << "8 ";
			for (int j = 0; j<8; j++)
			{
				fout << cnct[i][j] << ' ';
			}
			fout << '\n';
		}
		fout << "\nCELL_TYPES " << cnct.size() << '\n';
		for (unsigned int i = 0; i<cnct.size(); i++)
		{
			fout << "12\n";
		}
		fout.close();
	}
	else
	{
		cerr << "Cannot open " << fn2 << "!\n";
	}

}

void ReadRaw_hex(string fn, vector<array<double, 3>>& pts, vector<array<int, 8>>& cnct)
{
	int npt, nel;
	double tmp;
	string fn1(fn + ".raw");
	ifstream fin;
	fin.open(fn1);
	if (fin.is_open())
	{
		fin >> npt >> nel;
		pts.resize(npt);
		cnct.resize(nel);
		for (unsigned int i = 0; i < npt; i++)
		{
			fin >> pts[i][0] >> pts[i][1] >> pts[i][2] >> tmp;
		}
		for (unsigned int i = 0; i < nel; i++)
		{
			for (int j = 0; j < 8; j++) fin >> cnct[i][j];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fn1 << "!\n";
	}
}

void WriteRaw_hex(string fn, const vector<array<double, 3>>& pts, const vector<array<int, 8>>& cnct)
{
	string fn2(fn + "_hex.raw");
	ofstream fout;
	fout.open(fn2);
	if (fout.is_open())
	{
		fout << pts.size() << " " << cnct.size() << "\n";
		for (unsigned int i = 0; i<pts.size(); i++)
		{
			fout << pts[i][0] << " " << pts[i][1] << " " << pts[i][2] << " 0\n";
		}
		for (unsigned int i = 0; i<cnct.size(); i++)
		{
			for (int j = 0; j<8; j++)
			{
				fout << cnct[i][j] << ' ';
			}
			fout << '\n';
		}
		fout.close();
	}
	else
	{
		cerr << "Cannot open " << fn2 << "!\n";
	}

}

void Vtk2Raw_hex(string fn)
{
	vector<array<double, 3>> pts; 
	vector<array<int, 8>> cnct;
	ReadVtk_hex(fn, pts, cnct);
	WriteRaw_hex(fn, pts, cnct);
}

void ReadRaw_tri(string fn, vector<array<double, 3>>& pts, vector<array<int, 3>>& cnct)
{
	int npt, nel;
	double tmp;
	string fn1(fn + ".raw");
	ifstream fin;
	fin.open(fn1);
	if (fin.is_open())
	{
		fin >> npt >> nel;
		pts.resize(npt);
		cnct.resize(nel);
		for (unsigned int i = 0; i < npt; i++)
		{
			fin >> pts[i][0] >> pts[i][1] >> pts[i][2];
		}
		for (unsigned int i = 0; i < nel; i++)
		{
			for (int j = 0; j < 3; j++) fin >> cnct[i][j];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fn1 << "!\n";
	}
}

void WriteRaw_tri(string fn, const vector<array<double, 3>>& pts, const vector<array<int, 3>>& cnct)
{
	string fn2(fn + ".raw");
	ofstream fout;
	fout.open(fn2);
	if (fout.is_open())
	{
		fout << pts.size() << " " << cnct.size() << "\n";
		for (unsigned int i = 0; i<pts.size(); i++)
		{
			fout << pts[i][0] << " " << pts[i][1] << " " << pts[i][2] << "\n";
		}
		for (unsigned int i = 0; i<cnct.size(); i++)
		{
			for (int j = 0; j<3; j++)
			{
				fout << cnct[i][j] << ' ';
			}
			fout << '\n';
		}
		fout.close();
	}
	else
	{
		cerr << "Cannot open " << fn2 << "!\n";
	}

}

void ReadVtk_tri(string fn, vector<array<double, 3>>& pts, vector<array<int, 3>>& cnct)
{
	string fname(fn + ".vtk"), stmp;
	int npts, neles, itmp;
	ifstream fin;
	fin.open(fname);
	if (fin.is_open())
	{
		for (int i = 0; i<4; i++) getline(fin, stmp);//skip lines
		fin >> stmp >> npts >> stmp;
		pts.resize(npts);
		for (int i = 0; i<npts; i++)
		{
			fin >> pts[i][0] >> pts[i][1] >> pts[i][2];
		}
		getline(fin, stmp);
		fin >> stmp >> neles >> itmp;
		cnct.resize(neles);
		for (int i = 0; i<neles; i++)
		{
			fin >> itmp >> cnct[i][0] >> cnct[i][1] >> cnct[i][2];
		}
		fin.close();
	}
	else
	{
		cerr << "Cannot open " << fname << "!\n";
	}
}

void WriteVtk_tri(string fn, const vector<array<double, 3>>& pts, const vector<array<int, 3>>& cnct)
{
	string fn2(fn + ".vtk");
	ofstream fout;
	fout.open(fn2);
	if (fout.is_open())
	{
		fout << "# vtk DataFile Version 2.0\nHex test\nASCII\nDATASET UNSTRUCTURED_GRID\n";
		fout << "POINTS " << pts.size() << " float\n";
		for (unsigned int i = 0; i<pts.size(); i++)
		{
			fout << pts[i][0] << " " << pts[i][1] << " " << pts[i][2] << "\n";
		}
		fout << "\nCELLS " << cnct.size() << " " << 4 * cnct.size() << '\n';
		for (unsigned int i = 0; i<cnct.size(); i++)
		{
			fout << "3 ";
			for (int j = 0; j<3; j++)
			{
				fout << cnct[i][j] << ' ';
			}
			fout << '\n';
		}
		fout << "\nCELL_TYPES " << cnct.size() << '\n';
		for (unsigned int i = 0; i<cnct.size(); i++)
		{
			fout << "5\n";
		}
		fout.close();
	}
	else
	{
		cerr << "Cannot open " << fn2 << "!\n";
	}

}

void RepairConnect_tri(int npt, vector<array<int, 3>>& cnct, int ref)
{
	vector<vector<int>> p2f(npt);
	unsigned int i, j, k, k1;
	for (i = 0; i < cnct.size(); i++)
	{
		for (j = 0; j < 3; j++)
		{
			p2f[cnct[i][j]].push_back(i);
		}
	}
	vector<int> counted(cnct.size(), 0);
	counted[ref] = 1;
	vector<int> cur(1, ref);
	while (1)
	{
		vector<int> added(cnct.size(), 0);
		vector<int> next;//edge neighbors
		for (i = 0; i < cur.size(); i++)
		{
			int ednb[3] = { -1,-1,-1 };
			for (j = 0; j < 3; j++)
			{
				for (k = 0; k < p2f[cnct[cur[i]][j]].size(); k++)
				{
					int ic(p2f[cnct[cur[i]][j]][k]);
					if (ic != cur[i] && counted[ic] == 0)
					{
						for (k1 = 0; k1 < 3; k1++)
						{
							if (cnct[ic][k1] == cnct[cur[i]][(j + 1) % 3])
							{
								ednb[j] = ic;
								break;
							}
						}
					}
					if (ednb[j] != -1) break;
				}
			}
			for (j = 0; j < 3; j++)
			{
				if (ednb[j] != -1)
				{
					int loc(-1);
					for (k = 0; k < 3; k++)
					{
						if (cnct[ednb[j]][k] == cnct[cur[i]][j])
						{
							loc = k;
							break;
						}
					}
					if (loc != -1 && cnct[ednb[j]][(loc + 2) % 3] != cnct[cur[i]][(j + 1) % 3])
					{
						int itmp[3] = {cnct[ednb[j]][0],cnct[ednb[j]][2], cnct[ednb[j]][1]};
						for (k = 0; k < 3; k++)
						{
							cnct[ednb[j]][k] = itmp[k];
						}
					}
					counted[ednb[j]] = 1;
					if (added[ednb[j]] == 0)
					{
						next.push_back(ednb[j]);
						added[ednb[j]] = 1;
					}
				}
			}
		}
		cur.clear();
		if (next.size() == 0) break;
		cur = next;
	}
}

