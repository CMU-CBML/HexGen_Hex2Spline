#include "BSplineBasis.h"
#include <iostream>
//#include <fstream>
//#include <string>

using namespace std;

typedef unsigned int uint;

BSplineBasis::BSplineBasis()
{
	p=3; nbs=1;
	double tmp[5]={0.,1.,2.,3.,4.};
	kv.assign(tmp,tmp+5);
}

BSplineBasis::BSplineBasis(int deg,const vector<double>& v)
{
	if(Check(deg,v))
	{
		p=deg; kv=v; nbs=v.size()-deg-1;
	}
	else
	{
		p=3; nbs=1;
		double tmp[5]={0.,1.,2.,3.,4.};
		kv.assign(tmp,tmp+5);
	}
}

void BSplineBasis::Set(int deg,const vector<double>& v)
{
	if(Check(deg,v))
	{
		p=deg; kv=v; nbs=v.size()-deg-1;
	}
	else
	{
		p=3; nbs=1;
		double tmp[5]={0.,1.,2.,3.,4.};
		kv.assign(tmp,tmp+5);
	}
}

bool BSplineBasis::Check(int deg,const vector<double>& v)
{
	if(deg<0)
	{
		cerr<<"Wrong degree!\n";
		return false;
	}
	if(v.size()<=p+1)
	{
		cerr<<"Length of knot vector is wrong!\n";
		return false;
	}
	else
	{
		int flag(0);
		for(uint i=0;i<v.size()-1;i++)
		{
			if(v[i]>v[i+1])
			{
				flag=1;
				break;
			}
		}
		if(flag==1)
		{
			cerr<<"Knot vector is not in increasing order!\n";
			return false;
		}
		else
		{
			int rep(0);
			//check first knot
			for(uint i=0;i<v.size();i++)
			{
				if(v[i]==v[0]) rep++;
			}
			if(rep>p+1)
			{
				cerr<<"First knot repeated more than "<<p<<"+1 times\n";
				return false;
			}
			//check last knot
			rep=0;
			for(int i=v.size()-1;i>=0;i--)
			{
				if(v[i]==v.back()) rep++;
			}
			if(rep>p+1)
			{
				cerr<<"Last knot repeated more than "<<p<<"+1 times\n";
				return false;
			}
			//check interior knots
			for(uint i=1;i<v.size()-1;i++)
			{
				rep=0;
				for(uint j=i;j<v.size()-1;j++)
				{
					if(v[j]==v[i]) rep++;
				}
				if(rep>p)
				{
					cerr<<"Interior knots repeated more than "<<p<<" times\n";
					return false;
				}
			}
		}
	}
	return true;
}

//int BSplineBasis::FindSpan(double par)
//{
//	if(par<kv.front() || par>kv.back())
//	{
//		cerr<<"The parametric value is not valid!\n";
//		return 0;
//	}
//	int pos(-1);
//	if(par==kv.back())
//		pos=kv.size()-2;
//	else
//	{
//		for(uint i=0;i<kv.size()-1;i++)
//		{
//			if(par>=kv[i] && par<kv[i+1])
//			{
//				pos=i;
//				break;
//			}
//		}
//	}
//	if(pos<0)
//	{
//		cerr<<"Wrong span!\n";
//		return 0;
//	}
//	return pos;
//}
//
//void BSplineBasis::Evaluate(double par,int deriv,vector<double>& val)
//{
//	if(deriv<0 || deriv>p)
//	{
//		cerr<<"Wrong order of direvative\n";
//		return;
//	}
//	int i=FindSpan(par);
//	val.resize((deriv+1)*nbs,0.);
//
//	//evaluate
//	vector<vector<double>> ders(deriv+1,vector<double>(p+1,0.));
//	vector<vector<double>> ndu(p+1,vector<double>(p+1,0.));
//	vector<vector<double>> a(2,vector<double>(p+1,0.));
//	vector<double> left(p+1,0.);
//	vector<double> right(p+1,0.);
//	double saved=0.,temp=0.,d;
//	int j,r,s1,s2,k,rk,pk,j1,j2;
//	for(j=1;j<=p;j++)
//	{
//		left[j]=par-kv[i+1-j];
//		right[j]=kv[i+j]-par;
//		saved=0.;
//		for(r=0;r<j;r++)
//		{
//			ndu[j][r]=right[r+1]+left[j-r];
//			temp=ndu[r][j-1]/ndu[j][r];
//			ndu[r][j]=saved+right[r+1]*temp;
//			saved=left[j-r]*temp;
//		}
//		ndu[j][j]=saved;
//	}
//	for(j=0;j<=p;j++)
//	{
//		ders[0][j]=ndu[j][p];
//	}
//	for(r=0;r<=p;r++)
//	{
//		s1=0; s2=1;
//		a[0][0]=1.;
//		for(k=1;k<=deriv;k++)
//		{
//			d=0.;
//			rk=r-k; pk=p-k;
//			if(r>=k)
//			{
//				a[s2][0]=a[s1][0]/ndu[pk+1][rk];
//				d=a[s2][0]*ndu[rk][pk];
//			}
//			if(rk>=-1) j1=1;
//			else j1=-rk;
//			if(r-1<=pk) j2=k-1;
//			else j2=p-r;
//			for(j=j1;j<=j2;j++)
//			{
//				a[s2][j]=(a[s1][j]-a[s1][j-1])/ndu[pk+1][rk+j];
//				d+=a[s2][j]*ndu[rk+j][pk];
//			}
//			if(r<=pk)
//			{
//				a[s2][k]=-a[s1][k-1]/ndu[pk+1][r];
//				d+=a[s2][k]*ndu[r][pk];
//			}
//			ders[k][r]=d;
//			j=s1; s1=s2; s2=j;
//		}
//	}
//	r=p;
//	for(k=1;k<=deriv;k++)
//	{
//		for(j=0;j<=p;j++) ders[k][j]*=r;
//		r*=(p-k);
//	}
//}

void BSplineBasis::BasisFunction(int pos,double par,int deriv,vector<double>& val)
{
	val.clear();
	val.resize(deriv+1,0.);
	if(pos<0 || pos>=nbs)
	{
		cerr<<"Wrong basis functions ID!\n";
		return;
	}
	if(par<kv.front() || par>kv.back())
	{
		//cerr<<"Wrong parametric value!\n";
		return;
	}
	if(deriv>1)
	{
		cerr<<"Higher order of derivatives not available now!\n";
		return;
	}

	vector<double> N0(p+1,0.);
	double left,right,d1=0.;
	int i,j,q,flag;
	//find span
	//for(i=pos+1;i<=pos+p;i++)
	//{
	//	if(par==kv[i] && par>kv[i-1] && i+1<=pos+p && kv[i]==kv[i+1])
	//	{
	//		N0[i-pos-1]=1.;
	//	}
	//}
	//if(par==kv.back())
	//{
	//	int open(1);
	//	int end(kv.size()-1);
	//	for(int i0=0; i0<p; i0++)
	//	{
	//		if(kv[end-i0]!=kv[end-i0-1])
	//		{
	//			open=0; break;
	//		}
	//	}
	//	if(open==1)
	//	{
	//		val[0]=1.;
	//		if(deriv==1)
	//			val[1]=double(p)/(kv[pos+p]-kv[pos]);
	//		return;
	//	}
	//	else
	//	{
	//		return;
	//	}
	//}
	for(i=pos;i<=pos+p;i++)
	{
		if(par>=kv[i] && par<kv[i+1] || (par>kv[i] && par==kv[i+1] && par==kv.back()))
		{
			N0[i-pos]=1.;
			break;
		}
	}
	//recursive
	for(q=1;q<=p;q++)
	{
		j=p+1-q;
		for(i=0;i<j;i++)
		{
			if(kv[pos+i+q]==kv[pos+i]) left=0.;
			else left=(par-kv[pos+i])/(kv[pos+i+q]-kv[pos+i]);
			if(kv[pos+i+q+1]==kv[pos+i+1]) right=0.;
			else right=(kv[pos+i+q+1]-par)/(kv[pos+i+q+1]-kv[pos+i+1]);
			if(deriv==1 && q==p)//first derivative
			{
				double left1,right1;
				if(kv[pos+i+q]==kv[pos+i]) left1=0.;
				else left1=q/(kv[pos+i+q]-kv[pos+i]);
				if(kv[pos+i+q+1]==kv[pos+i+1]) right1=0.;
				else right1=q/(kv[pos+i+q+1]-kv[pos+i+1]);
				d1=left1*N0[0]-right1*N0[1];
			}
			//cout<<left<<" "<<right<<"\n";
			//cout<<N0[i]<<" "<<N0[i+1]<<"\n";
			N0[i]=left*N0[i]+right*N0[i+1];
			//cout<<N0[i]<<"\n";
			//getchar();
		}
	}
	val[0]=N0[0];
	if(deriv==1)
		val[1]=d1;
}

//void BSplineBasis::test()
//{
//	int nsmp=5,j;
//	vector<double> smp((kv.size()-1)*nsmp+1);
//	vector<double> fv((kv.size()-1)*nsmp+1);
//	unsigned int i,count=0;
//	for(i=0;i<kv.size()-1;i++)
//	{
//		for(j=0;j<nsmp;j++)
//		{
//			smp[count]=kv[i]+j*(kv[i+1]-kv[i])/nsmp;
//			count++;
//		}
//	}
//	smp[smp.size()-1]=kv.back();
//
//	vector<double> val;
//	for(int id=0;id<nbs;id++)
//	{
//		for(i=0;i<smp.size();i++)
//		{
//			BasisFunction(id,smp[i],1,val);
//			fv[i]=val[1];
//		}
//		ofstream fout;
//		fout.open("BasisFuncTest"+to_string(long long(id))+".txt");
//		if(fout.is_open())
//		{
//			for(i=0;i<smp.size();i++)
//			{
//				fout<<smp[i]<<' '<<fv[i]<<'\n';
//			}
//			fout.close();
//		}
//	}
//}
