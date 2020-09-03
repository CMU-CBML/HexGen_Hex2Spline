#ifndef BSPLINEBASIS_H
#define BSPLINEBASIS_H

#include <vector>

using namespace std;

class BSplineBasis
{
public:
	BSplineBasis();
	BSplineBasis(int deg,const vector<double>& v);
	void Set(int deg,const vector<double>& v);
	bool Check(int deg,const vector<double>& v);
	//int FindSpan(double par);
	//void Evaluate(double par,int deriv,vector<double>& val);
	void BasisFunction(int pos,double par,int deriv,vector<double>& val);
	//void test();
	
private:
	int p;//polynomial order
	int nbs;//# of basis functions
	vector<double> kv;//knot vector
};

#endif