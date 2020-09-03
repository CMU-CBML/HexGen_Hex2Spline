#ifndef KNOT_INSERTION_H
#define KNOT_INSERTION_H

#include <vector>
#include <array>

using namespace std;

void InitialT(const vector<double>& k1,const vector<double>& k2,vector<vector<double>>& T);
void KnotInsert(const vector<double>& k1,const vector<double>& k2,const vector<vector<double>>& T1,int q,vector<vector<double>>& T2);
void TMatrix(const vector<double>& k1,const vector<double>& k2,int p,vector<vector<double>>& T);

void InsertKnots(const vector<double>& kv,const vector<double>& knots,vector<double>& kv1);
void BezierInsertKnots(vector<double>& kv,array<double,2>& knots,vector<double>& kv1);
void BezierInsertKnots4(vector<double>& kv,array<double,2>& knots,vector<double>& kv1);

void DegreeElevate(vector<vector<double>>& demat);//from degree-3 to degree-4

void KnotVectorSubdiv(const vector<double>& k1, vector<double>& k2);

#endif