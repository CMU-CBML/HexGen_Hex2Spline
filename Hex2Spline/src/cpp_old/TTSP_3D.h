#ifndef TTSP_3D_H
#define TTSP_3D_H

#include <vector>
#include <utility>
#include <string>
#include "atlstr.h"
//#include "T_mesh.h"
#include "BasicDataStructure.h"
#include "Laplace.h"
using namespace std;

class TruncatedTspline_3D
{
public:
	TruncatedTspline_3D();
	void CreateUniformCube(string fn);
	void VisualizeControlMesh(string fn);
	double PartionOfUnity(int eid, const array<double, 3>& u);
	void CollectActives();

	void VisualizeTMesh(string fn);
	void VisualizeBezierMesh(vector<BezierElement3D> bzmesh, string fn);
	void VisualizeFaceMesh(string fn);
	bool CheckSubKnotVector(const double ku1[5], const double kv1[5], const double ku2[5], const double kv2[5]);//check whether v1 is a subvector of v2
	bool CheckSubKnotVector(const array<double,5>& ku1, const array<double,5>& kv1, const array<double,5>& ku2, const array<double,5>& kv2);

	void InitialConnect();
	void InitialRotate();
	void getElementRotate(vector<Matrix3d>& mat);
	void getCornerRotate(int hxid, int uvw_loc[3], int uvw_ref[3], Matrix3d& mat);
	void UpdateConnect();
	void FindEdgeTopoDirec();
	void EdgeConnect(int ed0, int pid0, int& type1, int& ed1);
	//void FindFaceTopoDirec();
	void FaceConnect(int fc0, int ed0, int& type1, int& fc1);
	void setReference();
	void FindKnotInterval();
	//void UpdateKnotInterval();
	void ShootRay(int pid, int edid, double kv[4], int& trun_flag);//start from point pid with the edge edid
	void SetLocalCoorSystem();
	void Find_Neighbor_Rot(int hxid, int pid, Matrix3d& rot);
	void getElementRotate_Unit(int loc, Matrix3d& rot);
	//void FindIEN_Unstruct();
	void FindIEN_PatchKV();
	void Update_IEN();
	void FindNextRing(const vector<int>& pr0, const vector<int>& er0, vector<int>& pr1, vector<int>& er1, vector<int>& pr1_pref, vector<int>& pr1_eref);
	void TranslateLCS(int pref,Matrix4d& lcs_ref,int eid,int pid,Matrix4d& lcs);
	void getLCS_inverse(const Matrix4d& mat_in, Matrix4d& mat_out);
	void FindLocalKnotVector(int pid, const Matrix4d& lcs, array<double, 5>& ku, array<double, 5>& kv, array<double, 5>& kw);
	bool CheckSupport(const array<double, 2>& u, const array<double, 2>& v, const array<double, 2>& w, const array<double, 5>& ku, const array<double, 5>& kv, const array<double, 5>& kw);
	
	void SetBezierMatIrrPatch(int eid);
	void AllBezierPatch();
	void Para2Physical(int eid, const array<double, 3>& u, array<double, 3>& pt);
	void UpdatePatchCP_Unstruct(int eid);

	void Truncation();
	bool CheckSubKnotVector(const array<double, 5>& ku1, const array<double, 5>& kv1, const array<double, 5>& kw1,
		const array<double, 5>& ku2, const array<double, 5>& kv2, const array<double, 5>& kw2);
	//void FindChildren();

	void ElementBasis(int eid, const array<double, 3>& u, vector<double>& Nt, vector<array<double, 3>>& dNdt);
	//void ElementBasis_Bezier(int eid, const array<double,3>& u, vector<double>& Nt, vector<array<double, 3>>& dNdt);
	void ElementBasis_Regular(int eid, const array<double, 3>& u, vector<double>& Nt, vector<array<double, 3>>& dNdt);
	void ElementBasis_Irregular(int eid, const array<double, 3>& u, vector<double>& Nt, vector<array<double, 3>>& dNdt);
	void FindPatchKnotVector_Irr(int eid, vector<array<double,5>>& patch_ku, vector<array<double,5>>& patch_kv);
	void SurfacePointMap(int eid, double u, double v, array<double,3>& pt, array<double,3>& norm);

	void ElementSubdivide_8(int eid);
	void ElementSubdivide_5(int eid, int dir);
	void ElementSubdivide_4(int eid, int dir);
	void ElementSubdivide_2(int eid, int dir);//add a yz plane if dir=0, xz for dir=1, xy for dir=2

	void FaceSubdivision_4(int fcid);
	void FaceSubdivision_2(int fcid, int dir);
	void FaceSubdivision_24(int fcid, int dir);//subdivide a face that has been split
	void EdgeSubdivision_2(int edid);

	void SolidFaceDirection(int eid, int fcid, int& dir, int& pos);
	void EdgeIndex_in_Face(Face3D& fc, const vector<int>& edid);
	void EdgeFaceIndex_in_Solid(Element3D& hex, const vector<int>& edid, const vector<int>& fcid);
	void Index_Direction(int eid, int vloc[8], int edloc[12], int fcloc[6], int dir);

	void Construct_NewEdge(int pt[2], int lev, double len);
	void Construct_NewFace(int cnct[4], int edid[4], double lv, double len);
	void Construct_BaseFace_Subdv(int fcid, int dir, int& ed_base, int fc_base[2]);
	void Construct_BaseFace_Bisct(int fcid, int dir);

	void IdentifyTarget(vector<int>& target);
	void IdentifyAddition(vector<int>& target);
	bool IdentifyAddition();
	void IdentifyPropagate(vector<int>& id_in, vector<int>& id_out);
	void RefineTopology();
	void RefineGeometry();
	void RefineGeometry_v1();

	void CalPatchCP_Regular(int eid);
	bool CheckFullRefine(const vector<array<double, 3>>& spt, const array<double, 5>& motu, const array<double, 5>& motv, const array<double, 5>& motw,
		const vector<array<double, 5>>& chdu, const vector<array<double, 5>>& chdv, const vector<array<double, 5>>& chdw, const vector<double>& coef);
	//void CheckTruncation();
	void UpdateGeometry(vector<int>& rfid);

	bool Identify_FaceFace(int eid);
	bool Identify_FaceEdge(int eid);
	bool Template_FaceFace();
	bool Template_FaceEdge();

	void RefineTest_1();
	void RefineTest_2();

	void Topo_Refine_Unstruct(const vector<int>& rfid, const vector<int>& rftype, vector<int>& rfid_more, vector<int>& rftype_more);
	void Geom_Refine_Unstruct(const vector<int>& rfid, const vector<int>& rftype);

	void BezierExtract_Unstruct(vector<BezierElement2D>& bzmesh);
	void BezierElementExtract_Unstruct(int eid,vector<BezierElement2D>& bzmesh);
	void BezierElementExtract_Unstruct_Irr(int eid,vector<BezierElement2D>& bzmesh);
	void BezierUnit_Unstruct(int eid,array<double,4>& kts,vector<BezierElement2D>& bzmesh);
	void BezierUnit_Unstruct_Irr(int eid,array<double,4>& kts,vector<BezierElement2D>& bzmesh);
	void BezierUnit_Unstruct_Trun(int eid,array<double,4>& kts,vector<BezierElement2D>& bzmesh);
	void BezierUnit_Unstruct_Irr_Trun(int eid,array<double,4>& kts,vector<BezierElement2D>& bzmesh);
	void BezierFinder_Unstruct(int eid,vector<array<double,4>>& be);
	void BezierVTK_Unstruct(string fn,vector<BezierElement2D>& bzmesh);
	void BezierControlMesh_Unstruct(string fn,vector<BezierElement2D>& bzmesh);

	//void run_surf_XP(string fn);
	//void runXP_Laplace();
	void run(string fn);
	void runh(string fn);

	void SetBezier3TranMat(int N, vector<vector<double>>& bmat);
	void VisualizeSolidVTK(string fn);

	void SetProblem_surf_XP(string fn);
	void SetLshapeProblem_XP(string fn);
	void SetProblem(string fn);
	void SetSharpFeature();
	void SetDomain();

	void Identify_Invalid_Elements(vector<int>& rid);
	void Identify_More_Elements(vector<int>& rid);
	void SetBounaryPoints(vector<int>& pid, vector<double>& disp, int& ncp);

	void Refine_Surf_Select_0(vector<int>& rfid, vector<int>& rftype);
	void Refine_Surf_Test_0();

	void ElementDomain();
	void PatchRefine_Regular(int lev, int eid);
	void BisectKnotInterval(const array<double, 5>& kv_in, vector<double>& kv_out);
	void PatchRefine_Irregular(int lev, int eid);
	void CatmullClark(int lev, int eid);
	void PatchRefine_Boundary(int lev, int eid);
	void ConstructBezierBasis(int lev, int eid);
	void ConstructBezierBasis_1(int lev, int eid);
	void ConstructBezierBasis_Boundary(int lev, int eid);
	void ConstructBezierBasis_Feature(int lev, int eid);
	void ConstructConnect(int lev);
	void ConstructFaceEdge(int lev, int eid);
	void Refine_Ghost(const vector<array<int, 2>>& rfid);

	void Selection(int lev);
	void SetSupport(int lev);
	void Select();
	void Truncate(int lev);
	void Basis_Regular(int lev, int eid, const array<double, 3>& u, vector<double>& Nt, vector<array<double, 3>>& dNdt);
	void Basis_Irregular(int lev, int eid, const array<double, 3>& u, vector<double>& Nt, vector<array<double, 3>>& dNdt);
	void Identify_Pseudo(vector<array<int,2>>& rfid);
	void Identify_Pseudo_1(vector<array<int, 2>>& rfid);
	void Identify_Test_1(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);
	void Identify_Test_2(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);
	void Identify_Test_3(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);
	void Identify_Test_4(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);//boundary
	void Refine(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);
	double BasisSum(int lev, int eid, const array<double,3>& u);

	void OutputCM(int lev, string fn);
	void OutputFace(int lev, string fn);
	void OutputEdge(int lev, string fn);
	void OutputGeom(int lev, string fn);
	void GeomMap(int lev, int eid, const array<double,3>& u, array<double,3>& pt);
	void GeomMap_Bezier(int lev, int eid, const array<double, 3>& u, array<double, 3>& pt);
	void OutputGeom_All(string fn);

	void GeomMap_Lev(int lev, int eid, const array<double, 3>& u, array<double, 3>& pt);
	double BasisSum_Lev(int lev, int eid, const array<double, 3>& u);
	void Basis_Regular_Lev(int lev, int eid, const array<double, 3>& u, vector<double>& Nt, vector<array<double, 3>>& dNdt);
	void Basis_Irregular_Lev(int lev, int eid, const array<double, 3>& u, vector<double>& Nt, vector<array<double, 3>>& dNdt);

	void AllBezierLev(int lev);
	void AnalysisInterface_Elastic(vector<BezierElement3D>& bzmesh, vector<int>& DrchBC, vector<double>& gh);
	void AnalysisInterface_Poisson(vector<BezierElement3D>& bzmesh, vector<int>& DrchBC, vector<double>& gh);//only for cube domain [0,1]^3
	void AnalysisInterface_Poisson_1(vector<BezierElement3D>& bzmesh, vector<int>& DrchBC, vector<double>& gh);//for arbitrary shapes
	void AnalysisInterface_Laplace(const vector<array<int, 2>>& pbc, const vector<double>& pdisp, vector<BezierElement3D>& bzmesh, vector<int>& DrchBC, vector<double>& gh);
	void AnalysisInterface_LeastSquare(vector<BezierElement3D>& bzmesh, vector<int>& DrchBC, vector<double>& gh);
	double SpecifyDirichBC(double x[3]);
	double SpecifyDirichBC_2(double x[3]);
	double SpecifyDirichBC_3(double x[3]);
	double SpecifyDirichBC_4(double x[3]);
	double SpecifyDirichBC_5(double x[3]);
	double SpecifyDirichBC_6(double x[3]);
	double SpecifyDirichBC_7(double x[3]);
	double SpecifyDirichBC_8(double x[3]);
	double SpecifyDirichBC_9(double x[3]);
	void Identify_Poisson(const vector<array<double, 2>>& ehid, const vector<double>& err, vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);
	void Identify_Poisson_1(const vector<array<double, 2>>& ehid, const vector<double>& err, vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);//this is new
	void Identify_Laplace(const vector<array<double, 2>>& ehid, const vector<double>& err, vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);
	void Identify_LeastSquare(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);//sphere
	void Identify_LeastSquare_Line(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);//line
	void SetInitialBC(vector<array<int, 2>>& ebc, vector<double>& edisp, vector<array<int, 2>>& pbc, vector<double>& pdisp);
	void SetBC(vector<array<int, 2>>& ebc, vector<double>& edisp, vector<array<int, 2>>& pbc, vector<double>& pdisp);
	void FittingBC(vector<int>& DrchBC, vector<double>& gh);

	void BezierPoints_ini();
	void BezierPoints_Refine(int lev, int pos, const vector<MatrixXd>& bsmat);//input is the father element
	void BezierSubdivMatrix(vector<MatrixXd>& bsmat);
	void OutputRefineID(string fn, const vector<array<int,2>>& rfid, const vector<array<int,2>>& gst);
	void InputRefineID(string fn, vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);
	void InputRefineID_Manual(string fn, vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);
	void GetRemoveRegion(double xy[3][2]);
	void OutputRemoveCM(string fn, double xy[3][2]);
	void CreateRemoveView(string fn_in, string fn_out);

	void OutputBasis(int lev, int pid, string fn);
	void PillowCube(string fn);
	void VisualizeBezier(const vector<BezierElement3D>& bzmesh, string fn);
	void InputCheck(string fn);
	void MeshRepair(string fn_in, string fn_out);
	void ReportXP();

	void GlobalRefine(int niter);//this is just using functions of local refinement, not efficient
	double MaxElementSize(const vector<BezierElement3D>& bzmesh);

	//true global refinement, can also be included in local refinement
	void InitializeMesh(string fn);
	void SetSharpFeature_1(double tol);//use tmesh, tmface, tmedge
	void SetSharpFeature_localrefine(double tol);//use hmesh, hface, hedge
	void Global_Subdivide();
	void Global_Subdivide_Simple();
	void BuildBasisFunction();
	void OutputCM(string fn);
	void OutputFace(string fn);
	void OutputFace_Boundary(string fn);
	void OutputEdge(string fn);
	void run_Global(int nref, string fn);

	void SetDomainRange(double rg[3][2], double nm[3], double& a);

	//reimplement local refinement for efficiency
	void Refine_eff(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);
	void PatchRefine_eff(int lev, int eid);

	//for pipeline paper
	void PipelineDataProcess(string fn_in, string fn_out);
	void PipelineBezierExtract(vector<BezierElement3D>& bzmesh);

	//B-splines for test
	void CreateBsplines(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void CreateBsplines(int nex[3]);
	void Bsplines_Refine();//global refinement
	void Bspline_BezierExtract(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void Bspline_BezierExtract_fit(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void FittingBC_Bsplines(const vector<BezierElement3D>& bzall, vector<int>& DrchBC, vector<double>& gh);
	void SetDomainRange_Bsplines(double rg[3][2], double nm[3], double& a);

	//solid T-splines with pillowed layers as open knot vectors
	void ReadHexMeshVTK(string fn);
	void DeleteDuplicatePoint(string fn);
	void BuildSolidTspline();//very limited for hex mesh
	void Build_Edge_Face();
	void MeshTopology();
	void BuildElementBasis();//IEN and knot vectors
	void SetRegularPatch(int eid);
	void UpdateIEN_Irr(int eid);
	void Refine_Global();
	void RecordTmeshInfo();
	void Refine_Element_Reg(int eid);
	//void Refine_Element_Bnd(int eid);
	void Refine_Element_Bnd_4(int eid, int dir);//knot interval in one direction is 0
	void Refine_Element_Bnd_2(int eid, int dir);//knot intervals in two directions are 0
	//void Refine_Element_Irr(int eid);
	void Find_FaceEdge_Index(int hxid[8], int fcid[36], int edid[54]);
	void Refine_Global_Topology();//connectivity, IEN, etc
	void Refine_Global_Geometry();//calculate coordinates of control points
	void Refine_Element_Int_Topology(int eid);//interior elements including regular and irregular
	void Refine_Element_Bnd4_Topology(int eid, int dir);//boundary into 4
	void Refine_Element_Bnd2_Topology(int eid, int dir);//boundary into 2
	void Refine_Element_Bnd0_Topology(int eid, int dir);//boundary unchanged, like one at corner
	void Refine_Element_Reg_Geometry(int eid);//including boundary elements
	void Refine_Element_Irr_Geometry(int eid);

	void Clear_current();
	void Clear_tmesh();
	void Clear_tmface();
	void Clear_tmedge();
	void Clear_tmp();
	void Clear_tmtmp();
	void Clear_tmfctmp();
	void Clear_tmedtmp();

	//3D convergence study by coupling
	void Run_Coupling(string fn, int nrf, vector<BezierElement3D>& bzmesh, vector<int>& IDBC);
	void BuildSplines_Unstruct();
	void BuildElementSplines_Interior(int eid);
	void BuildElementSplines_Boundary(int eid);
	void BezierExtract3D_Unstruct(vector<BezierElement3D>& bzmesh);
	void AddMeshfreeNodes(vector<BezierElement3D>& bzmesh);
	void Set_RKPM(vector<BezierElement3D>& bzmesh);
	void FindNeighbor_2r(int eid, vector<int>& hx2r);
	void FindMeshfreeNodes_glb(BezierElement3D& bzel);
	void Find_IDBC(vector<int>& IEN);
	void OutputMesh_Coupling(const vector<BezierElement3D>& bzmesh, string fn);
	//void GeomMap_Bezier(int eid, const array<double, 3>& u, array<double, 3>& pt);


	//3D coupling with Bezier
	void Run_Coupling_Bezier(string fn, int nrf, vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void Convert2Bezier();
	void Global_Subdivide_wBezier();
	void BezierPatch_Refine(vector<Element3D>& hxnew);
	void FindBezierIEN(const vector<Element3D>& hxnew);
	void ClearRefineMemory(vector<Element3D>& hxnew, vector<Face3D>& fcnew, vector<Edge3D>& ednew);
	void BezierExtract_BezierCoupling(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void RescaleDomain();
	void ReportDiffEleNumber();

	//Comparison
	void Run_Coupling_Comp(string fn, int ngrf, int nlrf, vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void BezierExtract_Comp(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void UserDefineIDBC(vector<int>& IDBC, vector<double>& gh);
	void BuildInitialEdges();
	void BezierExtract_Comp_Fitting(vector<BezierElement3D>& bzmesh, vector<int>& IDBC);
	void Fitting2InptHex(string fn_in, string fn_out);
	void Run_Coupling_Comp_Fit(string fn, int nrf, vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);

	void  Initialize_AdaptFit(string fn);
	void BezierExtract_AdaptFit(vector<BezierElement3D>& bzmesh);
	void StrongDirichletBC_Elasticity_AdaptFit(vector<int>& IDBC, vector<double>& gh);
	void Identify_Elasticity_AdaptFit(const vector<array<double, 2>>& eh, const vector<double>& vs, 
		vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);
	void VisualizeBC(string fn, const vector<int>& IDBC, const vector<double>& gh);

	//B-splines-like unstructured
	void Run_UBsplines(string fn, int nrf, vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void InitializeMesh_UB(string fn);
	void ReadVtk_Hex(string fn);
	void InitialConnect_UB();
	void SetType_UB();

	//quality
	int IsPillowNeeded();
	void Pillow(int nlayer = 2);
	void Smoothing(int nStep, double stepSize);//move interior points towards mass center

	/// Juelin's optimization
	double DeterminantOfThreePoint(double *p1, double *p2, double *p3);
	void SetHexQuality(int hexindex);
	void HexMeshQuality();
	void MaxMinHexJacobian(int vert, double *displacement);
	double Distance(double *A, double *B);
	double DotProduct(double p1[], double p2[]);
	void CrossProduct1(double *U, double *V, double *W);
	void CrossProduct2(double *u, double *v, double *w);
	
	void Jordan(double *A, double *f, double *b, int ir, int k, int m, int n);
	void SimplexMethod(double *A, double *b, double *c, int m, int n, double *sol);
	void ImproveHexJacobian(int Iter, double LBound);
	double VolumeOfHexahedron(double *point);
	double VolumeOfHex(int hexindex);
	double Compute_hexh();
	double Compute_Ji(int vert);
	double VolumeEnergy(double h, double lambda);
	void Volume_regularization(double Lambda, double Tau, int Iter);

	/// Angran's smoothing
	double InnerProduct(vector<double> a, vector<double>b);
	void ComputeMeanCurvature(int pid, double H[3]);
	void ComputeCurvatureSharp(int pid, double k[3]);
	void ComputePointNormal(int pid, double nm[3]);
	void GetQuadInfo_Angran(int fcid, int ploc, double& area, double center[3], double nm[3], double alpha, double beta, double gama, double dA[3]);
	int SmoothingPoint_Angran(int pid, double stepSize);
	int SmoothingPointBoundary_Angran(int pid, double stepSize);
	int SmoothingPointBoundarySharp_Angran(int pid, double stepSize);
	void Smoothing_Angran(int nSize, double stepSize);


	int SmoothingPoint(int pid, double stepSize);
	int LaplaceSmoothingPoint(int pid, double stepSize);
	int SmoothingPointBoundary(int pid, double stepSize);
	int SmoothingPointBoundarySharp(int pid, double stepSize);
	void GetHexVolAndCenter(int eid, double& vol, double center[3]);
	void GetQuadInfo(int fcid, int ploc, double& area, double center[3], double nm[3]);
	void GetHexMinJacob(int eid, double& minJacob, double GaussPos[3]);
	void GetGaussPoint(int ng, vector<double>& Gpt, vector<double>& wght);
	void JacobEval(int eid, double u, double v, double w, double& detJ);
	void JacobEval_Scale(int eid, double u, double v, double w, double& detJ);
	void GlobalMinJacob(double& minJacob_glb, int& min_pos, vector<int>& nBadEle);
	void Optimizing(int nStep, double stepSize);
	int OptimizingElement(int eid, double stepSize);
	int OptimizeElement(int eid, double stepSize);
	double AverageEdgeLength(int eid);
	void JacobEval_Grad(int eid, int ploc, double Jacob_Grad[3]);

	//void GlobalMinJacobBEXT(double& minJacob_glb, int& min_pos, vector<int>& nBadEle);
	//void GetHexMinJacobBEXT(int eid, double& minJacob_ele, double GaussPos[3]);
	//void JacobEval_ScaleBEXT(int eid, double u, double v, double w, double& detJ);


	void Smoothing_Adapt(int nSize, double stepSize);
	int SmoothingPoint_Adapt(int lev, int pid, double stepSize);
	int SmoothingPointBoundary_Adapt(int lev, int pid, double stepSize);
	void GetHexMinJacob_Adapt(int lev, int eid, double& minJacob_ele, double GaussPos[3]);
	void GetHexVolAndCenter_Adapt(int lev, int eid, double& vol, double center[3]);
	void GetQuadInfo_Adapt(int lev, int fcid, int ploc, double& area, double center[3], double nm[3]);
	//void JacobEval_Adapt(int eid, double u, double v, double w, double& detJ);
	void JacobEval_Scale_Adapt(int lev, int eid, double u, double v, double w, double detJ[2]);
	void GlobalMinJacob_Adapt(double& minJacob_glb, array<int, 2>& min_pos, vector<array<int, 2>>& BadEle);
	void BadElementFlag_Adapt(const vector<array<int,2>>& rfid);


	void Optimizing_1(int nStep, int pnitrMax, double stepSize);
	double MaxEdgeLength();
	double OptimizePoint(int pid, int nitrMax, double stepSize);
	void TranslateScale(int pid, double trs[3], vector<double>& scl);
	void TranslateScale_Reverse(int pid, double trs[3], const vector<double>& scl);
	void GetAdvanceDirection(int pid, double grad[3], double dir[3], double& f);
	double GetStepSize(int pid, double grad[3], double dir[3], double keta);
	void objGrad(int eid, int iloc, double delta, double& obj, double grad[3], double grad2[3][3]);
	void GetJacobMat(int eid, int iloc, double Jmat[3][3], double DJmat[3][3][3], double& detJ, double grad[3]);

	void Optimizing_glb(int nStep, double stepSize);

	//quality Laplace
	void LaplaceSmoothing(int nstep = 100);
	void LaplaceSmooth_Interior(int pid);
	void LaplaceSmooth_Boundary_NonSharp(int pid);
	void LaplaceSmooth_Boundary_Sharp(int pid);
	void SetSharpFeature_Manual(string fn);

	bool ComputeElementNormal_Quad(int cnct[4], array<double, 3>& normal_quad);
	bool CrossProduct(double vector_1[3], double vector_2[3], array<double, 3>& result);
	bool DotProduct(vector<array<double, 3>> check_normal_surface_from_edge, double & angle);
	void SetSharpFeature_Manual_Yu(string fn);

	void QualityImprove();
	void Preprocess();
	void BezierExtract_UBsplines(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);


	//hierarchical splines, identify
	void Identify_Refine_Test(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);
	void Identify_BadElement(vector<array<int, 2>>& rfid, vector<array<int, 2>>& gst);
	void StrongDirichletBC_Elasticity_PID_Adapt(string fn, vector<int>& IDBC, vector<double>& gh);


	//adaptive using open knot vectors
	void Initialize_THS3D(string fn);
	void SetZeroKnotInterval();
	void BuildSplines_Unstruct_Open();
	void BuildElementSplines_Interior_Open(int eid);
	void BuildElementSplines_Boundary_Open(int eid);


	//C0 Bezier + C1 B-splines + C2 B-splines
	void Run_C0C1Bezier(string fn, int nrf, vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void InitializeC0Component();
	void AddC0C1Bezier();
	void AddC0C1Bezier_Boundary();
	void TruncateC0C1();
	void BezierExtractC0C1_Trun(int eid, vector<vector<double>>& cmat);//for interior element
	void OutputGeom_C0C1(string fn);
	void GeomMap_C0C1(int eid, double u[3], double pt[3]);
	//void BezierExtractC1(int eid, vector<>);
	void Global_Subdivide_C0C1Bezier();
	void BezierExtract_C0C1Bezier(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void EnlargeRefineRegion();

	//C0 Bezier for comparison
	void Run_C0Bezier(string fn, int nrf, vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void AddC0Bezier();
	void OutputGeom_C0(string fn);
	void GeomMap_C0(int eid, double u[3], double pt[3]);
	//void Global_Subdivide_C0Bezier();
	void BezierExtract_C0Bezier(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void BezierExtract_C0Bezier_1(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);

	//c0 Bezier, improved version
	void Run_C02iBezier(string fn, int nrf, vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void InitializeC0Component_C02i();
	void AddC0C1Bezier_C02i();
	void TruncateC02i();

	//C012 construction with improved boundary
	//void Run_C012ib(string fn, int nrf, vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	//void InitializeC0Component_C012ib();

	//strongly imposed Dirichlet boundary condition
	void StrongDirichletBC(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void StrongDirichletBC_C0(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void StrongDirichletBC_Poisson_FromFile(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void WeakDirichletBC_Poisson_FromFile(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void StrongDirichletBC_Elasticity(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void StrongDirichletBC_Elasticity_FromFile(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void StrongDirichletBC_Elasticity_Manual(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void StrongDirichletBC_Elasticity_PID(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void StrongDirichletBC_Elasticity_EID(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void StrongDirichletBC_Elasticity_Coor(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void ReadBC(string fn, int nbc, vector<vector<int>>& ibc, vector<array<double,3>>& bc_disp);
	void ReadBC_PID(string fn, int nbc, vector<vector<int>>& pid, vector<array<int, 3>>& bc_dof, vector<array<double, 3>>& bc_disp);

	//patch test
	void PatchTest_BC(vector<int>& IDBC, vector<double>& gh);
	void PatchTest_OneElement(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);

	//LSDYNA file format
	void AngranOutputMesh(const vector<BezierElement3D>& bzmesh, string fn);

	void ReadBEXT3D(string fn, vector<BezierElement3D>& bzmesh, vector<int> &IDBC, vector<double> &gh);
	void WriteBezierInfo_Blend_LSDYNA(string fn, vector<BezierElement3D>& bzmesh);
	//void FindEffectiveDOF(vector<BezierElement3D>& bzmesh, vector<int>& aloc, vector<vector<int>>& eIEN);
	void WriteBezierInfo_AllSpline_LSDYNA(string fn, vector<BezierElement3D>& bzmesh);
	void WriteBezierInfo_AllSpline_LSDYNA_LocalRefine(string fn, vector<BezierElement3D>& bzmesh);
	void WriteBezier(const vector<BezierElement3D>& bzmesh, string fn1);


	//ABAQUS with GEM, plate with a hole
	void GeneratePlateHole3D(string fn);
	void ReadInputQuadVTK(string fn, vector<array<double,3>>& pts, vector<array<int,4>>& cnct);
	void OutputQuadVTK(string fn, vector<array<double, 3>>& pts, vector<array<int, 4>>& cnct);
	void Rescale(vector<array<double, 3>>& pts);
	void KnotVectorWeights(vector<vector<double>>& kv, vector<vector<double>>& wt);
	//void AdditionInput2D(vector<double>& ku, vector<double>& kv, vector<double>& wu, vector<double>& wv);
	void Mirror2D(int mr_flag, const vector<array<double, 3>>& pts0, vector<array<double, 3>>& pts);
	void ReorderPoints(int flag, vector<array<double, 3>>& pts);
	//void Rotate2D(int rt_flag, const vector<array<double, 3>>& pts0, vector<array<double, 3>>& pts, vector<int>& intfc);
	void Sweep(int s_flag, const vector<array<double, 3>>& pts2d, const vector<array<int, 4>>& cnct2d);
	void GetPatchConnect(const vector<vector<int>>& pid2d, vector<vector<int>>& pid3d);
	void NURBS_BezierExtraction(const vector<vector<int>>& pid, vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);
	void ReorderBezier(const vector<BezierElement3D>& bzmesh, vector<vector<int>>& bzabq);
	void WriteBezier_Abaqus(const vector<BezierElement3D>& bzmesh, string fn);
	void WriteBezier_Abaqus_LocalRefine(const vector<BezierElement3D>& bzmesh, string fn);
	void WriteBezier_h(const vector<BezierElement3D>& bzmesh, string fn);
	//void RefineBsplinePatch(int);
	void SetBSplinePatch(const vector<vector<int>>& pid);
	void SetBSPatchPointID();
	void BSPatch_BezierExtraction(vector<BezierElement3D>& bzmesh, vector<int>& IDBC, vector<double>& gh);

	void TTSP3D_Clear();

	void ReadInputTriMesh(string fn);
	void ProjectBoundary();
	void ProjectBoundary_NonSharp();
	void ProjectBoundary_Sharp();

	//used for polycube
	void SplitPolyCube(string fn);
	void IdentifyHex_Honda2(vector<int>& elist, vector<array<int,3>>& stype);
	void IdentifyHex_NAVAIR(vector<int>& elist, vector<array<int, 3>>& stype);
	void SplitHexOneDir(int eid, int dir, int np);
	void SplitHexTwoDir(int eid, int dir, int np);
	void SplitHexThrDir(int eid, int np);
	void RemoveRepeatPoint();
	void OutputCornerPoint(string fn);


	//vector<int> paid;//active control points
	vector<int> haid;//active hex elements
	vector<int> faid;//active faces
	vector<int> eaid;//active edges

private:
	vector<Vertex3D> cp;//control points
	vector<Vertex3D> cp_smooth; // for smoothing, store control points after smoothing
	vector<array<double, 3>> cpa;
	vector<Element3D> tmesh;//elements in T-mesh
	vector<Edge3D> tmedge;
	vector<Face3D> tmface;
	unsigned int npt_old;
	unsigned int nel_old;
	unsigned int nfc_old;
	unsigned int ned_old;
	vector<vector<Vertex3D>> hcp;
	vector<vector<Element3D>> hmesh;
	vector<vector<Edge3D>> hedge;
	vector<vector<Face3D>> hface;
	vector<vector<double>> kvec;

	vector<BPatch3D> bsp;

	vector<Element3D> tmtmp;//elements in T-mesh
	vector<Edge3D> tmedtmp;
	vector<Face3D> tmfctmp;

	vector<MFNode> mp;//meshfree nodes
	vector<array<double, 3>> bzcp;//Bezier control points
	//vector<array<double, 3>> bzcp_c01;//Bezier control points
	vector<int> bzcp_c01;

	vector<array<double, 3>> ptri;
	vector<array<int, 3>> etri;

	//used for solution
	double dmrg[3][2];
	double nmpl[3];
	double acoef;
	double dmlen[3];
};

void ManualPolycube_Honda2(string fn);

void PolycubeDeleteRepeatPoint(string fn);

void ManualAddCube(string fn);

void PolyCubeCorrespond(string fn_poly, string crsp);

void DeleteVoidElement(string fn);

void RepairConnectivity(string fn);


void ManualPolycube_NAVAIR(string fn_in, string fn_out);

void PolyCubeSubdiv(string fn_in);

#endif