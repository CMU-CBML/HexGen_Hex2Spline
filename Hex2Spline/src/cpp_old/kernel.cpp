#include "kernel.h"

void kernel::run()
{
	int niter(20);
	unsigned int i;
	TruncatedTspline_3D tt3;
	//tt3.SetProblem("../io/hex_input/cube_dense");
	//tt3.SetProblem("../io/hex_input/cube4");
	tt3.SetProblem("../io/hex_input/cube_coarse_0");
	string fld("../io/hex_out_2/"), fn("test16_");
	vector<int> dof_list(niter,0);
	vector<double> err_list(niter,0.);
	for (int itr = 0; itr < niter; itr++)
	{
		vector<BezierElement3D> bzmesh;
		vector<int> IDBC;
		vector<double> gh, err;
		//cout << "interface...\n";
		tt3.AnalysisInterface_Poisson(bzmesh, IDBC, gh);
		//tt3.AnalysisInterface_Laplace(bzmesh, IDBC, gh);
		//cout << "interface done!\n";
		//getchar();

		//cout << "npt: " << IDBC.size() << "\n";
		//getchar();

		Laplace lap;
		lap.SetProblem(IDBC, gh);
		stringstream ss;
		ss << itr;
		//cout << "before simulation...\n";
		//getchar();
		lap.Run(bzmesh, fld+fn+ss.str(), err);
		//cout << "simulation done!\n";
		//getchar();

		//tt3.VisualizeBezier(bzmesh, fld + fn + ss.str());

		double errL2(0.);
		for (i = 0; i < err.size(); i++) errL2 += err[i];
		errL2 = sqrt(errL2);
		dof_list[itr] = IDBC.size();
		err_list[itr] = errL2;
		//cout << "DOF: " << IDBC.size() << "\n";
		//cout << "L2-norm error: " << errL2 << "\n";
		//getchar();

		if (itr < niter - 1)
		{
			cout << itr << " refining...\n";
			//distribute error
			vector<array<double, 2>> eh(bzmesh.size());
			for (i = 0; i < bzmesh.size(); i++)
			{
				eh[i][0] = bzmesh[i].prt[0]; eh[i][1] = bzmesh[i].prt[1];
			}
			vector<array<int, 2>> rfid, gst;
			//tt3.Identify_Poisson(eh, err, rfid, gst);
			tt3.Identify_Poisson_1(eh, err, rfid, gst);
			//tt3.Identify_Laplace(eh, err, rfid, gst);
			tt3.Refine(rfid, gst);
			cout << "Refining done\n";
			tt3.OutputGeom_All(fld+fn + ss.str()+"_geom");
			//cout << "Output Geom done!\n";
			//getchar();
		}
	}

	//output error
	output_err(fld + fn + "err", dof_list, err_list);
}

void kernel::run_complex()
{
	int niter(5);
	unsigned int i;
	TruncatedTspline_3D tt3;
	//tt3.SetProblem("../io/hex_input/cube_dense");
	//tt3.SetProblem("../io/hex_input/cube9");
	//tt3.SetProblem("../io/hex_input/cube_coarse_0");
	//tt3.SetProblem("../io/hex_input/rod");
	//tt3.SetProblem("../io/hex_input/cad1");
	//tt3.SetProblem("../io/hex_input/crosshole1");
	//tt3.SetProblem("../io/hex_input/varco");
	tt3.SetProblem("../io/hex_input/statue");
	//tt3.SetProblem("../io/hex_input/hook");
	//tt3.SetProblem("../io/hex_input/gear");
	string fld("../io/complex2/");
	//string fld("../io/benchmark3/");
	//string fn("cube_THS1_");
	//string fn("cube9_");
	//string fn("rod1_");
	//string fn("cad1_");
	//string fn("crosshole2_");
	//string fn("varco1_");
	string fn("statue1_");
	//string fn("hook1_");
	//string fn("gear1_");
	//tt3.InputCheck(fld+fn);
	//tt3.MeshRepair("../io/hex_input/crosshole", "../io/hex_input/crosshole1");

	double remove[3][2];
	tt3.GetRemoveRegion(remove);
	//tt3.OutputRemoveCM(fld+fn,remove);
	//cout << "done\n";
	//getchar();

	vector<int> dof_list(niter, 0);
	vector<double> err_list(niter, 0.);
	vector<array<int, 2>> ebc, pbc;
	vector<double> edisp, pdisp;
	tt3.SetInitialBC(ebc, edisp, pbc, pdisp);
	for (int itr = 0; itr < niter; itr++)
	{
		vector<BezierElement3D> bzmesh;
		vector<int> IDBC;
		vector<double> gh, err;
		tt3.SetBC(ebc, edisp, pbc, pdisp);
		//cout << pbc.size() << "\n";
		//cout << "interface...\n";
		//tt3.AnalysisInterface_Poisson(bzmesh, IDBC, gh);
		//tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);
		tt3.AnalysisInterface_Laplace(pbc,pdisp,bzmesh, IDBC, gh);
		//cout << "interface done!\n";
		//getchar();

		//cout << "npt: " << IDBC.size() << "\n";
		//getchar();

		Laplace lap;
		lap.SetProblem(IDBC, gh);
		stringstream ss;
		ss << itr;
		//cout << "before simulation...\n";
		//getchar();
		lap.GetRemoveRegion(remove);
		lap.Run(bzmesh, fld + fn + ss.str(), err);
		//cout << "simulation done!\n";
		//getchar();

		//tt3.VisualizeBezier(bzmesh, fld + fn + ss.str());

		double errL2(0.);
		for (i = 0; i < err.size(); i++) errL2 += err[i];
		errL2 = sqrt(errL2);
		dof_list[itr] = IDBC.size();
		err_list[itr] = errL2;
		//cout << "DOF: " << IDBC.size() << "\n";
		//cout << "L2-norm error: " << errL2 << "\n";
		//getchar();

		output_err(fld + fn+ss.str() + "_err", dof_list, err_list);

		if (itr < niter - 1)
		{
			cout << itr << " refining...\n";
			//distribute error
			vector<array<double, 2>> eh(bzmesh.size());
			for (i = 0; i < bzmesh.size(); i++)
			{
				eh[i][0] = bzmesh[i].prt[0]; eh[i][1] = bzmesh[i].prt[1];
			}
			vector<array<int, 2>> rfid, gst;
			//tt3.Identify_Poisson_1(eh, err, rfid, gst);
			tt3.Identify_Laplace(eh, err, rfid, gst);
			//tt3.OutputRefineID(fld + fn + ss.str(), rfid, gst);
			//tt3.InputRefineID("../io/benchmark1/cube1_" + ss.str(), rfid, gst);
			tt3.Refine(rfid, gst);
			cout << "Refining done\n";
			//tt3.OutputGeom_All(fld + fn + ss.str() + "_geom");
			//cout << "Output Geom done!\n";
			//getchar();
		}
	}

	//output error
	//output_err(fld + fn + "err", dof_list, err_list);
}

void kernel::run_complex_fit()
{
	int niter(20);
	double thresh(0.008);//cube
	//double thresh(0.065);//rod
	//double thresh(0.0102);//base
	//double thresh(0.0704);//head
	unsigned int i;
	double xy[3][2], nm[3], a(50.);
	TruncatedTspline_3D tt3;

	//tt3.SetProblem("../io/hex_input/cube_dense");
	//tt3.SetProblem("../io/hex_input/cube4");
	//tt3.SetProblem("../io/hex_input/cube_coarse_0");
	//tt3.SetProblem("../io/hex_input/rod");
	//tt3.SetProblem("../io/hex_input/cad1");
	//tt3.SetProblem("../io/hex_input/crosshole1");
	//tt3.SetProblem("../io/hex_input/varco");
	//tt3.SetProblem("../io/hex_input/statue");
	//tt3.SetProblem("../io/hex_input/hook");
	//tt3.SetProblem("../io/hex_input/gear");
	//tt3.SetProblem("../io/yuxuan/rockerarm");
	tt3.SetProblem("../io/yuxuan/fertility");

	//string fld("../io/local2/");
	//string fld("../io/benchmark3/");
	//string fld("../io/yuxuan/rockerarm/");
	string fld("../io/yuxuan/fertility/");

	//string fn("cube_THS1_");
	//string fn("cube9_");
	//string fn("cube_");
	//string fn("cube_fit3_");
	//string fn("rod_");
	//string fn("base1_");
	//string fn("crosshole2_");
	//string fn("varco1_");
	//string fn("head_");
	//string fn("hook_");
	//string fn("gear1_");
	//string fn("rockerarm_");
	string fn("fertility_");

	//tt3.InputCheck(fld+fn);
	//tt3.MeshRepair("../io/hex_input/crosshole", "../io/hex_input/crosshole1");

	//double remove[3][2];
	//tt3.GetRemoveRegion(remove);
	tt3.SetDomainRange(xy, nm, a);

	//tt3.OutputRemoveCM(fld+fn,remove);
	//cout << "done\n";
	//getchar();

	vector<int> dof_list(niter, 0);
	//vector<double> h_max(niter, 0.);
	vector<double> err_list(niter, 0.);
	//vector<array<int, 2>> ebc, pbc;
	//vector<double> edisp, pdisp;
	//tt3.SetInitialBC(ebc, edisp, pbc, pdisp);

	double errL2(1.e6);
	for (int itr = 0; itr < niter; itr++)
	//int itr(0);
	//while (errL2 > thresh && itr < niter)
	{
		vector<BezierElement3D> bzmesh;
		vector<int> IDBC;
		vector<double> gh, err;
		//tt3.SetBC(ebc, edisp, pbc, pdisp);
		//cout << pbc.size() << "\n";
		//cout << "interface...\n";
		//tt3.AnalysisInterface_Poisson(bzmesh, IDBC, gh);
		tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);

		//if (itr == 1)
		//{
		//	tt3.WriteBezier_h(bzmesh, fld + "rockerarm_1");
		//	break;
		//}

		Laplace lap;
		lap.SetProblem(IDBC, gh);
		stringstream ss;
		ss << itr;
		//cout << "before simulation...\n";
		//getchar();
		//lap.GetRemoveRegion(remove);
		lap.GetEqParameter(xy, nm, a);
		lap.Run(bzmesh, fld + fn + ss.str(), err);
		//cout << "simulation done!\n";
		//getchar();

		//tt3.VisualizeBezier(bzmesh, fld + fn + ss.str());

		//double errL2(0.);
		errL2 = 0.;
		for (i = 0; i < err.size(); i++) errL2 += err[i];
		errL2 = sqrt(errL2);
		dof_list[itr] = IDBC.size();
		//h_max[itr] = tt3.MaxElementSize(bzmesh);
		err_list[itr] = errL2;
		//cout << "DOF: " << IDBC.size() << "\n";
		//cout << "L2-norm error: " << errL2 << "\n";
		//getchar();

		output_err(fld + fn + ss.str() + "_err", dof_list, err_list);
		//output_err(fld + fn + ss.str() + "_err_h", h_max, err_list);

		if (itr < niter - 1)
		//if (errL2 > thresh)
		{
			cout << itr << " refining...\n";
			//distribute error
			vector<array<double, 2>> eh(bzmesh.size());
			for (i = 0; i < bzmesh.size(); i++)
			{
				eh[i][0] = bzmesh[i].prt[0]; eh[i][1] = bzmesh[i].prt[1];
			}
			
			vector<array<int, 2>> rfid, gst;
			tt3.Identify_Poisson_1(eh, err, rfid, gst);
			//tt3.Identify_Laplace(eh, err, rfid, gst);
			//tt3.OutputRefineID(fld + fn + ss.str(), rfid, gst);
			//tt3.InputRefineID("../io/benchmark1/cube1_" + ss.str(), rfid, gst);
			tt3.Refine(rfid, gst);
			cout << "Refining done\n";

			//tt3.OutputGeom_All(fld + fn + ss.str() + "_geom");
			//cout << "Output Geom done!\n";
			//getchar();
		}
		
		//itr++;
	}

	//output error
	//output_err(fld + fn + "err", dof_list, err_list);
}

void kernel::run_complex_glb()
{
	int niter(1);
	unsigned int i;
	double xy[3][2], nm[3], a(50.);
	TruncatedTspline_3D tt3;
	tt3.InitializeMesh("../io/hex_input/cube4");
	//tt3.InitializeMesh("../io/hex_input/rod");
	//tt3.InitializeMesh("../io/hex_input/hook");
	//tt3.InitializeMesh("../io/hex_input/cad1");
	//tt3.InitializeMesh("../io/hex_input/statue");
	//tt3.InitializeMesh("../io/hex_input/gear");

	//tt3.SetProblem("../io/hex_input/cube_dense");
	//tt3.SetProblem("../io/hex_input/cube3");
	//tt3.SetProblem("../io/hex_input/gear");

	string fld("../io/global1/");
	//string fld("../io/benchmark3/");

	string fn("cube4_2");
	//string fn("rod1_2");
	//string fn("hook1_");
	//string fn("base0_");
	//string fn("head1_");
	//string fn("gear1_");

	tt3.run_Global(2, fld + fn);
	cout << "Done global refinement!\n";
	//getchar();

	//double remove[3][2];
	//tt3.GetRemoveRegion(remove);
	tt3.SetDomainRange(xy, nm, a);

	//tt3.OutputRemoveCM(fld+fn,remove);
	//cout << "done\n";
	//getchar();

	vector<int> dof_list(niter, 0);
	vector<double> h_max(niter, 0.);
	vector<double> err_list(niter, 0.);
	for (int itr = 0; itr < niter; itr++)
	{
		vector<BezierElement3D> bzmesh;
		vector<int> IDBC;
		vector<double> gh, err;
		tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);

		Laplace lap;
		lap.SetProblem(IDBC, gh);
		stringstream ss;
		ss << itr;
		//lap.GetRemoveRegion(remove);
		lap.GetEqParameter(xy, nm, a);
		lap.Run(bzmesh, fld + fn, err);

		//tt3.VisualizeBezier(bzmesh, fld + fn + ss.str());
		//cout << "Done visualize Bezier!\n";
		//getchar();

		double errL2(0.);
		for (i = 0; i < err.size(); i++) errL2 += err[i];
		errL2 = sqrt(errL2);
		dof_list[itr] = IDBC.size();
		h_max[itr] = tt3.MaxElementSize(bzmesh);
		err_list[itr] = errL2;
		//cout << "DOF: " << IDBC.size() << "\n";
		//cout << "L2-norm error: " << errL2 << "\n";
		//getchar();

		output_err(fld + fn + "_err", dof_list, err_list);
		output_err(fld + fn + "_err_h", h_max, err_list);
	}
}

void kernel::run_Bspline()
{
	int niter(1);
	unsigned int i;
	//string fld("../io/BSP/");
	string fld("../io/yuxuan/cube/");
	//string fn("ZBSP3_");
	string fn("BSP2");
	TruncatedTspline_3D tt3;
	vector<BezierElement3D> bzmesh;
	vector<int> IDBC;
	vector<double> gh, err;
	int nsp(1);
	int nex[3] = { nsp, nsp, nsp };
	//tt3.CreateBsplines(bzmesh, IDBC, gh);
	tt3.CreateBsplines(nex);
	int nref(1);
	for (int i = 0; i < nref; i++)
	{
		tt3.Bsplines_Refine();
	}
	tt3.Bspline_BezierExtract(bzmesh, IDBC, gh);

	tt3.WriteBezier_Abaqus(bzmesh, fld + fn);

	//LinearElasticity le;
	//le.Set_Problem_Bezier(IDBC, gh);
	//le.Run_Coupling_Bezier(bzmesh, fld + fn);

	//tt3.VisualizeBezier(bzmesh, fld + fn);
	//cout << "Done visualize Bezier!\n";
	//getchar();

	//Laplace lap;
	//lap.SetProblem(IDBC, gh);
	//lap.Run(bzmesh, fld + fn, err);

	//vector<int> dof_list(niter, 0);
	//vector<double> h_max(niter, 0.);
	//vector<double> err_list(niter, 0.);
	//double errL2(0.);
	//for (i = 0; i < err.size(); i++) errL2 += err[i];
	//errL2 = sqrt(errL2);
	//dof_list[0] = IDBC.size();
	//h_max[0] = tt3.MaxElementSize(bzmesh);
	//err_list[0] = errL2;

	//output_err(fld + fn + "_err", dof_list, err_list);
	//output_err(fld + fn + "_err_h", h_max, err_list);
}

void kernel::run_Bspline_fit()
{
	int niter(1);
	unsigned int i;
	double xy[3][2], nm[3], a(50.);
	string fld("../io/BSP/");
	string fn("BSP2_fit4");
	TruncatedTspline_3D tt3;
	vector<BezierElement3D> bzmesh;
	vector<int> IDBC;
	vector<double> gh, err;
	int nsp(2);
	int nex[3] = { nsp, nsp, nsp };
	//tt3.CreateBsplines(bzmesh, IDBC, gh);
	tt3.CreateBsplines(nex);
	int nref(4);
	for (int i = 0; i < nref; i++)
	{
		tt3.Bsplines_Refine();
	}
	tt3.SetDomainRange_Bsplines(xy, nm, a);//for solution 7
	tt3.Bspline_BezierExtract_fit(bzmesh, IDBC, gh);
	tt3.FittingBC_Bsplines(bzmesh, IDBC, gh);

	//tt3.VisualizeBezier(bzmesh, fld + fn);
	//cout << "Done visualize Bezier!\n";
	//getchar();

	Laplace lap;
	lap.SetProblem(IDBC, gh);
	lap.GetEqParameter(xy, nm, a);//for solution 7
	lap.Run(bzmesh, fld + fn, err);

	vector<int> dof_list(niter, 0);
	vector<double> h_max(niter, 0.);
	vector<double> err_list(niter, 0.);
	double errL2(0.);
	for (i = 0; i < err.size(); i++) errL2 += err[i];
	errL2 = sqrt(errL2);
	dof_list[0] = IDBC.size();
	h_max[0] = tt3.MaxElementSize(bzmesh);
	err_list[0] = errL2;

	output_err(fld + fn + "_err", dof_list, err_list);
	output_err(fld + fn + "_err_h", h_max, err_list);
}

void kernel::run_benchmark()
{
	int niter(4);
	unsigned int i;
	TruncatedTspline_3D tt3;
	//tt3.SetProblem("../io/hex_input/cube_dense");
	//tt3.SetProblem("../io/hex_input/cube9");
	//tt3.SetProblem("../io/hex_input/cube_coarse_0");
	tt3.SetProblem("../io/hex_input/cube_coarse_2");
	string fld("../io/benchmark4/");
	string fn("cube_THS1_");
	//string fn("cube_HS1_");

	//double remove[3][2];
	//tt3.GetRemoveRegion(remove);
	//tt3.OutputRemoveCM(fld+fn,remove);
	//cout << "done\n";
	//getchar();

	vector<int> dof_list(niter, 0);
	vector<double> err_list(niter, 0.);
	vector<array<int, 2>> ebc, pbc;
	vector<double> edisp, pdisp;
	//tt3.SetInitialBC(ebc, edisp, pbc, pdisp);
	for (int itr = 0; itr < niter; itr++)
	{
		vector<BezierElement3D> bzmesh;
		vector<int> IDBC;
		vector<double> gh, err;
		//tt3.SetBC(ebc, edisp, pbc, pdisp);
		//cout << pbc.size() << "\n";
		//cout << "interface...\n";
		tt3.AnalysisInterface_Poisson(bzmesh, IDBC, gh);
		//tt3.AnalysisInterface_Laplace(pbc,pdisp,bzmesh, IDBC, gh);
		//cout << "interface done!\n";
		//getchar();

		//cout << "npt: " << IDBC.size() << "\n";
		//getchar();

		Laplace lap;
		lap.SetProblem(IDBC, gh);
		stringstream ss;
		ss << itr;
		//cout << "before simulation...\n";
		//getchar();
		//lap.GetRemoveRegion(remove);
		lap.Run(bzmesh, fld + fn + ss.str(), err);
		//cout << "simulation done!\n";
		//getchar();

		//tt3.VisualizeBezier(bzmesh, fld + fn + ss.str());

		double errL2(0.);
		for (i = 0; i < err.size(); i++) errL2 += err[i];
		errL2 = sqrt(errL2);
		dof_list[itr] = IDBC.size();
		err_list[itr] = errL2;
		//cout << "DOF: " << IDBC.size() << "\n";
		//cout << "L2-norm error: " << errL2 << "\n";
		//getchar();

		output_err(fld + fn + ss.str() + "_err", dof_list, err_list);

		if (itr < niter - 1)
		{
			cout << itr << " refining...\n";
			//distribute error
			vector<array<double, 2>> eh(bzmesh.size());
			for (i = 0; i < bzmesh.size(); i++)
			{
				eh[i][0] = bzmesh[i].prt[0]; eh[i][1] = bzmesh[i].prt[1];
			}
			vector<array<int, 2>> rfid, gst;
			//tt3.Identify_Poisson_1(eh, err, rfid, gst);
			//tt3.Identify_Laplace(eh, err, rfid, gst);
			//tt3.OutputRefineID(fld + fn + ss.str(), rfid, gst);
			tt3.InputRefineID("../io/benchmark3/cube_THS1_" + ss.str(), rfid, gst);
			tt3.Refine(rfid, gst);
			cout << "Refining done\n";
			//tt3.OutputGeom_All(fld + fn + ss.str() + "_geom");
			//cout << "Output Geom done!\n";
			//getchar();
		}
	}

	//output error
	//output_err(fld + fn + "err", dof_list, err_list);
}

void kernel::run_leastsquare()
{
	int niter(2);
	unsigned int i;
	TruncatedTspline_3D tt3;
	//tt3.SetProblem("../io/hex_input/cube_dense");
	//tt3.SetProblem("../io/hex_input/cube4");
	tt3.SetProblem("../io/hex_input/cube_coarse_2");
	string fld("../io/least_square1/");
	//string fn("cube_THS2_2");
	string fn("cube_HS2_2");

	//tt3.GlobalRefine(2);

	for (int itr = 0; itr < niter; itr++)
	{
		cout << itr << " refining...\n";
		vector<array<int, 2>> rfid, gst;
		//tt3.Identify_Poisson_1(eh, err, rfid, gst);
		//tt3.Identify_Laplace(eh, err, rfid, gst);
		//tt3.Identify_LeastSquare(rfid, gst);
		tt3.Identify_LeastSquare_Line(rfid, gst);
		//tt3.OutputRefineID(fld + fn + ss.str(), rfid, gst);
		//tt3.InputRefineID("../io/benchmark3/cube_THS1_" + ss.str(), rfid, gst);
		tt3.Refine(rfid, gst);
		stringstream ss;
		ss << itr;
		cout << "Refining done\n";
		//tt3.OutputGeom_All(fld + fn + ss.str());
		//cout << "Output Geom done!\n";
		//getchar();
	}

	//cout << "Refining done\n";
	//tt3.OutputGeom_All(fld + fn + "_geom");
	//cout << "geom done!\n";
	//getchar();

	vector<BezierElement3D> bzmesh;
	vector<int> IDBC;
	vector<double> gh, err;
	//tt3.AnalysisInterface_LeastSquare(bzmesh, IDBC, gh);
	tt3.AnalysisInterface_Poisson(bzmesh, IDBC, gh);
	LeastSquare lap;
	lap.SetProblem(IDBC, gh);
	lap.Run(bzmesh, fld + fn, err);

	//output error
	//output_err(fld + fn + "err", dof_list, err_list);
}

void kernel::run_leastsquare_1()
{
	int niter(4);
	unsigned int i;
	TruncatedTspline_3D tt3;
	//tt3.SetProblem("../io/hex_input/cube_dense");
	//tt3.SetProblem("../io/hex_input/cube9");
	//tt3.SetProblem("../io/hex_input/cube_coarse_0");
	tt3.SetProblem("../io/hex_input/cube_coarse_2");
	string fld("../io/least_square1/");
	string fn("cube_THS1_");
	//string fn("cube_HS1_");

	//double remove[3][2];
	//tt3.GetRemoveRegion(remove);
	//tt3.OutputRemoveCM(fld+fn,remove);
	//cout << "done\n";
	//getchar();

	vector<int> dof_list(niter, 0);
	vector<double> err_list(niter, 0.);
	vector<array<int, 2>> ebc, pbc;
	vector<double> edisp, pdisp;
	//tt3.SetInitialBC(ebc, edisp, pbc, pdisp);
	for (int itr = 0; itr < niter; itr++)
	{
		vector<BezierElement3D> bzmesh;
		vector<int> IDBC;
		vector<double> gh, err;
		//tt3.SetBC(ebc, edisp, pbc, pdisp);
		//cout << pbc.size() << "\n";
		//cout << "interface...\n";
		//tt3.AnalysisInterface_Poisson(bzmesh, IDBC, gh);
		//tt3.AnalysisInterface_Laplace(pbc,pdisp,bzmesh, IDBC, gh);
		tt3.AnalysisInterface_LeastSquare(bzmesh, IDBC, gh);
		//cout << "interface done!\n";
		//getchar();

		//cout << "npt: " << IDBC.size() << "\n";
		//getchar();

		Laplace lap;
		lap.SetProblem(IDBC, gh);
		stringstream ss;
		ss << itr;
		//cout << "before simulation...\n";
		//getchar();
		//lap.GetRemoveRegion(remove);
		lap.Run(bzmesh, fld + fn + ss.str(), err);
		//cout << "simulation done!\n";
		//getchar();

		//tt3.VisualizeBezier(bzmesh, fld + fn + ss.str());

		double errL2(0.);
		for (i = 0; i < err.size(); i++) errL2 += err[i];
		errL2 = sqrt(errL2);
		dof_list[itr] = IDBC.size();
		err_list[itr] = errL2;
		//cout << "DOF: " << IDBC.size() << "\n";
		//cout << "L2-norm error: " << errL2 << "\n";
		//getchar();

		output_err(fld + fn + ss.str() + "_err", dof_list, err_list);

		if (itr < niter - 1)
		{
			cout << itr << " refining...\n";
			//distribute error
			vector<array<double, 2>> eh(bzmesh.size());
			for (i = 0; i < bzmesh.size(); i++)
			{
				eh[i][0] = bzmesh[i].prt[0]; eh[i][1] = bzmesh[i].prt[1];
			}
			vector<array<int, 2>> rfid, gst;
			//tt3.Identify_Poisson_1(eh, err, rfid, gst);
			//tt3.Identify_Laplace(eh, err, rfid, gst);
			//tt3.OutputRefineID(fld + fn + ss.str(), rfid, gst);
			tt3.InputRefineID("../io/benchmark3/cube_THS1_" + ss.str(), rfid, gst);
			tt3.Refine(rfid, gst);
			cout << "Refining done\n";
			//tt3.OutputGeom_All(fld + fn + ss.str() + "_geom");
			//cout << "Output Geom done!\n";
			//getchar();
		}
	}

	//output error
	//output_err(fld + fn + "err", dof_list, err_list);
}

void kernel::run_pipeline()
{
	int niter(5);
	unsigned int i;
	TruncatedTspline_3D tt3;
	tt3.SetProblem("../io/pipeline/base_disp");
	string fld("../io/pipeline/");
	string fn("base1_");

	//double remove[3][2];
	//tt3.GetRemoveRegion(remove);
	//tt3.OutputRemoveCM(fld+fn,remove);
	//cout << "done\n";
	//getchar();

	vector<BezierElement3D> bzmesh;
	tt3.PipelineBezierExtract(bzmesh);

	Laplace lap;
	lap.PipelineTmp(bzmesh,fld+fn);
}

void kernel::output_err(string fn, const vector<int>& dof, const vector<double>& err)
{
	string fname = fn + ".txt";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		for (i = 0; i < dof.size(); i++)
		{
			fout << dof[i] << " " << err[i] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void kernel::output_err(string fn, const vector<double>& dof, const vector<double>& err)
{
	string fname = fn + ".txt";
	ofstream fout;
	fout.open(fname.c_str());
	unsigned int i;
	if (fout.is_open())
	{
		for (i = 0; i < dof.size(); i++)
		{
			fout << dof[i] << " " << err[i] << "\n";
		}
		fout.close();
	}
	else
	{
		cout << "Cannot open " << fname << "!\n";
	}
}

void kernel::run_conv()
{
	int niter(1);
	unsigned int i;
	double xy[3][2], nm[3], a(50.);
	TruncatedTspline_3D tt3;
	//tt3.InitializeMesh("../io/hex_input/cube4");//this uses generalized Bezier elements
	//tt3.InitializeMesh("../io/hex_input/rod");
	//tt3.InitializeMesh("../io/hex_input/hook");
	//tt3.InitializeMesh("../io/hex_input/cad1");
	//tt3.InitializeMesh("../io/hex_input/statue");
	//tt3.InitializeMesh("../io/hex_input/gear");

	//tt3.SetProblem("../io/hex_input/cube_dense");
	//tt3.SetProblem("../io/hex_input/cube4");
	//tt3.SetProblem("../io/hex_input/gear");

	//tt3.InitializeWithTsplines("../io/hex_input/cube4");//this one uses T-spline definition
	//tt3.InitializeWithTsplines("../io/hex_input/cube_extraordinary");

	tt3.ReadHexMeshVTK("../io/hex_input/cubeex");
	//tt3.DeleteDuplicatePoint("../io/hex_input/cubeex");
	//cout << "done delete!\n"; getchar();
	//tt3.ReadHexMeshVTK("../io/hex_input/cube4");
	tt3.BuildSolidTspline();

	string fld("../io/global1/");
	//string fld("../io/benchmark3/");

	string fn("cube_ex_test1_0");
	//string fn("rod1_2");
	//string fn("hook1_");
	//string fn("base0_");
	//string fn("head1_");
	//string fn("gear1_");

	tt3.VisualizeSolidVTK(fld + fn);
	cout << "done solid vtk!\n"; getchar();

	//tt3.run_Global(2, fld + fn);
	//cout << "Done global refinement!\n";
	//getchar();

	//double remove[3][2];
	//tt3.GetRemoveRegion(remove);
	//tt3.OutputRemoveCM(fld+fn,remove);
	//cout << "done\n";
	//getchar();

	//tt3.SetDomainRange(xy, nm, a);

	//vector<int> dof_list(niter, 0);
	//vector<double> h_max(niter, 0.);
	//vector<double> err_list(niter, 0.);
	//for (int itr = 0; itr < niter; itr++)
	//{
	//	vector<BezierElement3D> bzmesh;
	//	vector<int> IDBC;
	//	vector<double> gh, err;
	//	tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);

	//	Laplace lap;
	//	lap.SetProblem(IDBC, gh);
	//	stringstream ss;
	//	ss << itr;
	//	//lap.GetRemoveRegion(remove);
	//	lap.GetEqParameter(xy, nm, a);
	//	lap.Run(bzmesh, fld + fn, err);

	//	//tt3.VisualizeBezier(bzmesh, fld + fn + ss.str());
	//	//cout << "Done visualize Bezier!\n";
	//	//getchar();

	//	double errL2(0.);
	//	for (i = 0; i < err.size(); i++) errL2 += err[i];
	//	errL2 = sqrt(errL2);
	//	dof_list[itr] = IDBC.size();
	//	h_max[itr] = tt3.MaxElementSize(bzmesh);
	//	err_list[itr] = errL2;
	//	//cout << "DOF: " << IDBC.size() << "\n";
	//	//cout << "L2-norm error: " << errL2 << "\n";
	//	//getchar();

	//	output_err(fld + fn + "_err", dof_list, err_list);
	//	output_err(fld + fn + "_err_h", h_max, err_list);
	//}
}

void kernel::run_coupling()
{
	string fn_in("../io/hex_input/cube4");
	//string fn_in("../io/hex_input/cube8");
	vector<BezierElement3D> bzmesh;
	vector<int> IDBC;
	TruncatedTspline_3D tt3;
	tt3.Run_Coupling(fn_in, 0, bzmesh, IDBC);

	//tt3.OutputMesh_Coupling(bzmesh, "../io/global1/cube8_couple");
	//cout << "Done output mesh!\n"; getchar();

	string fld("../io/global1/");
	//string fn("cube4_couple_1");
	string fn("cube4_RKPM_0");

	Laplace lap;
	lap.Set_Problem(IDBC);
	lap.Run_Coupling(bzmesh, fld + fn);
}

void kernel::run_coupling_Bezier()
{
	//string fn_in("../io/hex_input/cube4");
	//string fn_in("../io/hex_input/cube8");
	//string fn_in("../io/hex_input/cube_coarse");

	string fn_in("../io/hex_input/cube_ex0");
	//string fn_in("../io/hex_input/rod");
	//string fn_in("../io/hex_input/hook");

	vector<BezierElement3D> bzmesh;
	vector<int> IDBC;
	vector<double> gh;
	TruncatedTspline_3D tt3;
	tt3.Run_Coupling_Bezier(fn_in, 0, bzmesh, IDBC, gh);

	//tt3.Run_Coupling(fn_in, 2, bzmesh, IDBC);
	//tt3.OutputMesh_Coupling(bzmesh, "../io/global2/cube4");
	//tt3.OutputCM("../io/global2/cube_ex");
	//cout << "Done output mesh!\n"; getchar();

	//string fld("../io/global3/");
	string fld("../io/global4/");

	string fn("cubeex_Bezier_0");
	//string fn("rod_Bezier_3");
	//string fn("hook_Bezier_2");

	Laplace lap;
	lap.Set_Problem_Bezier(IDBC, gh);
	lap.Run_Coupling_Bezier(bzmesh, fld + fn);
}

void kernel::run_coupling_comp()
{
	//string fld("../io/navair_analysis/");
	//string fn_in(fld + "navair_input");
	//string fn_out(fld + "navair_0");


	/// Rod Model
	//string fld("../io/rod/");
	//string fn_in(fld + "rod_input");
	//string fn_out(fld + "rod_0");

	///// Cube test
	//string fld("../io/cube_test/");
	//string fn("cube_test");
	//string fn_in(fld + "cube_test");
	//string fn_out(fld + "cube_test_0");

	///// engine mount
	//string fld("../io/engine/");
	//string fn("engine_CM");
	//string fn_in(fld + fn);
	//string fn_out(fld + fn +"_output");

	///// BEXT3D
	//string fld("../io/BEXT3D/");
	//string fn("engine_CM");
	//string fn_in(fld + fn);
	//string fn_out(fld + fn + "_output");

	/// Plate with hole
	string fld("../io/PlateWithHole/");
	string fn_in(fld + "PlateWithHole_extrude2");
	string fn_out(fld + "PlateWithHole_0");

	vector<BezierElement3D> bzmesh;
	vector<int> IDBC;
	vector<double> gh;
	TruncatedTspline_3D tt3;

	tt3.InitializeMesh(fn_in);
	////tt3.SetSharpFeature_Manual(fld + "sharp.txt");
	tt3.SetSharpFeature_1(0.7);
	tt3.Run_Coupling_Comp(fn_in, 0, 0, bzmesh, IDBC, gh);

	//tt3.ReadBEXT3D("../io/BEXT3D/cube_test_BEXT.txt", bzmesh);
	//tt3.ReadBEXT3D("../io/BEXT3D/extrude_test_extrude_BEXT.txt", bzmesh, IDBC, gh);
	//tt3.ReadBEXT3D("../io/engine/engine_CM_output_BEXT.txt", bzmesh); 
	cout << "Read Done" << endl;
	//getchar();
	tt3.VisualizeBezierMesh(bzmesh, fn_out);

	//tt3.WriteBezier_Abaqus(bzmesh, fld + "_ABAQUS");
	tt3.WriteBezierInfo_AllSpline_LSDYNA(fn_out + "_BEXT", bzmesh);


	Laplace lap;
	//lap.Set_Problem_Bezier(IDBC, gh);
	//lap.Run_Coupling_Bezier(bzmesh, fn_out);
	lap.VisualizeDisp_Coupling_Bezier(bzmesh, "../io/BEXT3D/extrude_test_extrude_BEXT");
	lap.VisualizeVTK(bzmesh, "../io/PlateWithHole/PlateWithHole");
	//lap.VisualizeVTK(bzmesh, "../io/BEXT3D/engine_test");

	//LinearElasticity le;
	//le.Set_Problem_Bezier(IDBC, gh);
	//le.Run_Coupling_Bezier(bzmesh, fld + fn);
	//le.Run_MKL(bzmesh, fld + fn);
}


void kernel::run_UB_adapt_fit()
{
	int niter(3);
	double thresh(0.1);

	//string fn_in("../io/hex_input/heli_coarse_loc1");
	//string fn_in("../io/hex_input/cube_ex0");
	string fn_in("../io/hex_input/pillow/navair_coarse_CM");

	TruncatedTspline_3D tt3;
	tt3.Initialize_AdaptFit(fn_in);

	//string fld("../io/le_adapt/heli_coarse3/");
	//string fld("../io/le_adapt/cube_ex/");
	string fld("../io/le_adapt/navair_tmp/");

	//string fn("heli_");
	//string fn("cube_");
	string fn("navair_");

	vector<int> ncp_list(niter, 0);
	vector<int> nel_list(niter, 0);
	//vector<double> h_max(niter, 0.);
	vector<double> vsmax_list(niter, 0.);

	int itr(0);
	while (itr < niter)
	{
		stringstream ss;
		ss << itr;

		vector<BezierElement3D> bzmesh;
		vector<int> IDBC;
		vector<double> gh;
		//tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);
		tt3.BezierExtract_AdaptFit(bzmesh);
		tt3.StrongDirichletBC_Elasticity_AdaptFit(IDBC, gh);

		//tt3.VisualizeBC(fld + fn + ss.str(), IDBC, gh);
		//cout << "Done visualize BC\n";
		//getchar();

		vector<double> vs;
		double vsmax;
		LinearElasticity le;
		le.Set_Problem_Bezier(IDBC, gh);
		le.Run_UB_AdaptFit(bzmesh, fld + fn + ss.str(), vs, vsmax);

		ncp_list[itr] = IDBC.size() / 3;
		nel_list[itr] = bzmesh.size();
		vsmax_list[itr] = vsmax;

		if (itr > 0)
		{
			double vsmax_diff = fabs(vsmax - vsmax_list[itr - 1]);
			if (vsmax_diff < thresh)
			{
				cout << "Convergent!\n";
				//break;
			}
		}

		//Laplace lap;
		//lap.SetProblem(IDBC, gh);
		//stringstream ss;
		//ss << itr;
		////cout << "before simulation...\n";
		////getchar();
		////lap.GetRemoveRegion(remove);
		////lap.GetEqParameter(xy, nm, a);
		//lap.Run(bzmesh, fld + fn + ss.str(), err);
		////cout << "simulation done!\n";
		////getchar();

		////tt3.VisualizeBezier(bzmesh, fld + fn + ss.str());

		////double errL2(0.);
		//errL2 = 0.;
		//for (i = 0; i < err.size(); i++) errL2 += err[i];
		//errL2 = sqrt(errL2);
		//dof_list[itr] = IDBC.size();
		////h_max[itr] = tt3.MaxElementSize(bzmesh);
		//err_list[itr] = errL2;
		////cout << "DOF: " << IDBC.size() << "\n";
		////cout << "L2-norm error: " << errL2 << "\n";
		////getchar();

		output_err(fld + fn + ss.str() + "_cp", ncp_list, vsmax_list);
		output_err(fld + fn + ss.str() + "_el", nel_list, vsmax_list);
		//output_err(fld + fn + ss.str() + "_err_h", h_max, err_list);

		if (itr < niter - 1)
		{
			cout << itr << " refining...\n";
			//distribute error
			vector<array<double, 2>> eh(bzmesh.size());
			for (unsigned int i = 0; i < bzmesh.size(); i++)
			{
				eh[i][0] = bzmesh[i].prt[0]; eh[i][1] = bzmesh[i].prt[1];
			}

			vector<array<int, 2>> rfid, gst;
			tt3.Identify_Elasticity_AdaptFit(eh, vs, rfid, gst);
			//tt3.Identify_Refine_Test(rfid, gst);

			//tt3.OutputRefineID(fld + fn + ss.str(), rfid, gst);
			//tt3.InputRefineID("../io/benchmark1/cube1_" + ss.str(), rfid, gst);

			tt3.Refine(rfid, gst);
			cout << "Refining done\n";
			//tt3.OutputGeom_All(fld + fn + ss.str() + "_geom");
			//cout << "Output Geom done!\n";
			//getchar();
		}
		itr++;
	}
}

void kernel::run_THS3D()
{
	int niter(20);
	double thresh(0.1);//heli

	string fn_in("../io/hex_input/heli_coarse");
	//string fn_in("../io/hex_input/cube_coarse");

	TruncatedTspline_3D tt3;
	tt3.Initialize_THS3D(fn_in);

	string fld("../io/le_adapt/heli_coarse2/");
	//string fld("../io/le_adapt/cube_coarse/");

	string fn("heli_");
	//string fn("cube_");

	vector<int> ncp_list(niter, 0);
	vector<int> nel_list(niter, 0);
	//vector<double> h_max(niter, 0.);
	vector<double> vsmax_list(niter, 0.);

	int itr(0);
	while (itr < niter)
	{
		stringstream ss;
		ss << itr;

		vector<BezierElement3D> bzmesh;
		vector<int> IDBC;
		vector<double> gh;
		//tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);
		tt3.BezierExtract_AdaptFit(bzmesh);
		tt3.StrongDirichletBC_Elasticity_AdaptFit(IDBC, gh);
		tt3.VisualizeBC(fld + fn + ss.str(), IDBC, gh);
		cout << "Done visualize BC\n";
		getchar();

		vector<double> vs;
		double vsmax;
		LinearElasticity le;
		le.Set_Problem_Bezier(IDBC, gh);
		le.Run_UB_AdaptFit(bzmesh, fld + fn + ss.str(), vs, vsmax);

		ncp_list[itr] = IDBC.size() / 3;
		nel_list[itr] = bzmesh.size();
		vsmax_list[itr] = vsmax;

		if (itr > 0)
		{
			double vsmax_diff = fabs(vsmax - vsmax_list[itr - 1]);
			if (vsmax_diff < thresh)
			{
				cout << "Convergent!\n";
				//break;
			}
		}

		//Laplace lap;
		//lap.SetProblem(IDBC, gh);
		//stringstream ss;
		//ss << itr;
		////cout << "before simulation...\n";
		////getchar();
		////lap.GetRemoveRegion(remove);
		////lap.GetEqParameter(xy, nm, a);
		//lap.Run(bzmesh, fld + fn + ss.str(), err);
		////cout << "simulation done!\n";
		////getchar();

		////tt3.VisualizeBezier(bzmesh, fld + fn + ss.str());

		////double errL2(0.);
		//errL2 = 0.;
		//for (i = 0; i < err.size(); i++) errL2 += err[i];
		//errL2 = sqrt(errL2);
		//dof_list[itr] = IDBC.size();
		////h_max[itr] = tt3.MaxElementSize(bzmesh);
		//err_list[itr] = errL2;
		////cout << "DOF: " << IDBC.size() << "\n";
		////cout << "L2-norm error: " << errL2 << "\n";
		////getchar();

		output_err(fld + fn + ss.str() + "_cp", ncp_list, vsmax_list);
		output_err(fld + fn + ss.str() + "_el", nel_list, vsmax_list);
		//output_err(fld + fn + ss.str() + "_err_h", h_max, err_list);

		if (itr < niter - 1)
		{
			cout << itr << " refining...\n";
			//distribute error
			vector<array<double, 2>> eh(bzmesh.size());
			for (unsigned int i = 0; i < bzmesh.size(); i++)
			{
				eh[i][0] = bzmesh[i].prt[0]; eh[i][1] = bzmesh[i].prt[1];
			}

			vector<array<int, 2>> rfid, gst;
			tt3.Identify_Elasticity_AdaptFit(eh, vs, rfid, gst);
			//tt3.OutputRefineID(fld + fn + ss.str(), rfid, gst);
			//tt3.InputRefineID("../io/benchmark1/cube1_" + ss.str(), rfid, gst);
			tt3.Refine(rfid, gst);
			cout << "Refining done\n";
			//tt3.OutputGeom_All(fld + fn + ss.str() + "_geom");
			//cout << "Output Geom done!\n";
			//getchar();
		}
		itr++;
	}
}

void kernel::run_HS3D()
{
	int niter(4);
	//double thresh(0.1);//heli
	int nStep(60);
	double stepSize(0.3);

	//string fn_in("../io/hex_input/heli_coarse");
	//string fn_in("../io/hex_input/cube_coarse");
	//string fn_in("../io/hex_input/cube_ex0");
	string fn_in("../io/hex_input/smooth/honda1_CM");
	//string fn_in("../io/hex_input/smooth/honda2_CM");

	//string fld("../io/le_adapt/heli_coarse2/");
	//string fld("../io/le_adapt/cube_coarse/");
	//string fld("../io/le_adapt/cube_ex/");
	string fld("../io/le_adapt/honda1/");
	//string fld("../io/le_adapt/honda2/");

	//string fn("heli_");
	//string fn("cube_");
	string fn("honda1_tmp_");
	//string fn("honda2_");

	TruncatedTspline_3D tt3;
	tt3.Initialize_AdaptFit(fn_in);

	//vector<int> ncp_list(niter, 0);
	//vector<int> nel_list(niter, 0);
	//vector<double> vsmax_list(niter, 0.);

	int itr(0);
	bool badflag(true);
	while (itr < niter && badflag)
	{
		stringstream ss;
		ss << itr;
		cout << itr << "-th refining...\n";
		tt3.Smoothing_Adapt(nStep, stepSize);
		vector<array<int, 2>> rfid, gst;
		tt3.Identify_BadElement(rfid, gst);

		tt3.BadElementFlag_Adapt(rfid);
		tt3.OutputCM(itr, fld + fn + ss.str());

		if (rfid.size() == 0) badflag = false;
		if (itr < niter - 1)
		{
			tt3.Refine(rfid, gst);
		}
		cout << "\nRefining done\n";
		itr++;
	}

	//vector<BezierElement3D> bzmesh;
	//vector<int> IDBC;
	//vector<double> gh;
	////tt3.AnalysisInterface_Poisson_1(bzmesh, IDBC, gh);
	//tt3.BezierExtract_AdaptFit(bzmesh);
	////tt3.StrongDirichletBC_Elasticity_AdaptFit(IDBC, gh);
	//tt3.StrongDirichletBC_Elasticity_PID_Adapt("../io/le_adapt/honda1/bc", IDBC, gh);

	//vector<double> vs;
	//double vsmax;
	//LinearElasticity le;
	//le.Set_Problem_Bezier(IDBC, gh);
	//le.Run_UB_AdaptFit(bzmesh, fld + fn, vs, vsmax);
	////le.Run_MKL(bzmesh, fld + fn);
}





void kernel::run_C0C1Bezier()
{
	//string fn_in("../io/hex_input/cube4");
	//string fn_in("../io/hex_input/cube8");
	//string fn_in("../io/hex_input/cube_coarse");

	string fn_in("../io/hex_input/cube_ex0");
	//string fn_in("../io/hex_input/rod");
	//string fn_in("../io/hex_input/hook");
	//string fn_in("../io/hex_input/base");
	//string fn_in("../io/hex_input/heli");
	//string fn_in("../io/hex_input/heli_coarse");
	//string fn_in("../io/hex_input/rockerarm");

	vector<BezierElement3D> bzmesh;
	vector<int> IDBC;
	vector<double> gh;
	TruncatedTspline_3D tt3;
	tt3.Run_C0C1Bezier(fn_in, 0, bzmesh, IDBC, gh);

	//tt3.Run_Coupling(fn_in, 2, bzmesh, IDBC);
	//tt3.OutputMesh_Coupling(bzmesh, "../io/global2/cube4");
	//tt3.OutputCM("../io/global2/cube_ex");
	//cout << "Done output mesh!\n"; getchar();

	//string fld("../io/global3/");
	//string fld("../io/global4/");
	//string fld("../io/global5/");
	//string fld("../io/global8/");
	//string fld("../io/ptest/");
	string fld("../io/cube/");
	//string fld("../io/strongbc3/");//exp
	//string fld("../io/letest/");
	//string fld("../io/letest1/");
	//string fld("../io/poissontest/");
	//string fld("../io/Abaqus_inp/");

	//string fn("cubeex_Bezier_z");
	//string fn("rod_Bezier_3");
	//string fn("hook_Bezier_2");
	//string fn("base_Bezier_2");
	string fn("cubeex_c012_0");
	//string fn("heli_c012_0");
	//string fn("helicoarse_c012_0");
	//string fn("rockerarm_c012_0");

	//tt3.WriteBezier_Abaqus(bzmesh, fld + fn);

	Laplace lap;
	lap.Set_Problem_Bezier(IDBC, gh);
	lap.Run_Coupling_Bezier(bzmesh, fld + fn);

	//LinearElasticity le;
	//le.Set_Problem_Bezier(IDBC, gh);
	//le.Run_Coupling_Bezier(bzmesh, fld + fn);
}


void kernel::run_C0Bezier()
{
	//string fn_in("../io/hex_input/cube4");
	//string fn_in("../io/hex_input/cube8");
	//string fn_in("../io/hex_input/cube_coarse");

	string fn_in("../io/hex_input/cube_ex0");
	//string fn_in("../io/hex_input/rod");
	//string fn_in("../io/hex_input/hook");
	//string fn_in("../io/hex_input/base");

	vector<BezierElement3D> bzmesh;
	vector<int> IDBC;
	vector<double> gh;
	TruncatedTspline_3D tt3;
	tt3.Run_C0Bezier(fn_in, 0, bzmesh, IDBC, gh);

	//tt3.Run_Coupling(fn_in, 2, bzmesh, IDBC);
	//tt3.OutputMesh_Coupling(bzmesh, "../io/global2/cube4");
	//tt3.OutputCM("../io/global2/cube_ex");
	//cout << "Done output mesh!\n"; getchar();

	//string fld("../io/global3/");
	//string fld("../io/global4/");//sins
	//string fld("../io/global5/");//monomials
	//string fld("../io/global6/");//monomials, all C0
	//string fld("../io/global7/");//sin, all C0
	//string fld("../io/global8/");//sin, all C0, enlarge
	string fld("../io/strongbc4/");//exp
	//string fld("../io/weakbc1/");//exp

	string fn("cubeex_Bezier_0");
	//string fn("rod_Bezier_3");
	//string fn("hook_Bezier_2");
	//string fn("base_Bezier_2");

	Laplace lap;
	lap.Set_Problem_Bezier(IDBC, gh);
	lap.Run_Coupling_Bezier(bzmesh, fld + fn);
}


void kernel::run_C02iBezier()
{
	//string fn_in("../io/hex_input/cube_ex0");
	//string fn_in("../io/hex_input/rod");
	//string fn_in("../io/hex_input/hook");
	//string fn_in("../io/hex_input/base");
	//string fn_in("../io/yuxuan/heli");
	string fn_in("../io/hex_input/heli_coarse");
	//string fn_in("../io/yuxuan/navair");
	//string fn_in("../io/yuxuan/fertility");
	//string fn_in("../io/yuxuan/rockerarm");

	vector<BezierElement3D> bzmesh;
	vector<int> IDBC;
	vector<double> gh;
	TruncatedTspline_3D tt3;
	tt3.Run_C02iBezier(fn_in, 0, bzmesh, IDBC, gh);

	//tt3.PatchTest_OneElement(bzmesh, IDBC, gh);

	//tt3.WriteBezierInfo_Blend_LSDYNA("../io/LSDYNA/cube_ds", bzmesh);
	//tt3.WriteBezier_Abaqus(bzmesh,"../io/letest/rockerarm");
	//tt3.WriteBezier_Abaqus(bzmesh, "../io/letest/navair");
	//tt3.WriteBezier(bzmesh, "../io/letest/cube");//BEXT format
	//tt3.WriteBezier(bzmesh, "../io/yuxuan/heli/heli_bezier.txt");//BEXT format
	//tt3.WriteBezier(bzmesh, "../io/yuxuan/fertility/fertility_bezier.txt");
	//tt3.WriteBezier(bzmesh, "../io/yuxuan/rockerarm/rockerarm_bezier.txt");
	//cout << "Done output Abaqus!\n"; getchar();

	//tt3.Run_Coupling(fn_in, 2, bzmesh, IDBC);
	//tt3.OutputMesh_Coupling(bzmesh, "../io/global2/cube4");
	//tt3.OutputCM("../io/global2/cube_ex");
	//cout << "Done output mesh!\n"; getchar();

	//string fld("../io/global3/");
	//string fld("../io/global4/");
	//string fld("../io/global5/");
	//string fld("../io/global8/");
	//string fld("../io/global9/");
	//string fld("../io/ptest/");
	//string fld("../io/strongbc1/");//sin, C02i
	//string fld("../io/strongbc2/");//exp, C02i
	//string fld("../io/yuxuan/heli/");
	//string fld("../io/yuxuan/navair/");
	//string fld("../io/yuxuan/fertility/");
	//string fld("../io/yuxuan/rockerarm/");
	//string fld("../io/letest/");
	//string fld("../io/poissontest/");
	string fld("../io/letest2/");

	//string fn("cubeex_Bezier_3");
	//string fn("rod_Bezier_3");
	//string fn("hook_Bezier_2");
	//string fn("base_Bezier_2");

	//string fn("cubeex_c02i_1");
	//string fn("heli_c02i_1");
	string fn("heli_coarse_c02i_0");
	//string fn("navair_c02i_1");
	//string fn("fertility_c02i_0");
	//string fn("rockerarm_c02i_1");

	//tt3.WriteBezier_Abaqus(bzmesh, fld + fn);

	//Laplace lap;
	//lap.Set_Problem_Bezier(IDBC, gh);
	//lap.Run_Coupling_Bezier(bzmesh, fld + fn);

	LinearElasticity le;
	le.Set_Problem_Bezier(IDBC, gh);
	le.Run_Coupling_Bezier(bzmesh, fld + fn);
}






//////////////////////////////////////////////////////////////////////

void kernel::run_UBsplines()
{
	//string fn_in("../io/hex_input/cube4");
	//string fn_in("../io/hex_input/cube8");
	//string fn_in("../io/hex_input/cube_coarse");
	//string fn_in("../io/hex_input/cube_ex0");
	//string fn_in("../io/hex_input/rod");
	//string fn_in("../io/hex_input/hook");
	//string fn_in("../io/hex_input/base");
	//string fn_in("../io/LSDYNA/statue");
	//string fn_in("../io/LSDYNA/gear");
	//string fn_in("../io/hex_input/fertility");
	//string fn_in("../io/hex_input/heli");
	string fn_in("../io/hex_input/rockerarm");
	//string fn_in("../io/hex_input/navair");

	vector<BezierElement3D> bzmesh;
	vector<int> IDBC;
	vector<double> gh;
	TruncatedTspline_3D tt3;
	tt3.Run_Coupling_Comp(fn_in, 0, 0, bzmesh, IDBC, gh);

	//string fld("../io/LSDYNA/");
	//string fn_out("statue_ds");
	//string fn_out("gear_ds_nofeature");
	//string fn_out("gear_ds_feature");
	//tt3.WriteBezierInfo_AllSpline_LSDYNA(fld + fn_out, bzmesh);

	//tt3.Run_Coupling(fn_in, 2, bzmesh, IDBC);
	//tt3.OutputMesh_Coupling(bzmesh, "../io/global2/cube4");
	//tt3.OutputEdge("../io/global4/cube");
	//tt3.OutputFace("../io/global4/cube");
	//tt3.OutputCM("../io/global2/cube_ex");
	//cout << "Done output mesh!\n"; getchar();

	//string fld("../io/global3/");
	//string fld("../io/global4/");
	//string fld("../io/ptest/");
	//string fld("../io/models/");
	//string fld("../io/yuxuan/fertility/");
	//string fld("../io/yuxuan/heli/");
	string fld("../io/usp1/rockerarm/");
	//string fld("../io/yuxuan/navair/");

	//string fn("cubeex_Spline_4");
	//string fn("rod_Spline_3");
	//string fn("hook_Spline_1");
	string fn("rockerarm_Spline_0");
	//string fn("fertility_Spline_1");
	//string fn("heli_Spline_1");
	//string fn("navair_Spline_1");

	//tt3.OutputCM(fld + fn);
	//tt3.WriteBezier_Abaqus(bzmesh, fld + fn);

	//Laplace lap;
	//lap.Set_Problem_Bezier(IDBC, gh);
	//lap.Run_Coupling_Bezier(bzmesh, fld + fn);

	//LinearElasticity le;
	//le.Set_Problem_Bezier(IDBC, gh);
	//le.Run_Coupling_Bezier(bzmesh, fld + fn);
}









//////////////////////////////////////////////////////////////////////////////
void kernel::run_MeshQualityImprove()
{
	string fld("../io/rod_smooth/");	
	string fn_in(fld + "rod_pillow");
	////string fn_in(fld + "navair_lap_CM");
	////string fn_in(fld + "navair_smooth_CM");
	////string fn_in(fld + "navair_pillow_CM");
	//
	////string fn_out(fld + "navair_lap");
	string fn_out(fld + "rod_smooth");
	//string fn_out(fld + "navair_opt");
	////string fn_out(fld + "navair_pillow");

	//string fld("../io/central_pillar/");
	//string fn_in(fld + "central_pillar_volume");
	//string fn_out(fld + "central_pillar_volume_smooth");

	//string fld("../io/engine_bad/");	
	//string fn_in(fld + "engine_input");
	//string fn_in(fld + "engine_good_input");
	//string fn_in(fld + "engine_lap");
	//string fn_in(fld + "engine_smooth");
	//string fn_in(fld + "engine_pillow");

	//string fn_out(fld + "engine_lap");
	//string fn_out(fld + "engine_pillow");
	//string fn_out(fld + "engine_smooth");
	//string fn_out(fld + "engine_good_opt");
	//string fn_out(fld + "engine_opt");

	TruncatedTspline_3D tt3;
	tt3.InitializeMesh(fn_in);

	//tt3.LaplaceSmoothing(50);
	//tt3.OutputCM(fn_out);
	//cout << "done lap!\n";
	//getchar();



	//tt3.SetSharpFeature_1(0.7);
	tt3.SetSharpFeature_Manual(fld + "sharp_pillow.txt");


	////tt3.SetSharpFeature_1();//automatic way may not work
	string fn1(fld + "navair_sharp");
	tt3.OutputEdge(fn1);
	tt3.OutputCM(fn1);
	cout << "done setting sharp feature\n";
	//getchar();

	//tt3.CheckJacobBEXT();

	//mesh quality improvement includes three parts (Pillow, Smoothing and Optimizizing) and they should be run one by one

	//tt3.Pillow(1);//run Pillow first if needed, which outputs an mesh with an outer layer (input argument "1" means one layer) added to the input mesh.
	//			//the name of the outpue mesh can be changed within the function definition
	//tt3.OutputCM(fn_out);

	//tt3.Smoothing(100, 0.1);//then run Smoothing on the pillowed mesh by specifying the
							   //number of steps (e.g. 350) and the step size (e.g. 0.5), which
							   //iteratively moves each point in the mesh to its mass center.
							   //The input "fn_in" should be changed to the pillowed mesh or leave it be if smoothing is applied without pillowing
							   //It outputs a smoothed mesh, see the commented code "OutputCM" in the definition
	tt3.Smoothing_Angran(1, 0.0001);
	tt3.OutputCM(fn_out);




	//tt3.Optimizing(50, 0.1);//last run Optimizing if there is still negative Jacobian elements after smoothing
	//	                      //specifying the number of steps (e.g. 100) and the step size (e.g. 0.01),
	//						//which in each iteration quality improves the element of the worst Jacobian.
	//						//The input "fn_in" should be changed to the smoothed mesh.
	//						//It outputs an optimized mesh, check the commented code "OutputCM" in the function definition.

	//tt3.OutputCM(fn_out);

}


void kernel::run_coupling_comp(int flag_sharp, int ngref, int nlref, double tol_sharp, string fn)
{
	char *buffer1 = strdup(fn.c_str());
	char *buffer2 = strdup(fn.c_str());
	char *fld_char = buffer1;
	char *fn_wo_extension_char = buffer2;
	string fld;
	string fn_in;
	PathRemoveFileSpecA(fld_char);
	PathRemoveExtensionA(fn_wo_extension_char);
	stringstream ss1;
	ss1 << fld_char;
	ss1 >> fld;
	stringstream ss2;
	ss2 << fn_wo_extension_char;
	ss2 >> fn_in;

	stringstream ss_out, ss_out_GEM;
	string fn_out;
	string fn_out_GEM = fld + "\\Abq_Spline";

	// * CMU naming rule: refinement level
	ss_out << fn_in << "_gref_lev" << ngref << "_lref_lev" << nlref;
	ss_out >> fn_out;
	// * GEM naming rele: abaqus
	//ss_out_GEM << fld << "_gref_lev" << ngref << "_lref_lev" << nlref;
	//ss_out_GEM >> fn_out_GEM;

	cout << "Input Mesh: "<< fn << endl;

	vector<BezierElement3D> bzmesh;
	vector<int> IDBC;
	vector<double> gh;
	TruncatedTspline_3D tt3;

	tt3.InitializeMesh(fn_in);

	switch (flag_sharp)
	{
	default:
		cout << "Sharp feature OFF" << endl;
		break;
	case 1:
		cout << "Automatically detect sharp feature ON: TOL = " << tol_sharp << endl;
		if (nlref == 0)
		{
			tt3.SetSharpFeature_1(tol_sharp);
		}
		else
		{
			tt3.SetSharpFeature_localrefine(tol_sharp);
		}
		
		break;
	case 2:
		cout << "Manually apply sharp feature ON" << endl;
		if (nlref == 0)
		{
			//tt3.SetSharpFeature_Manual(fn_in + "_sharp.txt");
			tt3.SetSharpFeature_Manual(fld + "\\sharp.txt");
			//tt3.SetSharpFeature_Manual_Yu("sharp.txt");
		}
		else
		{
			//tt3.SetSharpFeature_Manual_Yu(fn_in + "_sharp.txt");
			tt3.SetSharpFeature_Manual_Yu(fld + "\\sharp.txt");
		}

		break;
	}

	if (nlref != 0)
	{
		tt3.AllBezierLev(0);
		tt3.BezierPoints_ini();//represent geometry using initial control points and blending functions
		tt3.SetSupport(0);
		tt3.Run_Coupling_Comp(fn_in, ngref, nlref, bzmesh, IDBC, gh);
		//tt3.VisualizeBezierMesh(bzmesh, fn_out);
		tt3.WriteBezierInfo_AllSpline_LSDYNA_LocalRefine(fn_out + "_BEXT", bzmesh);
		//tt3.WriteBezier_Abaqus_LocalRefine(bzmesh, fn_out + "_ABAQUS");
	}


	if (nlref == 0)
	{
		tt3.Run_Coupling_Comp(fn_in, ngref, nlref, bzmesh, IDBC, gh);
		//tt3.VisualizeBezierMesh(bzmesh, fn_out);
		tt3.OutputCM(fn_out + "_CM");
		tt3.WriteBezierInfo_AllSpline_LSDYNA(fn_out + "_BEXT", bzmesh);
		tt3.WriteBezier_Abaqus(bzmesh, fn_out_GEM);
	}


	//tt3.WriteBezier_Abaqus(bzmesh, fld + "_ABAQUS");
	//


	//tt3.AngranOutputMesh(bzmesh, fn_out + "_Angran_");

	Laplace lap;
	////lap.Set_Problem_Bezier(IDBC, gh);
	////lap.Run_Coupling_Bezier(bzmesh, fn_out);
	lap.VisualizeVTK(bzmesh, fn_out);

	//LinearElasticity le;
	//le.Set_Problem_Bezier(IDBC, gh);
	////le.Run_Coupling_Bezier(bzmesh, fld + fn);
	//le.Run_MKL(bzmesh, fld + fn);
}

void kernel::run_MeshQualityImprove(int mode, int flag_sharp, double tol_sharp, int opt_par1, double opt_par2, string fn)
{
	//string fld("../io/navair/");	
	////string fn_in(fld + "navair_input");
	////string fn_in(fld + "navair_lap_CM");
	////string fn_in(fld + "navair_smooth_CM");
	////string fn_in(fld + "navair_pillow_CM");
	//
	////string fn_out(fld + "navair_lap");
	////string fn_out(fld + "navair_smooth");
	//string fn_out(fld + "navair_opt");
	////string fn_out(fld + "navair_pillow");

	char *buffer1 = strdup(fn.c_str());
	char *buffer2 = strdup(fn.c_str());
	char *fld_char=buffer1;
	char *fn_wo_extension_char= buffer2;
	string fld;
	string fn_out;
	PathRemoveFileSpecA(fld_char);
	PathRemoveExtensionA(fn_wo_extension_char);
	stringstream ss1;
	ss1 << fld_char;
	ss1 >> fld;
	stringstream ss2;
	ss2 << fn_wo_extension_char;
	ss2 >> fn_out;

	cout << "Input Mesh: " << fn << endl;
	//string fld("../io/central_pillar/");
	//string fn_in(fld + "central_pillar_volume");
	//string fn_out(fld + "central_pillar_volume_smooth");




	TruncatedTspline_3D tt3;
	tt3.InitializeMesh(fn_out);

	//tt3.LaplaceSmoothing(50);
	//tt3.OutputCM(fn_out);
	//cout << "done lap!\n";
	//getchar();

	switch (flag_sharp)
	{
	default:
		cout << "Sharp feature OFF" << endl;
		break;
	case 1:
		cout << "Automatically detect sharp feature ON: TOL = "<< tol_sharp  << endl;
		tt3.SetSharpFeature_1(tol_sharp);
		break;
	case 2:
		cout << "Manually apply sharp feature ON" << endl;
		tt3.SetSharpFeature_Manual(fld + "\\sharp.txt");		
		break;
	}
		
	////tt3.SetSharpFeature_1();//automatic way may not work
	//string fn1(fld + "navair_sharp");
	//tt3.OutputEdge(fn1);
	//tt3.OutputCM(fn1);
	//cout << "done setting sharp feature\n";
	//getchar();

	//tt3.CheckJacobBEXT();

	//mesh quality improvement includes three parts (Pillow, Smoothing and Optimizizing) and they should be run one by one
	
	switch (mode)
	{
	case 0: //Laplace Smoothing
		tt3.LaplaceSmoothing(opt_par1);
		tt3.OutputCM(fn_out + "_lap");
		cout << "done lap!\n";
		break;
	case 1: //Pillowing
		///tt3.Pillow(1);//run Pillow first if needed, which outputs an mesh with an outer layer (input argument "1" means one layer) added to the input mesh.
					  ///the name of the outpue mesh can be changed within the function definition
		tt3.Pillow(opt_par1);
		tt3.OutputCM(fn_out + "_pillow");
		break;
	case 2: //Smoothing
		///tt3.Smoothing(50, 0.4);//then run Smoothing on the pillowed mesh by specifying the
		///					   //number of steps (e.g. 350) and the step size (e.g. 0.5), which
		///					   //iteratively moves each point in the mesh to its mass center.
		///					   //The input "fn_in" should be changed to the pillowed mesh or leave it be if smoothing is applied without pillowing
		///					   //It outputs a smoothed mesh, see the commented code "OutputCM" in the definition
		//tt3.Smoothing(opt_par1, opt_par2);
		tt3.Smoothing_Angran(opt_par1, opt_par2);
		tt3.OutputCM(fn_out + "_smooth");
		break;
	case 3: //Optimizing
		///tt3.Optimizing(50, 0.1);//last run Optimizing if there is still negative Jacobian elements after smoothing
		///                        //specifying the number of steps (e.g. 100) and the step size (e.g. 0.01),
		///						//which in each iteration quality improves the element of the worst Jacobian.
		///						//The input "fn_in" should be changed to the smoothed mesh.
		///						//It outputs an optimized mesh, check the commented code "OutputCM" in the function definition.
		tt3.Optimizing(opt_par1, opt_par2);
		tt3.OutputCM(fn_out + "_opt");
		break;
	}

	

	

	
}







//////////////////////////////////////////////////////////////////

void kernel::run_AbaqusGEM()
{
	string fn_in("../io/PlateHole/input_CM_3");
	TruncatedTspline_3D tt3;
	tt3.GeneratePlateHole3D(fn_in);
}