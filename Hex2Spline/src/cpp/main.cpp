#include <iostream>
//#include <vector>
//#include <array>
//#include <algorithm>
//#include "KnotInsertion.h"
//#include "TTSP_3D.h"
//#include "Elasiticity_3D.h"
//#include "Laplace.h"
//#include "LBO.h"
//#include <sstream>
//#include <iomanip>
#include "kernel.h"
#include "cxxopts.hpp"
using namespace std;

void Commandline(int argc, char* argv[]);

int main(int argc, char* argv[])
{
	//this piece of code is used to subdivide coarse polycube
	//TruncatedTspline_3D tt3;
	//tt3.InitializeMesh("../io/honda/honda1/honda1_paraHex_hex");
	//for(int i=0;i<2;i++)
	//	tt3.Global_Subdivide_Simple();
	//tt3.OutputCM("../io/honda/honda1/honda1_paraHex_hex1");
	//Vtk2Raw_hex("../io/honda/honda1/honda1_paraHex_hex1_CM");
	//Raw2Vtk_hex("../io/honda/honda1/honda1_paraHex_hex1_hex");

	kernel app;

	Commandline(argc, argv);

	//app.run_MeshQualityImprove();//mesh quality improvement

	//app.run_coupling_comp();//construction on the input hex mesh using TH-spline3D


	//////////////////////////////////////////////////////////////////////////////////////
	//app.run();
	//app.run_complex();
	//app.run_complex_fit();
	//app.run_complex_glb();
	//app.run_Bspline();
	//app.run_Bspline_fit();
	//app.run_benchmark();
	//app.run_leastsquare();
	//app.run_pipeline();

	//app.run_conv();
	//app.run_coupling();

	//app.run_coupling_Bezier();

	//app.run_UB_adapt_fit();
	//app.run_THS3D();//not done yet
	//app.run_HS3D();

	//app.run_C0C1Bezier();
	//app.run_C0Bezier();//do not support patch test, need modification
	//app.run_C02iBezier();

	//TruncatedTspline_3D tt3;
	//tt3.CreateRemoveView("../io/complex/rod1_3","../io/complex1/rod1_3");
	//tt3.PipelineDataProcess("../io/pipeline/complex_Enrich_Disp", "../io/pipeline/base_disp");
	
	//app.run_AbaqusGEM();
	//////////////////////////////////////////////////////////////////////////////////////////////

	/*cout<<"DONE!\n";
	getchar();
	return 0;*/
}

void Commandline(int argc, char* argv[])
{
	kernel app;

	try
	{
		cxxopts::Options options(argv[0], "CMU Solid Software");
		options
			.positional_help("[optional args]")
			.show_positional_help();

		bool flag_quality = false;

		int ngref = 0;
		int nlref = 0;
		int flag_quality_option = 0;
		int flag_sharp = 0;
		double tol_sharp;
		int quality_par1;
		double quality_par2;

		bool flag_spline = false;
		bool flag_analysis = false;
		bool flag_lref = false;

		string fn_in;
		string fn_out;

		options.add_options("General Settings")
			("h,help", "Print help")
			("Q,quality", "Mesh quality improvement mode", cxxopts::value<bool>(flag_quality))
			("S,spline", "Spline construction mode", cxxopts::value<bool>(flag_spline))
			("s,sharp", "0-No sharp feature, 1-Automatic sharp feature, 2-Manual sharp feature", cxxopts::value<int>(flag_sharp))
			("t,stol", "Tolerance for automatically detecting sharp feature", cxxopts::value<double>(tol_sharp))
			("I,input", "Input file", cxxopts::value<std::string>(fn_in))
#ifdef CXXOPTS_USE_UNICODE
			("unicode", u8"A help option with non-ascii: ид. Here the size of the"
				" string should be correct")
#endif
			;
		options.add_options("Mesh Quality Improvement")
			("m,method",
				"Improvement methods: 0-Laplacian Smoothing interior points (Give iternation number -n)\; 1-Pillowing (Give pillow layer number -n)\; 2-Smoothing (Give iteration number -n and smooth step -p)\; 3-Optimization (Give iteration number -n and optimization step -p)", cxxopts::value<int>(flag_quality_option))
				("n,number", "Pillowing layer number, Smoothing and Optimization number of steps", cxxopts::value<int>(quality_par1))
			("p,parameter", "Smoothing / Optimization step size", cxxopts::value<double>(quality_par2))
#ifdef CXXOPTS_USE_UNICODE
			("unicode", u8"A help option with non-ascii: ид. Here the size of the"
				" string should be correct")
#endif
			;
		options.add_options("Spline Construction")
			("g,globalref", "Set the level of global refinement, default is 0", cxxopts::value<int>(ngref))
			//("l,localref", "Set the level of local refinement, default is 0", cxxopts::value<int>(nlref))
			("l,localref", "Enable local refinement", cxxopts::value<bool>(flag_lref))
#ifdef CXXOPTS_USE_UNICODE
			("unicode", u8"A help option with non-ascii: ид. Here the size of the"
				" string should be correct")
#endif
			;

		auto result = options.parse(argc, argv);
		if (result.count("help"))
		{
			//cout << options.help({ "MeshQualityImprovement", "SplineConstruction", "Analysis" }) << std::endl;
			cout << options.help({ "General Settings", "Mesh Quality Improvement", "Spline Construction" }) << endl;
			exit(0);
		}

		if (flag_quality)
		{
			app.run_MeshQualityImprove(flag_quality_option, flag_sharp, tol_sharp, quality_par1, quality_par2, fn_in);//mesh quality improvement
		}

		if (flag_spline)
		{
			if (flag_lref)
				nlref = 1;
			app.run_coupling_comp(flag_sharp, ngref, nlref, tol_sharp, fn_in);//construction on the input hex mesh using TH-spline3D
		}
			
	}
	catch (const cxxopts::OptionException& e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}
}

