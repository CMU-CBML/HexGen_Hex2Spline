#include <iostream>
#include "kernel.h"
#include "cxxopts.hpp"
using namespace std;

void CommandlineHex2Spline(int argc, char *argv[]);

int main(int argc, char *argv[])
{
	kernel app;

	CommandlineHex2Spline(argc, argv);
}

void CommandlineHex2Spline(int argc, char *argv[])
{
	kernel app;

	try
	{
		cxxopts::Options options(argv[0], "Hex2Spline");
		options
			.positional_help("[optional args]")
			.show_positional_help();

		int ngref = 0;
		int nlref = 0;
		int flag_sharp = 0;
		double tol_sharp;

		bool flag_spline = false;
		bool flag_lref = false;

		string fn_in;
		string fn_out;

		options.add_options("General Settings")("h,help", "Print help")("S,spline", "Spline construction mode", cxxopts::value<bool>(flag_spline))("s,sharp", "0-No sharp feature, 1-Automatic sharp feature, 2-Manual sharp feature", cxxopts::value<int>(flag_sharp))("t,stol", "Tolerance for automatically detecting sharp feature", cxxopts::value<double>(tol_sharp))("I,input", "Input file", cxxopts::value<std::string>(fn_in))
#ifdef CXXOPTS_USE_UNICODE
			("unicode", u8"A help option with non-ascii: ��. Here the size of the"
						" string should be correct")
#endif
			;
		options.add_options("Spline Construction")("g,globalref", "Set the level of global refinement, default is 0", cxxopts::value<int>(ngref))
			//("l,localref", "Set the level of local refinement, default is 0", cxxopts::value<int>(nlref))
			("l,localref", "Enable local refinement", cxxopts::value<bool>(flag_lref))
#ifdef CXXOPTS_USE_UNICODE
				("unicode", u8"A help option with non-ascii: ��. Here the size of the"
							" string should be correct")
#endif
			;

		auto result = options.parse(argc, argv);
		if (result.count("help"))
		{
			//cout << options.help({ "MeshQualityImprovement", "SplineConstruction", "Analysis" }) << std::endl;
			cout << options.help({"General Settings", "Spline Construction"}) << endl;
			exit(0);
		}

		if (flag_spline)
		{
			if (flag_lref)
				nlref = 1;
			app.run_coupling_comp(flag_sharp, ngref, nlref, tol_sharp, fn_in); //construction on the input hex mesh using TH-spline3D
		}
	}
	catch (const cxxopts::OptionException &e)
	{
		std::cout << "error parsing options: " << e.what() << std::endl;
		exit(1);
	}
}