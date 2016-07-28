#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/tokenizer.hpp>
#include <boost/token_functions.hpp>

using namespace boost;
using namespace boost::program_options;

#include <iostream>
#include <fstream>
#include <exception>
#include <cstdlib>
#include "Caller.h"

using namespace std;

int main( int argc, char** argv)
{
	// Add descriptive text for display when help argument is supplied
	options_description desc( "\nSiNVICT: Ultra Sensitive Detection of Single Nucleotide Variants and Indels in Circulating Tumour DNA.\n\nAllowed arguments");

	// Define command line arguments using on of the following formats:
	// ("long-name, short-name", "Description of argument") for flag values
	// ("long-name, short-name", <data-type>, "Description of argument") for arguments with values
	
	desc.add_options()
		( "help,h", "Print help message.")
		( "error-rate,e", value<double>(), "Error rate for the sequencing platform used.")
		( "min-depth,m", value<int>(), "Minimum Read Depth required for high confidence in a call.")
		( "left-strand-bias,l", value<double>(), "Lower limit for the strand bias value interval to be used in assessing the confidence of a call.")
		( "right-strand-bias,r", value<double>(), "Upper limit for the strand bias value interval to be used in assessing the confidence of a call.")
		( "read-end-fraction,f", value<double>(), "Average position of the called base on the reads supporting the call as a fraction. End values such as 0.01 as useful for filtering read end artifacts.")
		( "qscore-cutoff,q", value<int>(), "Cutoff value for the qScore assigned to each call by the Poisson model used.")
		( "tumor-directory-path,t", value< vector<string> >(), "Specifies directory for the input files.")
		( "output-directory-path,o", value< vector<string> >(), "Specifies directory for the output files.");
	
	// Parse the command line, catching and displaying any parser errors
	variables_map vm;
	try
	{
		store(command_line_parser(argc, argv).options(desc).run(), vm);
		notify(vm);
	} catch(std::exception &e)
	{
		cout << endl << e.what() << endl;
		cout << desc << endl;
	}

	double errorRate = 0.01;
	int minDepth = 100;
	double leftStrandBias = 0.3;
	double rightStrandBias = 0.7;
	double readEndFraction = 0.01;
	int qScoreCutoff = 95;
	string tumorDirectoryPath;
	string outputDirectoryPath;
	bool tumorDirectorySpecified = false;
	bool outputDirectorySpecified = false;

	// Display the state of the arguments supplied
	if(vm.count("help"))
	{
		cout << desc << endl;
		cout << "  --tumor-directory-path and --output-directory-path must be specified" << endl << endl;
		std::cout << "  Usage: ./sinvict -e [error-rate] -m [min-depth] -l [left-strand-bias] -r [right-strand-bias] ";
		std::cout << "-f [read-end-fraction] -q [qscore-cutoff] -t <tumor-directory-path> -o <output-directory-path>" << std::endl;
		exit(1);
	}

	if(vm.count("tumor-directory-path"))
	{
		vector<string> inputPath = vm["tumor-directory-path"].as< vector<string> >();
		cout << "--tumor-directory-path specified with value = " << inputPath[0] << endl;
		tumorDirectoryPath.assign(inputPath[0]);
		tumorDirectorySpecified = true;
	}

	if(vm.count("output-directory-path"))
	{
		vector<string> outputPath = vm["output-directory-path"].as< vector<string> >();
		cout << "--output-directory-path specified with value = " << outputPath[0] << endl;
		outputDirectoryPath.assign(outputPath[0]);
		outputDirectorySpecified = true;
	}

	if( !tumorDirectorySpecified || !outputDirectorySpecified)
	{
		cout << desc << endl;
		cout << "  --tumor-directory-path and --output-directory-path must be specified" << endl << endl;
		std::cout << "  Usage: ./sinvict -e [error-rate] -m [min-depth] -l [left-strand-bias] -r [right-strand-bias] ";
		std::cout << "-f [read-end-fraction] -q [qscore-cutoff] -t <tumor-directory-path> -o <output-directory-path>" << std::endl;
		exit( 1);
	}
	else
	{
		if(vm.count("error-rate"))
		{
			errorRate = vm["error-rate"].as<double>();
			cout << "--error-rate set to: " << errorRate << endl;
		}
		else
		{
			cout << "--error-rate not set, using default value: " << errorRate << endl;
		}

		if(vm.count("min-depth"))
		{
			minDepth = vm["min-depth"].as<int>();
			cout << "--min-depth set to: " << minDepth << endl;
		}
		else
		{
			cout << "--min-depth not set, using default value: " << minDepth << endl;
		}

		if(vm.count("left-strand-bias"))
		{
			leftStrandBias = vm["left-strand-bias"].as<double>();
			cout << "--left-strand-bias set to: " << leftStrandBias << endl;
		}
		else
		{
			cout << "--left-strand-bias not set, using default value: " << leftStrandBias << endl;
		}

		if(vm.count("right-strand-bias"))
		{
			rightStrandBias = vm["right-strand-bias"].as<double>();
			cout << "--right-strand-bias set to: " << rightStrandBias << endl;
		}
		else
		{
			cout << "--right-strand-bias not set, using default value: " << rightStrandBias << endl;
		}

		if(vm.count("read-end-fraction"))
		{
			readEndFraction = vm["read-end-fraction"].as<double>();
			cout << "--read-end-fraction set to: " << readEndFraction << endl;
		}
		else
		{
			cout << "--read-end-fraction not set, using default value: " << readEndFraction << endl;
		}

		if(vm.count("qscore-cutoff"))
		{
			qScoreCutoff = vm["qscore-cutoff"].as<int>();
			cout << "--qscore-cutoff set to: " << qScoreCutoff << endl;
		}
		else
		{
			cout << "--qscore-cutoff not set, using default value: " << qScoreCutoff << endl;
		}

		Caller caller( errorRate, minDepth, leftStrandBias, rightStrandBias, readEndFraction, qScoreCutoff, tumorDirectoryPath.c_str(), outputDirectoryPath.c_str());
		caller.callLocationsMixture();
	}

	return 0;
}
