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
		( "errorRate,e", value<double>(), "Error rate for the sequencing platform used.")
		( "minDepth,m", value<int>(), "Minimum Read Depth required for high confidence in a call.")
		( "leftStrandBias,l", value<double>(), "Lower limit for the strand bias value interval to be used in assessing the confidence of a call.")
		( "rightStrandBias,r", value<double>(), "Upper limit for the strand bias value interval to be used in assessing the confidence of a call.")
		( "readEndFraction,f", value<double>(), "Average position of the called base on the reads supporting the call as a fraction. End values such as 0.01 as useful for filtering read end artifacts.")
		( "QScoreCutoff,q", value<int>(), "Cutoff value for the qScore assigned to each call by the Poisson model used.")
		( "tumorDirectoryPath,t", value< vector<string> >(), "Specifies directory for the input files.")
		( "outputDirectoryPath,o", value< vector<string> >(), "Specifies directory for the output files.");
	
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
		cout << "  --tumorDirectoryPath and --outputDirectoryPath must be specified" << endl << endl;
		std::cout << "  Usage: ./sinvict -e [errorRate] -m [minDepth] -l [leftStrandBias] -r [rightStrandBias] ";
		std::cout << "-f [readEndFraction] -q [QScoreCutoff] -t <tumorDirectoryPath> -o <outputDirectoryPath>" << std::endl;
		exit(1);
	}

	if(vm.count("tumorDirectoryPath"))
	{
		vector<string> inputPath = vm["tumorDirectoryPath"].as< vector<string> >();
		cout << "--tumorDirectoryPath specified with value = " << inputPath[0] << endl;
		tumorDirectoryPath.assign(inputPath[0]);
		tumorDirectorySpecified = true;
	}

	if(vm.count("outputDirectoryPath"))
	{
		vector<string> outputPath = vm["outputDirectoryPath"].as< vector<string> >();
		cout << "--outputDirectoryPath specified with value = " << outputPath[0] << endl;
		outputDirectoryPath.assign(outputPath[0]);
		outputDirectorySpecified = true;
	}

	if( !tumorDirectorySpecified || !outputDirectorySpecified)
	{
		cout << desc << endl;
		cout << "  --tumorDirectoryPath and --outputDirectoryPath must be specified" << endl << endl;
		std::cout << "  Usage: ./sinvict -e [errorRate] -m [minDepth] -l [leftStrandBias] -r [rightStrandBias] ";
		std::cout << "-f [readEndFraction] -q [QScoreCutoff] -t <tumorDirectoryPath> -o <outputDirectoryPath>" << std::endl;
		exit( 1);
	}
	else
	{
		if(vm.count("errorRate"))
		{
			errorRate = vm["errorRate"].as<double>();
			cout << "--errorRate set to: " << errorRate << endl;
		}
		else
		{
			cout << "--errorRate not set, using default value: " << errorRate << endl;
		}

		if(vm.count("minDepth"))
		{
			minDepth = vm["minDepth"].as<int>();
			cout << "--minDepth set to: " << minDepth << endl;
		}
		else
		{
			cout << "--minDepth not set, using default value: " << minDepth << endl;
		}

		if(vm.count("leftStrandBias"))
		{
			leftStrandBias = vm["leftStrandBias"].as<double>();
			cout << "--leftStrandBias set to: " << leftStrandBias << endl;
		}
		else
		{
			cout << "--leftStrandBias not set, using default value: " << leftStrandBias << endl;
		}

		if(vm.count("rightStrandBias"))
		{
			rightStrandBias = vm["rightStrandBias"].as<double>();
			cout << "--rightStrandBias set to: " << rightStrandBias << endl;
		}
		else
		{
			cout << "--rightStrandBias not set, using default value: " << rightStrandBias << endl;
		}

		if(vm.count("readEndFraction"))
		{
			readEndFraction = vm["readEndFraction"].as<double>();
			cout << "--readEndFraction set to: " << readEndFraction << endl;
		}
		else
		{
			cout << "--readEndFraction not set, using default value: " << readEndFraction << endl;
		}

		if(vm.count("QScoreCutoff"))
		{
			qScoreCutoff = vm["QScoreCutoff"].as<int>();
			cout << "--QScoreCutoff set to: " << qScoreCutoff << endl;
		}
		else
		{
			cout << "--QScoreCutoff not set, using default value: " << qScoreCutoff << endl;
		}

		Caller caller( errorRate, minDepth, leftStrandBias, rightStrandBias, readEndFraction, qScoreCutoff, tumorDirectoryPath.c_str(), outputDirectoryPath.c_str());
		caller.callLocationsMixture();
	}

	return 0;
}
