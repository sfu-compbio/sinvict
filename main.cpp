#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "argh.h"
#include "Caller.h"

void printHelp();

int main(int argc, char** argv)
{
	// Initialize the "argh" command-line parser
	argh::parser cmdl(argc, argv);

	// If the help flag is set, print help message and exit
	if(cmdl["--help"] || cmdl["-h"])
	{
		printHelp();
		return 0;
	}

	// Default parameters
	double errorRate = 0.01;
	int minDepth = 100;
	double leftStrandBias = 0.3;
	double rightStrandBias = 0.7;
	double readEndFraction = 0.01;
	int qScoreCutoff = 95;
	std::string tumorDirectoryPath;
	std::string outputDirectoryPath;
	bool tumorDirectorySpecified = false;
	bool outputDirectorySpecified = false;

	// Now go through all command arguments to overwrite default arguments when applicable
	for(auto& param : cmdl.params())
	{
		if(param.first == "tumor-directory-path" || param.first == "t")
		{
			if(!(cmdl(param.first) >> tumorDirectoryPath))
			{
				std::cerr << "Invalid value entered for tumor-directory-path. Must be string." << std::endl;
				printHelp();
				return 0;
			}
			tumorDirectorySpecified = true;
		}
		else if(param.first == "output-directory-path" || param.first == "o")
		{
			if(!(cmdl(param.first) >> outputDirectoryPath))
			{
				std::cerr << "Invalid value entered for output-directory-path. Must be string." << std::endl;
				printHelp();
				return 0;
			}
			outputDirectorySpecified = true;
		}
		else if(param.first == "error-rate" || param.first == "e")
		{
			if(!(cmdl(param.first) >> errorRate))
			{
				std::cerr << "Invalid value entered for error-rate. Must be double >= 0.0" << std::endl;
				printHelp();
				return 0;
			}
		}
		else if(param.first == "min-depth" || param.first == "m")
		{
			if(!(cmdl(param.first) >> minDepth))
			{
				std::cerr << "Invalid value entered for min-depth. Must be integer > 0." << std::endl;
				printHelp();
				return 0;
			}
		}
		else if(param.first == "left-strand-bias" || param.first == "l")
		{
			if(!(cmdl(param.first) >> leftStrandBias))
			{
				std::cerr << "Invalid value entered for left-strand-bias. Must be double >= 0.0" << std::endl;
				printHelp();
				return 0;
			}
		}
		else if(param.first == "right-strand-bias" || param.first == "r")
		{
			if(!(cmdl(param.first) >> rightStrandBias))
			{
				std::cerr << "Invalid value entered for right-strand-bias. Must be double >= 0.0" << std::endl;
				printHelp();
				return 0;
			}
		}
		else if(param.first == "read-end-fraction" || param.first == "f")
		{
			if(!(cmdl(param.first) >> readEndFraction))
			{
				std::cerr << "Invalid value entered for read-end-fraction. Must be double >= 0.0" << std::endl;
				printHelp();
				return 0;
			}
		}
		else if(param.first == "qscore-cutoff" || param.first == "q")
		{
			if(!(cmdl(param.first) >> qScoreCutoff))
			{
				std::cerr << "Invalid value entered for qscore-cutoff. Must be integer > 0" << std::endl;
				printHelp();
				return 0;
			}
		}
	}

	// Input and Output directories MUST be specified
	// If not, print help message and exit
	if(!tumorDirectorySpecified || !outputDirectorySpecified)
	{
		printHelp();
		return 0;
	}

	// Print final set of parameters
	std::cout << "--error-rate set to : " << errorRate << std::endl;
	std::cout << "--min-depth set to : " << minDepth << std::endl;
	std::cout << "--left-strand-bias set to : " << leftStrandBias << std::endl;
	std::cout << "--right-strand-bias set to: " << rightStrandBias << std::endl;
	std::cout << "--read-end-fraction set to: " << readEndFraction << std::endl;
	std::cout << "--qscore-cutoff set to: " << qScoreCutoff << std::endl;
	std::cout << "--tumor-directory-path specified with value = " << tumorDirectoryPath << std::endl;
	std::cout << "--output-directory-path specified with value = " << outputDirectoryPath << std::endl;

	Caller caller(errorRate, minDepth, leftStrandBias, rightStrandBias, readEndFraction, qScoreCutoff, tumorDirectoryPath.c_str(), outputDirectoryPath.c_str());
	caller.callLocationsMixture();

	return 0;
}

void printHelp()
{
	// Help message body
	std::string desc = "";
	desc = desc + "\t--help or -h : Print help message.\n";
	desc = desc + "\t--error-rate or -e : Error rate for the sequencing platform used.\n";
	desc = desc + "\t--min-depth or -m : Minimum Read Depth required for high confidence in a call.\n";
	desc = desc + "\t--left-strand-bias or -l : Lower limit for the strand bias value interval to be used in assessing the confidence of a call.\n";
	desc = desc + "\t--right-strand-bias or -r : Upper limit for the strand bias value interval to be used in assessing the confidence of a call.\n";
	desc = desc + "\t--read-end-fraction or -f : Average position of the called base on the reads supporting the call as a fraction.";
	desc = desc + "End values such as 0.01 as useful for filtering read end artifacts.\n";
	desc = desc + "\t--qscore-cutoff or -q : Cutoff value for the qScore assigned to each call by the Poisson model used.\n";
	desc = desc + "\t--tumor-directory-path or -t : Specifies directory for the input files.\n";
	desc = desc + "\t--output-directory-path or -o : Specifies directory for the output files.\n";

	// Print help message
	std::cout << ("\nSiNVICT: Ultra Sensitive Detection of Single Nucleotide Variants and Indels in Circulating Tumour DNA.");
	std::cout << std::endl << std::endl;
	std::cout << ("Allowed arguments:");
	std::cout << std::endl;
	std::cout << desc;
	std::cout << std::endl;

	std::cout << "\t--tumor-directory-path and --output-directory-path must be specified\n" << std::endl;
	std::cout << "\tUsage: ./sinvict -e=[error-rate] -m=[min-depth] -l=[left-strand-bias] -r=[right-strand-bias] ";
	std::cout << "-f=[read-end-fraction] -q=[qscore-cutoff] -t=<tumor-directory-path> -o=<output-directory-path>\n" << std::endl;
}
