#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <getopt.h>
#include "Caller.h"

void printHelp();

int main(int argc, char** argv)
{
	int opt, opt_index;
	
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
	int usePoissonGermline = 1;
	int disableLvl5Filter = 0;

	static struct option long_options[] = {
		{"error-rate",				required_argument, 0,  'e' },
		{"min-depth",				required_argument, 0,  'm' },
		{"left-strand-bias",		required_argument, 0,  'l' },
		{"right-strand-bias",		required_argument, 0,  'r' },
		{"read-end-fraction",		required_argument, 0,  'f' },
		{"qscore-cutoff",			required_argument, 0,  'q' },
		{"tumor-directory-path",	required_argument, 0,  't' },
		{"output-directory-path",	required_argument, 0,  'o' },
		{"use-poisson-germline",	required_argument, 0,  's' },
		{"disable-lvl5-filter",		required_argument, 0,  '5' },
		{0, 0, 0, 0 }
	};

	while ( -1 !=  (opt = getopt_long( argc, argv, "he:m:l:r:f:q:t:o:s:5:", long_options, &opt_index )  ) )
	{
		switch(opt)
		{
			case 'h':
				printHelp();
				return 0;
			case 'e':
				errorRate 		= atof( optarg );
				break;
			case 'm':
				minDepth  		= atoi( optarg );
				break;
			case 'l':
				leftStrandBias  = atof( optarg );
				break;
			case 'r':
				rightStrandBias	= atof( optarg );
				break;
			case 'f':
				readEndFraction	= atof( optarg );
				break;
			case 'q':
				qScoreCutoff  	= atoi( optarg );
				break;
			case 't':
				tumorDirectoryPath = std::string( optarg );
				tumorDirectorySpecified = true;
				break;
			case 'o':
				outputDirectoryPath = std::string( optarg );
				outputDirectorySpecified = true;
				break;
			case 's':
				usePoissonGermline = atoi(optarg);
				break;
			case '5':
				disableLvl5Filter = atoi(optarg);
				break;
			default:
				printHelp();
				return 0;
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
	std::cout << "--use-poisson-germline specified with value = " << usePoissonGermline << std::endl;
	std::cout << "--disable-lvl5-filter specified with value = " << disableLvl5Filter << std::endl;
	std::cout << "--tumor-directory-path specified with value = " << tumorDirectoryPath << std::endl;
	std::cout << "--output-directory-path specified with value = " << outputDirectoryPath << std::endl;

	Caller caller(errorRate, minDepth, leftStrandBias, rightStrandBias, readEndFraction, qScoreCutoff, tumorDirectoryPath.c_str(), outputDirectoryPath.c_str(), usePoissonGermline, disableLvl5Filter);
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
	desc = desc + "\t--use-poisson-germline or -s : Use a more robust poisson model to guess somatic/germline status.\n";
	desc = desc + "\t--disable-lvl5-filter or -5 : If the input value is 1, SNR filter is disabled.\n";
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
	std::cout << "-f=[read-end-fraction] -q=[qscore-cutoff] -s=<use-poisson-germline> -5=<disable-lvl5-filter> -t=<tumor-directory-path> -o=<output-directory-path> \n" << std::endl;
}
