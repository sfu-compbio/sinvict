#ifndef CALLER_H
#define CALLER_H

#include <string>
#include <unordered_map>
#include "Location.h"

class Caller
{
	private:
		// Hash Table for Location objects
		std::unordered_map<std::string, Location> locationTable;

		// Paths to the top level directories of the input and output readcount files
		std::string tumorDirectoryPath;
		std::string benignDirectoryPath;
		std::string outputDirectoryPath;

		// Paths to individual input and output files
		std::vector<std::string> tumorPaths;
		std::vector<std::string> benignPaths;
		std::vector<std::string> outputPaths;

		// Filter parameters
		double poissonLambda;
		int minDepth;
		double strandBiasLeft;
		double strandBiasRight;
		int minQScore;
		double readEndFraction;
		int usePoissonGermline;

		// Enable disabling certain filters
		int disableLvl5Filter;

	public:
		Caller( const double poissonLambda, const int minDepth, const double leftStrandBias, const double rightStrandBias, const double readEndFraction, const int qCutoff, const char* tumorDirectoryPath, const char* benignDirectoryPath, const char* outputDirectoryPath, const int usePoissonGermline, const int disableLvl5Filter);
		Caller( const double poissonLambda, const int minDepth, const double leftStrandBias, const double rightStrandBias, const double readEndFraction, const int qCutoff, const char* tumorDirectoryPath, const char* outputDirectoryPath, const int usePoissonGermline, const int disableLvl5Filter);
		int loadEntries( const std::string path);
		void calculateStatistics();
		int callLocationsMixture();
		std::vector<Location> callPoissonDist( double poissonLambda, int minQScore);
		std::vector<Location> callDepthFilter( std::vector<Location> unfilteredCalls, int minDepth);
		std::vector<Location> callStrandBiasFilter( std::vector<Location> unfilteredCalls, double strandBiasLeft, double strandBiasRight);
		std::vector<Location> callAmpliconEndFilter( std::vector<Location> unfilteredCalls, double readEndFraction);
		std::vector<Location> callAverageFilter( std::vector<Location> unfilteredCalls);
		std::vector<Location> callHomopolymerFilter( std::vector<Location> unfilteredCalls);
		void printCaller();
		void printLocations();
		void printUCSC( std::vector<Location> locations, std::ofstream& outputFile);
		void printCITUP( std::vector<Location> locations, std::ofstream& outputFile);
};
#endif
