#ifndef ALLELE_H
#define ALLELE_H

#include <fstream>
#include <string>
#include <vector>

class Allele
{
	private:
		std::string base;
		int count;
		double avgMappingQuality;
		double avgBaseQuality;
		double avgSEMappingQuality;
		int numPlusStrand;
		int numMinusStrand;
		double avgPosAsFraction;
		double avgNumMismatchesAsFraction;
		double avgSumMismatchQualities;
		int numQ2ContainingReads;
		double avgDistanceToQ2StartInQ2Reads;
		double avgClippedLength;
		double avgDistanceToEffective3pEnd;

		// Inferred variables
		bool variant;
		double percentage;
		double pValue;
		double qScore;

	public:
		Allele( const std::string base, const int count, const double avgMappingQuality,
				const double avgBaseQuality, const double avgSEMappingQuality,
				const int numPlusStrand, const int numMinusStrand, const double avgPosAsFraction,
				const double avgNumMismatchesAsFraction, const double avgSumMismatchQualities,
				const int numQ2ContainingReads, const double avgDistanceToQ2StartInQ2Reads,
				const double avgClippedLength, const double avgDistanceToEffective3pEnd,
				const double percentage, const bool variant);
		std::string getBase();
		int getCount();
		double getAvgMappingQuality();
		double getAvgBaseQuality();
		double getAvgSEMappingQuality();
		int getNumPlusStrand();
		int getNumMinusStrand();
		double getAvgPosAsFraction();
		double getAvgNumMismatchesAsFraction();
		double getAvgSumMismatchQualities();
		int getNumQ2ContainingReads();
		double getAvgDistanceToQ2StartInQ2Reads();
		double getAvgClippedLength();
		double getAvgDistanceToEffective3pEnd();

		// Methods for inferred variables
		bool isVariant();
		double getPercentage();
		double getPValue();
		double getQScore();
		void setPValue( double pValue);
		void setQScore( double qScore);
		void setMostFreqVariant( bool mostFreq);

		// Print allele information
		void printAllele();
		void printAllele( std::ofstream& out, int readDepth, int usePoissonGermline);
		void printAlleleUCSC( std::ofstream& out);
};
#endif
