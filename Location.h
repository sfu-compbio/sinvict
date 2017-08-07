#ifndef LOCATION_H
#define LOCATION_H

#include <fstream>
#include <string>
#include <vector>
#include "Sample.h"

class Location
{
	private:
		std::string chr;
		int position;
		std::vector<Sample> samples;
		std::string mutatedBase;

		// The mean, variance, standard deviation, and Coefficient of Variation for the Variant Allele Percentages
		double meanVAP;
		double varianceVAP;
		double stdVAP;
		double cov;

	public:
		Location( const std::string chr, const int position);
		std::string getChr();
		int getPosition();
		std::vector<Sample> getSamples();
		void addSample( Sample newSample);
		void clearSamples();
		double getMeanVAP();
		double getVarianceVAP();
		double getStdVAP();
		double getCOV();
		void setMeanVAP( double newMean);
		void setVarianceVAP( double newVariance);
		void setStdVAP( double newStd);
		void setCOV( double newCOV);
		void setMutatedBase( std::string base);
		void printLocation();
		void printLocation( std::ofstream& out, int usePoissonGermline);
		void printVCFHeader( std::ofstream& out);
		void printLocationVCF( std::ofstream& out);
		void printUCSCHeader( std::ofstream& out);
		void printLocationUCSC( std::ofstream& out);
		void printLocationANNOVAR();
		void printLocationANNOVAR( std::ofstream& out);
		void printLocationCITUP( std::ofstream& out);
};
#endif
