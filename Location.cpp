#include <iostream>
#include "Location.h"

Location::Location( const std::string chr, const int position)
{
	this->chr = chr;
	this->position = position;

	// Initialize location statistics
	meanVAP = -1.0;
	varianceVAP = -1.0;
	stdVAP  = -1.0;
	cov = -1.0;
}

std::string Location::getChr()
{
	return chr;
}

int Location::getPosition()
{
	return position;
}

std::vector<Sample> Location::getSamples()
{
	return samples;
}

void Location::addSample( Sample newSample)
{
	samples.push_back( newSample);
}

void Location::clearSamples()
{
	samples.clear();
}

double Location::getMeanVAP()
{
	return meanVAP;
}

double Location::getVarianceVAP()
{
	return varianceVAP;
}

double Location::getStdVAP()
{
	return stdVAP;
}

double Location::getCOV()
{
	return cov;
}

void Location::setMeanVAP( double newMean)
{
	this->meanVAP = newMean;
}

void Location::setVarianceVAP( double newVariance)
{
	this->varianceVAP = newVariance;
}

void Location::setStdVAP( double newStd)
{
	this->stdVAP = newStd;
}

void Location::setCOV( double newCOV)
{
	this->cov = newCOV;
}

void Location::setMutatedBase( std::string base)
{
	this->mutatedBase = base;
}

void Location::printLocation()
{
//	std::cout << chr << "\t" << position << "\t" << meanVAP << "\t" << varianceVAP << "\t" << stdVAP << "\t" << cov << "\n";
	
	std::cout << chr << "\t" << position << "\n";
	for( int i = 0; i < samples.size(); i++)
	{
		std::cout << "\t";
		samples[i].printSample();
		std::cout << "\n";
	}	
}

void Location::printLocation( std::ofstream& out, int usePoissonGermline)
{
/*	out << chr << "\t" << position << "\t" << meanVAP << "\t" << varianceVAP << "\t" << stdVAP << "\t" << cov << "\t";
	for( int i = 0; i < samples.size(); i++)
	{
		out << samples[i].getSampleName() << "\t";
	}
	out << "\n";
*/
	for( int i = 0; i < samples.size(); i++)
	{
		out << chr << "\t" << position << "\t";
		samples[i].printSample( out, usePoissonGermline);
		out << "\n";
	}
}

void Location::printVCFHeader( std::ofstream& out)
{
	// Meta-information lines
	out << "##fileformat=VCFv4.2" << "\n";
	out << "##source=sinvict" << "\n";

	// Header line
	out << "#CHROM" << "\t" << "POS" << "\t" << "ID" << "\t" << "REF" << "\t";
	out << "ALT" << "\t" << "QUAL" << "\t" << "FILTER" << "\t" << "FORMAT" << "\n";
}

void Location::printLocationVCF( std::ofstream& out)
{
	ReadcountEntry currentEntry = samples[0].getReadcountEntry();
	Allele variantAllele = currentEntry.getMostFreqVariantAllele();

	std::string refBase = currentEntry.getRefBase();
	std::string altBase = variantAllele.getBase();
	if( altBase == "A" || altBase == "C" || altBase == "G" || altBase == "T" || altBase == "N")
	{
		out << chr << "\t" << position << "\t" << "." << "\t";
		out << refBase << "\t" << altBase << "." << "\t" << "PASS" << "\t" << "." << "\n";
	}
}

void Location::printUCSCHeader( std::ofstream& out)
{
	out << "UCSC_Coordinate" << "\t" << "Chr" << "\t" << "Pos" << "\t" << "RefBase" << "\t" << "AltBase" << "\t";
	for( int i = 0; i < samples.size(); i++)
	{
		out << samples[i].getSampleName() << "\t";
	}

	out << "Mean" << "\t" << "Var" << "\t" << "StdDev" << "\t" << "CoV" << "\n";
}

void Location::printLocationUCSC( std::ofstream& out)
{
	ReadcountEntry currentEntry = samples[0].getReadcountEntry();

	// A
	out << chr << ":" << position << "\t" << chr << "\t" << position << "\t";
	out << currentEntry.getRefBase() << "\t" << "A" << "\t";

	for( int i = 0; i < samples.size(); i++)
	{
		currentEntry = samples[i].getReadcountEntry();
		std::vector<Allele> alleles = currentEntry.getAlleles();

		out << alleles[0].getPercentage() << "\t";
	}

	out << meanVAP << "\t" << varianceVAP << "\t" << stdVAP << "\t" << cov << "\n";

	// C
	out << chr << ":" << position << "\t" << chr << "\t" << position << "\t";
	out << currentEntry.getRefBase() << "\t" << "C" << "\t";

	for( int i = 0; i < samples.size(); i++)
	{
		currentEntry = samples[i].getReadcountEntry();
		std::vector<Allele> alleles = currentEntry.getAlleles();

		out << alleles[1].getPercentage() << "\t";
	}

	out << meanVAP << "\t" << varianceVAP << "\t" << stdVAP << "\t" << cov << "\n";

	// G
	out << chr << ":" << position << "\t" << chr << "\t" << position << "\t";
	out << currentEntry.getRefBase() << "\t" << "G" << "\t";

	for( int i = 0; i < samples.size(); i++)
	{
		currentEntry = samples[i].getReadcountEntry();
		std::vector<Allele> alleles = currentEntry.getAlleles();

		out << alleles[2].getPercentage() << "\t";
	}

	out << meanVAP << "\t" << varianceVAP << "\t" << stdVAP << "\t" << cov << "\n";

	// T
	out << chr << ":" << position << "\t" << chr << "\t" << position << "\t";
	out << currentEntry.getRefBase() << "\t" << "T" << "\t";

	for( int i = 0; i < samples.size(); i++)
	{
		currentEntry = samples[i].getReadcountEntry();
		std::vector<Allele> alleles = currentEntry.getAlleles();

		out << alleles[3].getPercentage() << "\t";
	}

	out << meanVAP << "\t" << varianceVAP << "\t" << stdVAP << "\t" << cov << "\n";
}

void Location::printLocationANNOVAR( std::ofstream& out)
{
	ReadcountEntry currentEntry = samples[0].getReadcountEntry();
	Allele variantAllele = currentEntry.getMostFreqVariantAllele();

	std::string altBase = variantAllele.getBase();
	if( altBase == "A" || altBase == "C" || altBase == "G" || altBase == "T")
	{
		out << chr << "\t" << position << "\t" << position << "\t";
		out << '0' << "\t" << altBase << "\n";
	}
}

void Location::printLocationCITUP( std::ofstream& out)
{
	ReadcountEntry currentEntry = samples[0].getReadcountEntry();
	Allele variantAllele = currentEntry.getMostFreqVariantAllele();

	std::string altBase = variantAllele.getBase();
	if( altBase == "A" || altBase == "C" || altBase == "G" || altBase == "T")
	{
		out << chr << "\t" << position << "\t";
		for( int i = 0; i < samples.size(); i++)
		{
			currentEntry = samples[i].getReadcountEntry();
			variantAllele = currentEntry.getMostFreqVariantAllele();

			// Do not multiply by 2, since we don't simulate diploid genome
			double frequency = ( double)  ( variantAllele.getCount()) / ( double) currentEntry.getReadDepth();
			out << frequency << "\t";
		}
		out << "\n";
	}
}
