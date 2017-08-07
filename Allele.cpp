#include <iostream>
#include "Allele.h"
#include "Filter.h"

Allele::Allele( const std::string base, const int count, const double avgMappingQuality,
		const double avgBaseQuality, const double avgSEMappingQuality,
		const int numPlusStrand, const int numMinusStrand, const double avgPosAsFraction,
		const double avgNumMismatchesAsFraction, const double avgSumMismatchQualities,
		const int numQ2ContainingReads, const double avgDistanceToQ2StartInQ2Reads,
		const double avgClippedLength, const double avgDistanceToEffective3pEnd,
		const double percentage, const bool variant)
{
	this->base = base;
	this->count = count;
	this->avgMappingQuality = avgMappingQuality;
	this->avgBaseQuality = avgBaseQuality;
	this->avgSEMappingQuality = avgSEMappingQuality;
	this->numPlusStrand = numPlusStrand;
	this->numMinusStrand = numMinusStrand;
	this->avgPosAsFraction = avgPosAsFraction;
	this->avgNumMismatchesAsFraction = avgNumMismatchesAsFraction;
	this->avgSumMismatchQualities = avgSumMismatchQualities;
	this->numQ2ContainingReads = numQ2ContainingReads;
	this->avgDistanceToQ2StartInQ2Reads = avgDistanceToQ2StartInQ2Reads;
	this->avgClippedLength = avgClippedLength;
	this->avgDistanceToEffective3pEnd = avgDistanceToEffective3pEnd;

	// Inferred Variables
	this->variant = variant;
	this->percentage = percentage;
	this->pValue = -1;
	this->qScore = -1;
}

std::string Allele::getBase()
{
	return base;
}

int Allele::getCount()
{
	return count;
}

double Allele::getAvgMappingQuality()
{
	return avgMappingQuality;
}

double Allele::getAvgBaseQuality()
{
	return avgBaseQuality;
}

double Allele::getAvgSEMappingQuality()
{
	return avgSEMappingQuality;
}

int Allele::getNumPlusStrand()
{
	return numPlusStrand;
}

int Allele::getNumMinusStrand()
{
	return numMinusStrand;
}

double Allele::getAvgPosAsFraction()
{
	return avgPosAsFraction;
}

double Allele::getAvgNumMismatchesAsFraction()
{
	return avgNumMismatchesAsFraction;
}

double Allele::getAvgSumMismatchQualities()
{
	return avgSumMismatchQualities;
}

int Allele::getNumQ2ContainingReads()
{
	return numQ2ContainingReads;
}

double Allele::getAvgDistanceToQ2StartInQ2Reads()
{
	return avgDistanceToQ2StartInQ2Reads;
}

double Allele::getAvgClippedLength()
{
	return avgClippedLength;
}

double Allele::getAvgDistanceToEffective3pEnd()
{
	return avgDistanceToEffective3pEnd;
}

bool Allele::isVariant()
{
	return variant;
}

double Allele::getPercentage()
{
	return percentage;
}

double Allele::getPValue()
{
	return pValue;
}

double Allele::getQScore()
{
	return qScore;
}

void Allele::setPValue( const double pValue)
{
	this->pValue = pValue;
}

void Allele::setQScore( const double qScore)
{
	this->qScore = qScore;
}

void Allele::printAllele()
{
	std::cout << base << "\t" << count << "\t" << "isVar:" << variant << "\t" << percentage << "\t" << "+:" <<  numPlusStrand << "\t" << "-:" << numMinusStrand << "\t";
	std::cout << avgPosAsFraction;
}

void Allele::printAllele( std::ofstream& out, int readDepth, int usePoissonGermline)
{
	out << base << "\t" << count << "\t" << percentage << "\t" << "+:" << numPlusStrand << "\t" << "-:" << numMinusStrand << "\t";
	out << avgPosAsFraction << "\t";

	if(usePoissonGermline)
	{
		double lambda2 = readDepth * 0.5;
		double pval2 = Filter::illuminaPoissonFilter(count, lambda2);
		if(pval2 < 0.05)
		{
			out << "Germline";
		}
		else
		{
			out << "Somatic";
		}
	}
	else
	{
		if(percentage >= 40)
		{
			out << "Germline";
		}
		else
		{
			out << "Somatic";
		}
	}
}

void Allele::printAlleleUCSC( std::ofstream& out)
{
	out << base << "\t" << variant << "\t" << percentage;
}
