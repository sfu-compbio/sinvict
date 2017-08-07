#include <iostream>
#include "ReadcountEntry.h"

ReadcountEntry::ReadcountEntry( const std::string refBase, const int readDepth)
{
	this->refBase = refBase;
	this->readDepth = readDepth;
	this->indexMostFreqVariantAllele = -1;
}

ReadcountEntry::ReadcountEntry()
{
	this->refBase = "E";
	this->readDepth = -1;
	this->indexMostFreqVariantAllele = -1;
}

std::string ReadcountEntry::getRefBase()
{
	return refBase;
}

int ReadcountEntry::getReadDepth()
{
	return readDepth;
}

std::vector<Allele> ReadcountEntry::getAlleles()
{
	return alleles;
}

void ReadcountEntry::addAllele( Allele newAllele)
{
	alleles.push_back( newAllele);
}

void ReadcountEntry::setMostFreqVariantAllele()
{
	double mostFreqVariantPercentage = -1;
	int numberOfAlleles = alleles.size();
	int index = -1;
	
	int i;
	for( i = 0; i < numberOfAlleles; i++)
	{
		if( alleles[i].isVariant() && alleles[i].getPercentage() > mostFreqVariantPercentage)
		{
			mostFreqVariantPercentage = alleles[i].getPercentage();
			index = i;
		}
	}

	indexMostFreqVariantAllele = index;
}

Allele ReadcountEntry::getMostFreqVariantAllele()
{
	return alleles[indexMostFreqVariantAllele];
}

int ReadcountEntry::getIndexMostFreqVariantAllele()
{
	return indexMostFreqVariantAllele;
}

void ReadcountEntry::printReadcountEntry()
{
	std::cout << refBase << "\t" << readDepth;
	for( int i = 0; i < alleles.size(); i++)
	{
		std::cout << "\n" << "\t" << "\t" << "\t";
		alleles[i].printAllele();
		if( i == indexMostFreqVariantAllele)
		{
			std::cout << "\t" << "MFVA";
		}
	}
}

void ReadcountEntry::printReadcountEntry( std::ofstream& out, int usePoissonGermline)
{
	out << refBase << "\t" << readDepth << "\t";
	for( int i = 0; i < alleles.size(); i++)
	{
		if( i == indexMostFreqVariantAllele)
		{
			alleles[i].printAllele( out, readDepth, usePoissonGermline);
		}
	}
}

void ReadcountEntry::printReadcountEntryUCSC( std::ofstream& out)
{
	out << refBase << "\t" << readDepth;
	for( int i = 0; i < alleles.size(); i++)
	{
		if( i == indexMostFreqVariantAllele)
		{
			alleles[i].printAlleleUCSC( out);
		}
	}
}
