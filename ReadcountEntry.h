#ifndef READCOUNT_ENTRY_H
#define READCOUNT_ENTRY_H

#include <fstream>
#include <string>
#include <vector>
#include "Allele.h"

class ReadcountEntry
{
	private:
		std::string refBase;
		int readDepth;
		std::vector<Allele> alleles;
		int indexMostFreqVariantAllele;

	public:
		ReadcountEntry();
		ReadcountEntry( const std::string refBase, const int readDepth);
		std::string getRefBase();
		int getReadDepth();
		std::vector<Allele> getAlleles();
		void addAllele( Allele newAllele);
		void setMostFreqVariantAllele();
		Allele getMostFreqVariantAllele();
		int getIndexMostFreqVariantAllele();
		void printReadcountEntry();
		void printReadcountEntry( std::ofstream& out, int usePoissonGermline);
		void printReadcountEntryUCSC( std::ofstream& out);
};
#endif
