#ifndef SAMPLE_H
#define SAMPLE_H

#include <fstream>
#include <string>
#include "ReadcountEntry.h"

class Sample
{
	private:
		std::string path;
		std::string sampleName;
		ReadcountEntry readcountEntry;

	public:
		Sample( const std::string path, const ReadcountEntry& readcountEntry);
		ReadcountEntry getReadcountEntry();
		std::string getPath();
		std::string getSampleName();
		void printSample();
		void printSample( std::ofstream& out, int usePoissonGermline);
		void printSampleUCSC( std::ofstream& out, int usePoissonGermline);
};
#endif
