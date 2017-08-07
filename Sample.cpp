#include <iostream>
#include "Sample.h"

Sample::Sample( const std::string path, const ReadcountEntry& readcountEntry)
{
	this->path = path;
	
	// Strip the path from the filename to form the sample name
	int sampleNameBegin = path.find_last_of( "/");
	int sampleNameEnd = path.find_last_of( ".");
	int sampleNameLength = sampleNameEnd - sampleNameBegin - 1;
	sampleName = path.substr( sampleNameBegin + 1, sampleNameLength);

	this->readcountEntry = readcountEntry;
}

ReadcountEntry Sample::getReadcountEntry()
{
	return readcountEntry;
}

std::string Sample::getSampleName()
{
	return sampleName;
}

std::string Sample::getPath()
{
	return path;
}

void Sample::printSample()
{
	std::cout << sampleName << "\t";
	readcountEntry.printReadcountEntry();
}

void Sample::printSample( std::ofstream& out, int usePoissonGermline)
{
	out << sampleName << "\t";
	readcountEntry.printReadcountEntry( out, usePoissonGermline);
}

void Sample::printSampleUCSC( std::ofstream& out, int usePoissonGermline)
{
	out << sampleName << "\t";
	readcountEntry.printReadcountEntry( out, usePoissonGermline);
}
