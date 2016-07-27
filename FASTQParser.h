#ifndef FASTQ_PARSER_H
#define FASTQ_PARSER_H

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

class FASTQParser
{
	private:
		std::fstream inputFile;
		std::fstream outputFile;

	public:
		FASTQParser( const std::string& path);
		~FASTQParser();
		void trimReadEnds( std::fstream& inputFile, std::fstream& outputFile, int trimLength);
};
#endif
