#include "FASTQParser.h"

FASTQParser::FASTQParser( const std::string& path)
{
	inputFile.open( path, std::fstream::in);

	// Prepare the output path
	int pathEnd = path.find_last_of( ".");
	std::string basePath = path.substr( 0, pathEnd);
	std::string outputPath = basePath + "_trimmed.fastq";

	outputFile.open( outputPath, std::fstream::out);
}

FASTQParser::~FASTQParser()
{
	inputFile.close();
	outputFile.close();
}

void FASTQParser::trimReadEnds( std::fstream& in, std::fstream& out, int trimLength)
{
	std::string readName;
	std::string seq;
	std::string optional;
	std::string qual;
	std::string newSeq;
	std::string newQual;
	int seqSize;
	int totalTrimSize;

	int totalNumberOfLines = 0;
	int lineCount = 1;
	std::string line;
	while( std::getline( in, line))
	{
		if( lineCount == 1)
		{
			readName.assign( line);
		}
		else if( lineCount == 2)
		{
			seqSize = seq.size();
			totalTrimSize = seqSize - ( 2 * trimLength);

			seq.assign( line);
			newSeq = seq.substr( trimLength, totalTrimSize);
		}
		else if( lineCount == 3)
		{
			optional.assign( line);
		}
		else if( lineCount == 4)
		{
			qual.assign( line);
			newQual = qual.substr( trimLength, totalTrimSize);
			
			// Write out the next trimmed entry
			out << readName << '\n' << newSeq << '\n' << optional << '\n' << newQual << '\n';
		}

		// Increment the total number of lines parsed
		totalNumberOfLines++;

		// Increment line count, reset the count if all four lines of an entry have been parsed
		lineCount++;
		if( lineCount == 5)
		{
			lineCount = 1;
		}
	
	}

	if( totalNumberOfLines % 4 != 0)
	{
		std::cerr << "Number of lines not a multiple of 4. Check FASTQ file correctness." << std::endl;
		exit( 1);
	}
}
