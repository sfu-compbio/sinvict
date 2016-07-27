#include <iostream>
#include <dirent.h>
#include <cerrno>
#include <algorithm>
#include "Common.h"

std::vector<std::string> Common::split( const std::string &str, const std::string &delimiter, bool keepEmpty)
{
	std::vector<std::string> result;
	if( delimiter.empty())
	{
		result.push_back( str);
		return result;
	}

	std::string::const_iterator substart = str.begin(), subend;
	while( true)
	{
		subend = std::search( substart, str.end(), delimiter.begin(), delimiter.end());
		std::string temp( substart, subend);
		if( keepEmpty || !temp.empty())
		{
			result.push_back( temp);
		}

		if( subend == str.end())
		{
			break;
		}
		substart = subend + delimiter.size();
	}
	
	return result;
}

std::string Common::revComp( const std::string &str)
{
	std::string reverseComplement( str);
	std::reverse( reverseComplement.begin(), reverseComplement.end());
	int i;

	for( i = 0; i < reverseComplement.size(); i++)
	{
		if( reverseComplement[i] == 'A')
		{
			reverseComplement[i] = 'T';
		}
		else if( reverseComplement[i] == 'C')
		{
			reverseComplement[i] = 'G';
		}
		else if( reverseComplement[i] == 'G')
		{
			reverseComplement[i] = 'C';
		}
		else if( reverseComplement[i] == 'T')
		{
			reverseComplement[i] = 'A';
		}
	}

	return reverseComplement;
}

std::map<std::string, std::string> Common::loadDNACodonTable()
{
	std::map<std::string, std::string> dnaCodonTable;

	// NOTE: The first M should be interpreted as START instead
	dnaCodonTable["ATG"] = "M";
	dnaCodonTable["GCT"] = "A"; dnaCodonTable["GCC"] = "A"; dnaCodonTable["GCA"] = "A"; dnaCodonTable["GCG"] = "A";
	dnaCodonTable["CGT"] = "R"; dnaCodonTable["CGC"] = "R"; dnaCodonTable["CGA"] = "R"; dnaCodonTable["CGG"] = "R";
	dnaCodonTable["AGA"] = "R"; dnaCodonTable["AGG"] = "R"; dnaCodonTable["AAT"] = "N"; dnaCodonTable["AAC"] = "N";
	dnaCodonTable["GAT"] = "D"; dnaCodonTable["GAC"] = "D"; dnaCodonTable["TGT"] = "C"; dnaCodonTable["TGC"] = "C";
	dnaCodonTable["CAA"] = "Q"; dnaCodonTable["CAG"] = "Q"; dnaCodonTable["GAA"] = "E"; dnaCodonTable["GAG"] = "E";
	dnaCodonTable["GGT"] = "G"; dnaCodonTable["GGC"] = "G"; dnaCodonTable["GGA"] = "G"; dnaCodonTable["GGG"] = "G";
	dnaCodonTable["CAT"] = "H"; dnaCodonTable["CAC"] = "H"; dnaCodonTable["ATT"] = "I"; dnaCodonTable["ATC"] = "I";
	dnaCodonTable["ATA"] = "I"; dnaCodonTable["TTA"] = "L"; dnaCodonTable["TTG"] = "L"; dnaCodonTable["CTT"] = "L";
	dnaCodonTable["CTC"] = "L"; dnaCodonTable["CTA"] = "L"; dnaCodonTable["CTG"] = "L"; dnaCodonTable["AAA"] = "K";
	dnaCodonTable["AAG"] = "K"; dnaCodonTable["TTT"] = "F"; dnaCodonTable["TTC"] = "F"; dnaCodonTable["CCT"] = "P";
	dnaCodonTable["CCC"] = "P"; dnaCodonTable["CCA"] = "P"; dnaCodonTable["CCG"] = "P"; dnaCodonTable["TCT"] = "S";
	dnaCodonTable["TCC"] = "S"; dnaCodonTable["TCA"] = "S"; dnaCodonTable["TCG"] = "S"; dnaCodonTable["AGT"] = "S";
	dnaCodonTable["AGC"] = "S"; dnaCodonTable["ACT"] = "T"; dnaCodonTable["ACC"] = "T"; dnaCodonTable["ACA"] = "T";
	dnaCodonTable["ACG"] = "T"; dnaCodonTable["TGG"] = "W"; dnaCodonTable["TAT"] = "Y"; dnaCodonTable["TAC"] = "Y";
	dnaCodonTable["GTT"] = "V"; dnaCodonTable["GTC"] = "V"; dnaCodonTable["GTA"] = "V"; dnaCodonTable["GTG"] = "V";
	dnaCodonTable["TAA"] = "STOP"; dnaCodonTable["TGA"] = "STOP"; dnaCodonTable["TAG"] = "STOP";

	return dnaCodonTable;
}

std::string Common::getAminoAcids( const std::string &str, const std::map<std::string, std::string> &table)
{
	std::string aminoAcids;
	int i;
	int j;
	int limit;
	
	limit = str.size() / 3;
	aminoAcids.resize( limit + 1);

	limit = limit * 3;
	
	j = 0;
	for( i = 0; i < limit; i = i + 3)
	{
		std::string key = str.substr( i, 3);
		std::map<std::string, std::string>::const_iterator it = table.find( key);
		std::string aaChar = it->second;
		aminoAcids[j] = ( aaChar.c_str())[0];
		j++;
	}
	aminoAcids[j] = '\0';

	return aminoAcids;
}

int Common::max3( const int a, const int b, const int c)
{
	if( a >= b && a >= c)
	{
		return a;
	}
	else if( b >= a && b >= c)
	{
		return b;
	}
	else if( c >= a && c >= b)
	{
		return c;
	}
}

void Common::getFilesInDir( const std::string& path, std::vector<std::string>& fileNames)
{
	DIR* dir;
	struct dirent* dirEntry;
	std::string nextFilename;

	dir = opendir( path.c_str());
	if( dir != NULL)
	{
		dirEntry = readdir( dir);
		while( dirEntry != NULL)
		{
			nextFilename.assign( dirEntry->d_name);
			if( nextFilename != "." && nextFilename != "..")
			{
				fileNames.push_back( nextFilename);
			}
			dirEntry = readdir( dir);
		}
		closedir( dir);
	}
	else
	{
		perror( "Error opening directory");
		exit( 1);
	}
}
