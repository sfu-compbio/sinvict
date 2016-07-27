#ifndef COMMON_H
#define COMMON_H

#include <string>
#include <vector>
#include <map>

class Common
{
	public:
		static std::vector<std::string> split( const std::string &str, const std::string &delimiter, bool keepEmpty);
		static std::string revComp( const std::string &str);
		static std::map<std::string, std::string> loadDNACodonTable();
		static std::string getAminoAcids( const std::string &str, const std::map<std::string, std::string> &table);
		static int max3( const int a, const int b, const int c);
		static void getFilesInDir( const std::string& path, std::vector<std::string>& filenames);
};
#endif
