#ifndef BATCH_H
#define BATCH_H

#include <fstream>
#include <string>
#include <vector>
#include "Location.h"

class Batch
{
	private:
		std::string path;
		std::vector<Location> locations;

	public:
		Batch( const std::string path);
		std::vector<Location> getLocations();
		void addLocation( Location newLocation);
		void printBatch();
		void printBatch( std::ofstream& out, int usePoissonGermline);
};
#endif
