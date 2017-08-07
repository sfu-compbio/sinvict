#include <iostream>
#include "Batch.h"

Batch::Batch( const std::string path)
{
	this->path = path;
}

std::vector<Location> Batch::getLocations()
{
	return locations;
}

void Batch::addLocation( Location newLocation)
{
	locations.push_back( newLocation);
}

void Batch::printBatch()
{
	std::cout << "Batch at directory: " << path << "\n\n";

	std::cout << "List of Locations: " << "\n";
	std::cout << "================== " << "\n";
	for( int i = 0; i < locations.size(); i++)
	{
		locations[i].printLocation();
	}
}

void Batch::printBatch( std::ofstream& out, int usePoissonGermline)
{
	out << "Batch at directory: " << path << "\n\n";

	out << "List of Locations: " << "\n";
	out << "================== " << "\n";
	for( int i = 0; i < locations.size(); i++)
	{
		locations[i].printLocation( out, usePoissonGermline);
	}
}
