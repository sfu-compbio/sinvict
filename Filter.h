#ifndef FILTER_H
#define FILTER_H

#include <string>
#include <vector>
#include <unordered_map>
#include "Location.h"

class Filter
{
	public:
		static std::vector<Location> applyBaseThresholds( std::unordered_map<std::string, Location> candidateLocations);
		static double illuminaPoissonFilter( const int count, const double lambda);
};
#endif
