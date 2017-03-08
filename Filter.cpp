#include <iostream>
#include <fstream>
#include <cmath>
#include "Filter.h"
#include "Statistics.h"

#define MIN_DEPTH 20
#define MIN_NONREF_COUNT 10

std::vector<Location> Filter::applyBaseThresholds( std::unordered_map<std::string, Location> candidateLocations)
{
	std::vector<Location> newCandidateLocations;
	std::unordered_map<std::string, Location>::iterator iter;
	for( iter = candidateLocations.begin(); iter != candidateLocations.end(); ++iter)
	{
		Location currentLocation = iter->second;
		
		std::vector<Sample> samples = currentLocation.getSamples();
		bool passFilter = true;
		
		for( int i = 0; i < samples.size(); i++)
		{
			ReadcountEntry currentReadcountEntry = samples[i].getReadcountEntry();
			Allele mostFreqVariantAllele = currentReadcountEntry.getMostFreqVariantAllele();
			if( currentReadcountEntry.getReadDepth() < MIN_DEPTH || mostFreqVariantAllele.getCount() < MIN_NONREF_COUNT)
			{
				passFilter = false;
				break;
			}
		}

		if( passFilter)
		{
			newCandidateLocations.push_back( currentLocation);
		}
		
	}

	return newCandidateLocations;
}

double Filter::illuminaPoissonFilter(const int count, const double lambda)
{
	double pval;
	double qscore;
	if(count >= 1)
	{
		pval = (double) 1 - Statistics::poissonCDF((double)(count - 1), lambda);
		if(pval < 0)
		{
			pval = 0;
		}
	}
	else
	{
		pval = 1;
	}

	return pval;
}
