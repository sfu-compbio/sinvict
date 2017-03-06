#ifndef STATISTICS_H
#define STATISTICS_H

#include <vector>

class Statistics
{
	public:
		static double hypergeometricProbability( int a, int b, int c, int d);
		static double leftTailedFisher( int a, int b, int c, int d);
		static double poissonPDF(const double k, const double lambda);
		static double poissonCDF(const double k, const double lambda);
		static double mean( const std::vector<double> values);
		static double variance( const std::vector<double> values, const double mean);
		static double standardDeviation( const double variance);
		static double coefficientOfVariation( const double mean, const double standardDeviation);
};
#endif
