#include <cmath>
#include <boost/math/distributions/poisson.hpp>
#include "Statistics.h"

double Statistics::leftTailedFisher( int a, int b, int c, int d)
{
	int n = a + b + c + d;
	double p = 0;
	int min;

	// Calculate the p-value for the initial table
	p = Statistics::hypergeometricProbability( a, b, c, d);

	if( a < d)
	{
		min = a;
	}
	else
	{
		min = d;
	}

	// For every other valid table, calculate the p-value for that table and add it to the original p-value
	for( int i = 0; i < min; i++)
	{
		a = a - 1;
		d = d - 1;
		b = b + 1;
		c = c + 1;

		double pNew = Statistics::hypergeometricProbability( a, b, c, d);
		p = p + pNew;
	}

	return p;
}

double Statistics::hypergeometricProbability( int a, int b, int c, int d)
{
	int n = a + b + c + d;
	double p = 1.0;
	int j = n;

	for( int i = a + 1; i <= a + c; i++)
	{
		p = p * i;
		if( p > 1)
		{
			p = p / ( double) j;
			j--;
		}
	}

	for( int i = c + 1; i <= c + d; i++)
	{
		p = p * i;
		if( p > 1)
		{
			p = p / ( double) j;
			j--;
		}
	}

	for( int i = b + 1; i <= b + d; i++)
	{
		p = p * i;
		if( p > 1)
		{
			p = p / ( double) j;
			j--;
		}
	}

	while( j > a + b + 1)
	{
		p = p / ( double) j;
		j--;
	}

	for( int i = 1; i <= c; i++)
	{
		p = p / ( double) i;
	}

	return p;
}

double Statistics::poissonCdf( int x, double lambda)
{
	boost::math::poisson_distribution<> poissonDist( lambda);
	return boost::math::cdf( poissonDist, x);
}

double Statistics::mean( const std::vector<double> values)
{
	double sum = 0;
	int numberOfElements = values.size();

	for( int i = 0; i < numberOfElements; i++)
	{
		sum = sum + values[i];
	}

	double mean = ( double) sum / ( double) numberOfElements;
	return mean;
}

double Statistics::variance( const std::vector<double> values, const double mean)
{
	double sum = 0;
	int numberOfElements = values.size();

	for( int i = 0; i < numberOfElements; i++)
	{
		double difference = values[i] - mean;
		double squaredDifference = pow( difference, 2.0);
		sum = sum + squaredDifference;
	}

	double variance = ( double) sum / ( double) numberOfElements;
	return variance;
	
}

double Statistics::standardDeviation( const double variance)
{
	double standardDeviation = sqrt( variance);
	return standardDeviation;
}

double Statistics::coefficientOfVariation( const double mean, const double standardDeviation)
{
	double cov = standardDeviation / mean;
	return cov;
}
