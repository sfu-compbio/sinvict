#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <cmath>
#include "Caller.h"
#include "Statistics.h"
#include "Filter.h"
#include "Common.h"

Caller::Caller( const double poissonLambda, const int minDepth, const double leftStrandBias, const double rightStrandBias, const double readEndFraction, const int qCutoff, const char* tumorDirectoryPath, const char* benignDirectoryPath, const char* outputDirectoryPath, const int usePoissonGermline, const int disableLvl5Filter)
{
	this->poissonLambda = poissonLambda;
	this->minDepth = minDepth;
	this->strandBiasLeft = leftStrandBias;
	this->strandBiasRight = rightStrandBias;
	this->minQScore = qCutoff;
	this->readEndFraction = readEndFraction;
	this->usePoissonGermline = usePoissonGermline;
	this->disableLvl5Filter = disableLvl5Filter;

	( this->tumorDirectoryPath).assign( tumorDirectoryPath);
	( this->benignDirectoryPath).assign( benignDirectoryPath);
	( this->outputDirectoryPath).assign( outputDirectoryPath);

	Common::getFilesInDir( this->tumorDirectoryPath, tumorPaths);
	std::sort( tumorPaths.begin(), tumorPaths.end());

	Common::getFilesInDir( this->benignDirectoryPath, benignPaths);
	std::sort( benignPaths.begin(), benignPaths.end());

	for( int i = 0; i < tumorPaths.size(); i++)
	{
		std::string nextTumorSamplePath = this->tumorDirectoryPath + "/" + tumorPaths[i];
		std::string nextBenignSamplePath = this->benignDirectoryPath + "/" + benignPaths[i];

		loadEntries( nextTumorSamplePath);
		loadEntries( nextBenignSamplePath);
	}
}

Caller::Caller( const double poissonLambda, const int minDepth, const double leftStrandBias, const double rightStrandBias, const double readEndFraction, const int qCutoff, const char* tumorDirectoryPath, const char* outputDirectoryPath, const int usePoissonGermline, const int disableLvl5Filter)
{
	this->poissonLambda = poissonLambda;
	this->minDepth = minDepth;
	this->strandBiasLeft = leftStrandBias;
	this->strandBiasRight = rightStrandBias;
	this->minQScore = qCutoff;
	this->readEndFraction = readEndFraction;
	this->usePoissonGermline = usePoissonGermline;
	this->disableLvl5Filter = disableLvl5Filter;

	( this->tumorDirectoryPath).assign( tumorDirectoryPath);
	( this->outputDirectoryPath).assign( outputDirectoryPath);

	// Setup output files
	( this->outputPaths).push_back( this->outputDirectoryPath + "/calls_level1.sinvict");
	( this->outputPaths).push_back( this->outputDirectoryPath + "/calls_level2.sinvict");
	( this->outputPaths).push_back( this->outputDirectoryPath + "/calls_level3.sinvict");
	( this->outputPaths).push_back( this->outputDirectoryPath + "/calls_level4.sinvict");
	( this->outputPaths).push_back( this->outputDirectoryPath + "/calls_level5.sinvict");
	( this->outputPaths).push_back( this->outputDirectoryPath + "/calls_level6.sinvict");

	Common::getFilesInDir( this->tumorDirectoryPath, tumorPaths);
	std::sort( tumorPaths.begin(), tumorPaths.end());

	for( int i = 0; i < tumorPaths.size(); i++)
	{
		std::string nextTumorSamplePath = this->tumorDirectoryPath + "/" + tumorPaths[i];

		loadEntries( nextTumorSamplePath);
	}

	calculateStatistics();
}

int Caller::loadEntries( const std::string path)
{
	std::string nextLine;
	std::string key;
	std::string chr;
	std::string refBase;
	int readDepth;
	int pos;

	// Open the sample file
	std::ifstream inputFile( path.c_str());
	if( !inputFile.is_open())
	{
		perror( "Error opening input readcount file");
		exit( 1);
	}

	// For each line in the sample file (which will correspond to a genomic location)
	while( std::getline( inputFile, nextLine))
	{
		// Split the line into tokens separated by whitespace (for columns, since this is a tab delimited file)
		std::istringstream strStream( nextLine);
		std::istream_iterator<std::string> begin( strStream), end;
		std::vector<std::string> stringTokens( begin, end);

		// Get all main fields
		chr = stringTokens[0];
		pos = atoi( stringTokens[1].c_str());
		refBase = stringTokens[2];
		refBase[0] = toupper( refBase[0]);
		readDepth = atoi( stringTokens[3].c_str());

		// Generate the key, (chr:pos)
		key = stringTokens[0] + ":" + stringTokens[1];
		if( key == "")
		{
			std::cout << "Empty key" << std::endl;
		}

		// Create the base ReadcountEntry object
		ReadcountEntry nextReadcountEntry( refBase, readDepth);

		// Get all subfields for each allele, the 5th column (stringTokens[4]) is garbage due to a bug with the bam-readcount program, ignore it
		for( int i = 5; i < stringTokens.size(); i++)
		{
			std::vector<std::string> nextSubTokens = Common::split( stringTokens[i], ":", true);

			// Create the Allele objects and add them to the current ReadcountEntry object
			std::string base = nextSubTokens[0];
			int count = atoi( nextSubTokens[1].c_str());
			double avgMappingQuality = atof( nextSubTokens[2].c_str());
			double avgBaseQuality = atof( nextSubTokens[3].c_str());
			double avgSEMappingQuality = atof( nextSubTokens[4].c_str());
			int numPlusStrand = atoi( nextSubTokens[5].c_str());
			int numMinusStrand = atoi( nextSubTokens[6].c_str());
			double avgPosAsFraction = atof( nextSubTokens[7].c_str());
			double avgNumMismatchesAsFraction = atof( nextSubTokens[8].c_str());
			double avgSumMismatchQualities = atof( nextSubTokens[9].c_str());
			int numQ2ContainingReads = atoi( nextSubTokens[10].c_str());
			double avgDistanceToQ2StartInQ2Reads = atof( nextSubTokens[11].c_str());
			double avgClippedLength = atof( nextSubTokens[12].c_str());
			double avgDistanceToEffective3pEnd = atof( nextSubTokens[13].c_str());

			bool variant = false;
			if( base != refBase)
			{
				variant = true;
			}

			double percentage;
			if( readDepth != 0)
			{
				percentage = ( double) count / ( double) readDepth * 100;
			}
			else
			{
				percentage = 0;
			}

			Allele nextAllele( base, count, avgMappingQuality, avgBaseQuality, avgSEMappingQuality, numPlusStrand, numMinusStrand,
							   avgPosAsFraction, avgNumMismatchesAsFraction, avgSumMismatchQualities, numQ2ContainingReads,
							   avgDistanceToQ2StartInQ2Reads, avgClippedLength, avgDistanceToEffective3pEnd, percentage, variant);

			nextReadcountEntry.addAllele( nextAllele);
		}

		// Now, the ReadcountEntry object is filled, so we can create the Sample object
		nextReadcountEntry.setMostFreqVariantAllele();
		Sample nextSample( path, nextReadcountEntry);

		// Finally, add the Sample object to the Location object,
		// Check if the Location object with the current key exists in the hash table
		std::unordered_map<std::string, Location>::iterator iter = locationTable.find( key);
		if( iter == locationTable.end())
		{
			// If it does not exist, create the Location object
			Location newLocation( chr, pos);

			// Add the new Sample to the Location object
			newLocation.addSample( nextSample);

			// Insert the new key-Location pair to the hash table
			std::pair<std::string, Location> newKeyPair( key, newLocation);
			locationTable.insert( newKeyPair);
		}
		else
		{
			bool sampleExists = false;
			std::vector<Sample> samples = ( iter->second).getSamples();
			for( int j = 0; j < samples.size(); j++)
			{
				if( samples[j].getSampleName() == nextSample.getSampleName())
				{
					sampleExists = true;
				}
			}

			if( !sampleExists)
			{
				( iter->second).addSample( nextSample);
			}
		}
	}

	// Check if the file was read correctly
	if( inputFile.bad())
	{
		perror( "Error reading input readcount file");
	}

	// Close the input sample file
	inputFile.close();
}

void Caller::calculateStatistics()
{
	std::unordered_map<std::string, Location>::iterator iter;
	for( iter = locationTable.begin(); iter != locationTable.end(); ++iter)
	{
		std::vector<double> variantPercentages;
		std::vector<Sample> sampleList = ( iter->second).getSamples();
		for( int i = 0; i < sampleList.size(); i++)
		{
			ReadcountEntry re = sampleList[i].getReadcountEntry();
			Allele mostFreqVariant = re.getMostFreqVariantAllele();
			variantPercentages.push_back( mostFreqVariant.getPercentage());
		}

		// Calculate mean
		double mean = Statistics::mean( variantPercentages);

		// Calculate variance
		double variance = Statistics::variance( variantPercentages, mean);

		// Calculate std
		double std = Statistics::standardDeviation( variance);

		// Calculate snr
		double cov = Statistics::coefficientOfVariation( mean, std);

		// Set statistics for the current Location
		( iter->second).setMeanVAP( mean);
		( iter->second).setVarianceVAP( variance);
		( iter->second).setStdVAP( std);
		( iter->second).setCOV( cov);
	}
}

std::vector<Location> Caller::callPoissonDist( double poissonLambda, int minQScore)
{
	std::vector<Location> newCandidateLocations;
	std::unordered_map<std::string, Location>::iterator iter;
	std::string altBase;
	for( iter = locationTable.begin(); iter != locationTable.end(); ++iter)
	{
		Location newLocation = iter->second;

		// Clear the Sample list of the copy of the location
		newLocation.clearSamples();
		bool keepLocation = false;

		std::vector<Sample> sampleList = ( iter->second).getSamples();
		for( int i = 0; i < sampleList.size(); i++)
		{
			ReadcountEntry readcountEntry = sampleList[i].getReadcountEntry();
			Allele mostFreqVariantAllele = readcountEntry.getMostFreqVariantAllele();

			int mostFreqNonRefCount = mostFreqVariantAllele.getCount();
			double lambda = readcountEntry.getReadDepth() * poissonLambda;

			// call illuminaPoissonFilter
			double pValue = Filter::illuminaPoissonFilter( mostFreqNonRefCount, lambda);
			double qScore = -10 * std::log10( pValue);

			// if at least one Sample passes through the filter, keep the location
			if( qScore > minQScore)
			{
				//mostFreqVariantAllele.setPValue( pValue);
				//mostFreqVariantAllele.setQScore( qScore);

				// Add only the called Samples to the emptied list
				newLocation.addSample( sampleList[i]);
				keepLocation = true;
			}
		}

		std::vector<Sample> newSamples = newLocation.getSamples();
		double highestVAP = -1;
		for( int i = 0; i < newSamples.size(); i++)
		{
			ReadcountEntry readcountEntry = newSamples[i].getReadcountEntry();
			Allele variantAllele = readcountEntry.getMostFreqVariantAllele();

			if( variantAllele.getPercentage() > highestVAP)
			{
				highestVAP = variantAllele.getPercentage();
				altBase = variantAllele.getBase();
			}
		}

		( iter->second).setMutatedBase( altBase);
		if( keepLocation)
		{
			newCandidateLocations.push_back( newLocation);
		}
	}
	return newCandidateLocations;
}

std::vector<Location> Caller::callAverageFilter( std::vector<Location> unfilteredCalls)
{
	std::vector<Location> newCandidateLocations;
	for( int i = 0; i < unfilteredCalls.size(); i++)
	{
		bool keepLocation = false;

		// Make a copy of the Location for the new list of candidate locations
		Location newLocation = unfilteredCalls[i];

		// Clear the Sample list of the copy of the location
		newLocation.clearSamples();

		std::vector<Sample> sampleList = unfilteredCalls[i].getSamples();
		int numSamples = sampleList.size();
		for( int j = 0; j < numSamples; j++)
		{
			// Get the mean of the current location
			double mean = unfilteredCalls[i].getMeanVAP();

			double currentSamplePercentage = sampleList[j].getReadcountEntry().getMostFreqVariantAllele().getPercentage();

			// If the current Sample's non reference percentage is  5 times or more greater than the new average, pass
			if( currentSamplePercentage >= 3 * mean)
			{
				// If at least one Sample passes through the filter, keep the location
				// Add only the called Samples to the emptied list
				newLocation.addSample( sampleList[j]);
				keepLocation = true;
			}
		}

		if( keepLocation)
		{
			newCandidateLocations.push_back( newLocation);
		}
	}
	return newCandidateLocations;
}

std::vector<Location> Caller::callDepthFilter( std::vector<Location> unfilteredCalls, int minDepth)
{
	std::vector<Location> newCandidateLocations;
	for( int i = 0; i < unfilteredCalls.size(); i++)
	{
		bool keepLocation = false;

		// Make a copy of the Location for the new list of candidate locations
		Location newLocation = unfilteredCalls[i];

		// Clear the Sample list of the copy of the location
		newLocation.clearSamples();

		std::vector<Sample> sampleList = unfilteredCalls[i].getSamples();
		int numSamples = sampleList.size();
		for( int j = 0; j < numSamples; j++)
		{
			// Get the read depth
			int readDepth = sampleList[j].getReadcountEntry().getReadDepth();

			if( readDepth >=  minDepth)
			{
				// If at least one Sample passes through the filter, keep the location
				// Add only the called Samples to the emptied list
				newLocation.addSample( sampleList[j]);
				keepLocation = true;
			}
		}

		if( keepLocation)
		{
			newCandidateLocations.push_back( newLocation);
		}
	}
	return newCandidateLocations;
}

std::vector<Location> Caller::callStrandBiasFilter( std::vector<Location> unfilteredCalls, double strandBiasLeft, double strandBiasRight)
{
	std::vector<Location> newCandidateLocations;
	for( int i = 0; i < unfilteredCalls.size(); i++)
	{
		bool keepLocation = false;

		// Make a copy of the Location for the new list of candidate locations
		Location newLocation = unfilteredCalls[i];

		// Clear the Sample list of the copy of the location
		newLocation.clearSamples();

		std::vector<Sample> sampleList = unfilteredCalls[i].getSamples();
		int numSamples = sampleList.size();
		for( int j = 0; j < numSamples; j++)
		{
			// Calculate the strand-bias
			int numReadsForward = sampleList[j].getReadcountEntry().getMostFreqVariantAllele().getNumPlusStrand();
			int numReadsReverse = sampleList[j].getReadcountEntry().getMostFreqVariantAllele().getNumMinusStrand();
			double strandBias = ( double) numReadsForward / ( double) ( numReadsForward + numReadsReverse);

			if( strandBias >= strandBiasLeft && strandBias <= strandBiasRight)
			{
				// If at least one Sample passes through the filter, keep the location
				// Add only the called Samples to the emptied list
				newLocation.addSample( sampleList[j]);
				keepLocation = true;
			}
		}

		if( keepLocation)
		{
			newCandidateLocations.push_back( newLocation);
		}
	}
	return newCandidateLocations;
}

std::vector<Location> Caller::callAmpliconEndFilter( std::vector<Location> unfilteredCalls, double readEndFraction)
{
	std::vector<Location> newCandidateLocations;
	for( int i = 0; i < unfilteredCalls.size(); i++)
	{
		bool keepLocation = false;

		// Make a copy of the Location for the new list of candidate locations
		Location newLocation = unfilteredCalls[i];

		// Clear the Sample list of the copy of the location
		newLocation.clearSamples();

		std::vector<Sample> sampleList = unfilteredCalls[i].getSamples();
		int numSamples = sampleList.size();
		for( int j = 0; j < numSamples; j++)
		{
			// Get average position on the reads as a fraction
			double avgPosAsFraction = sampleList[j].getReadcountEntry().getMostFreqVariantAllele().getAvgPosAsFraction();

			if( avgPosAsFraction >= readEndFraction)
			{
				// If at least one Sample passes through the filter, keep the location
				// Add only the called Samples to the emptied list
				newLocation.addSample( sampleList[j]);
				keepLocation = true;
			}
		}

		if( keepLocation)
		{
			newCandidateLocations.push_back( newLocation);
		}
	}
	return newCandidateLocations;
}

std::vector<Location> Caller::callHomopolymerFilter( std::vector<Location> unfilteredCalls)
{
	std::vector<Location> newCandidateLocations;
	for( int i = 0; i < unfilteredCalls.size(); i++)
	{
		// Get the chr and pos of the current location
		std::string chr = unfilteredCalls[i].getChr();
		int pos = unfilteredCalls[i].getPosition();

		// Create keys for the adjacent positions
		std::vector<std::string> keys;
		keys.push_back( chr + ":" + std::to_string( pos - 1));
		keys.push_back( chr + ":" + std::to_string( pos - 2));
		keys.push_back( chr + ":" + std::to_string( pos - 3));
		keys.push_back( chr + ":" + std::to_string( pos + 1));
		keys.push_back( chr + ":" + std::to_string( pos + 2));
		keys.push_back( chr + ":" + std::to_string( pos + 3));

		// Get the Location objects from the lookup table with the keys
		std::vector<Location> adjacentLocations;
		bool allLocationsExist = true;
		for( int j = 0; j < keys.size(); j++)
		{
			std::unordered_map<std::string, Location>::iterator iter = locationTable.find( keys[j]);
			if( iter == locationTable.end())
			{
				allLocationsExist = false;
				break;
			}
			else
			{
				adjacentLocations.push_back( iter->second);
			}
		}

		bool leftNotHomopolymer = true;
		bool rightNotHomopolymer = true;
		if( allLocationsExist)
		{
			// Check if the right and left three bases are identical to each other
			std::vector<std::string> bases;
			for( int j = 0; j < adjacentLocations.size(); j++)
			{
				std::vector<Sample> samples = adjacentLocations[j].getSamples();
				bases.push_back( samples[0].getReadcountEntry().getRefBase());
			}

			if( bases[0] == bases[1] && bases[1] == bases[2])
			{
				leftNotHomopolymer = false;
			}

			if( bases[3] == bases[4] && bases[4] == bases[5])
			{
				rightNotHomopolymer = false;
			}
		}

		// Keep the location if there are at least 3 adjacent positions on each side of the current location,
		// AND these 3 bases are not identical at any side
		if( allLocationsExist && leftNotHomopolymer &&  rightNotHomopolymer)
		{
			newCandidateLocations.push_back( unfilteredCalls[i]);
		}
	}
	return newCandidateLocations;
}

int Caller::callLocationsMixture()
{
	std::ofstream outputFile1;
	std::ofstream outputFile2;
	std::ofstream outputFile3;
	std::ofstream outputFile4;
	std::ofstream outputFile5;
	std::ofstream outputFile6;

	outputFile1.open( outputPaths[0].c_str());
	outputFile2.open( outputPaths[1].c_str());
	outputFile3.open( outputPaths[2].c_str());
	outputFile4.open( outputPaths[3].c_str());
	outputFile5.open( outputPaths[4].c_str());
	outputFile6.open( outputPaths[5].c_str());
	
	// Obtain the initial set of calls using a Poisson Model
	std::vector<Location> firstLevelPass = callPoissonDist( poissonLambda, minQScore);

	// Apply basic read depth filter
	std::vector<Location> secondLevelPass = callDepthFilter( firstLevelPass, minDepth);

	// Print out the initial calls made by the Poisson Model
//	printUCSC( firstLevelPass, outputFile1);
//	printCITUP( firstLevelPass, outputFile1);
	for( int i = 0; i < firstLevelPass.size(); i++)
	{
		firstLevelPass[i].printLocation( outputFile1, usePoissonGermline);
	}

	// Clear out the first level pass vector
	firstLevelPass.clear();

	// Apply strand bias filter
	std::vector<Location> thirdLevelPass = callStrandBiasFilter( secondLevelPass, strandBiasLeft, strandBiasRight);

	// Print out the next level of calls
//	printUCSC( secondLevelPass, outputFile2);
//	printCITUP( secondLevelPass, outputFile2);
	for( int i = 0; i < secondLevelPass.size(); i++)
	{
		secondLevelPass[i].printLocation( outputFile2, usePoissonGermline);
	}

	// Clear out the second level pass vector
	secondLevelPass.clear();

	// Apply read end filter
	std::vector<Location> fourthLevelPass = callAmpliconEndFilter( thirdLevelPass, readEndFraction);

	// Print out next level of calls
//	printUCSC( thirdLevelPass, outputFile3);
//	printCITUP( thirdLevelPass, outputFile3);
	for( int i = 0; i < thirdLevelPass.size(); i++)
	{
		thirdLevelPass[i].printLocation( outputFile3, usePoissonGermline);
	}

	// Clear out the third level pass vector
	thirdLevelPass.clear();

	// Apply average filter
	std::vector<Location> fifthLevelPass;
	if(disableLvl5Filter == 0)
	{
		fifthLevelPass = callAverageFilter( fourthLevelPass);
	}
	else
	{
		fifthLevelPass = fourthLevelPass;
	}

	// Print out next level of calls
//	printUCSC( fourthLevelPass, outputFile4);
//	printCITUP( fourthLevelPass, outputFile4);
	for( int i = 0; i < fourthLevelPass.size(); i++)
	{
		fourthLevelPass[i].printLocation( outputFile4, usePoissonGermline);
	}

	// Clear out the fourth level pass vector
	fourthLevelPass.clear();

	// Apply homopolymer region filter
	std::vector<Location> sixthLevelPass = callHomopolymerFilter( fifthLevelPass);

	// Print out next level of calls
//	printUCSC( fifthLevelPass, outputFile5);
//	printCITUP( fifthLevelPass, outputFile5);
	for( int i = 0; i < fifthLevelPass.size(); i++)
	{
		fifthLevelPass[i].printLocation( outputFile5, usePoissonGermline);
	}

	// Clear out the fifth level pass vector
	fifthLevelPass.clear();

	// Print out the next level of calls
	for( int i = 0; i < sixthLevelPass.size(); i++)
	{
		sixthLevelPass[i].printLocation( outputFile6, usePoissonGermline);
		//sixthLevelPass[i].printLocationVCF( outputFile6);
	}

//	printUCSC( sixthLevelPass, outputFile6);
//	printCITUP( sixthLevelPass, outputFile6);

	// Clear out the sixth level pass vector
	sixthLevelPass.clear();

	outputFile1.close();
	outputFile2.close();
	outputFile3.close();
	outputFile4.close();
	outputFile5.close();
	outputFile6.close();

	return 0;
}

void Caller::printUCSC( std::vector<Location> locations, std::ofstream& outputFile)
{
	if( locations.size() > 0)
	{
		std::string chr = locations[0].getChr();
		int pos = locations[0].getPosition();
		std::string key = chr + ":" + std::to_string( pos);

		std::unordered_map<std::string, Location>::iterator iter = locationTable.find( key);
		Location nextLocation = iter->second;
	//	nextLocation.printUCSCHeader( outputFile);

		for( int i = 0; i < locations.size(); i++)
		{
			// Get the chr and pos of the current location
			chr = locations[i].getChr();
			pos = locations[i].getPosition();
			key = chr + ":" + std::to_string( pos);

			iter = locationTable.find( key);
			nextLocation = iter->second;

//			nextLocation.printLocationUCSC( outputFile);
			nextLocation.printLocationANNOVAR( outputFile);
		}
	}
}

void Caller::printCITUP( std::vector<Location> locations, std::ofstream& outputFile)
{

	outputFile << "Num_mutations: " << locations.size() << "\n";
	outputFile << "Num_samples: 1" << "\n";
	outputFile << "Error_rate: 0.001" << "\n";

	for( int i = 0; i < locations.size(); i++)
	{
		locations[i].printLocationCITUP( outputFile);
	}
}

void Caller::printCaller()
{
	std::cout << "Tumor Readcount Directory: " << tumorDirectoryPath << "\t" << "Output Directory: " << outputDirectoryPath << "\n\n";
	std::cout << "Tumor Readcount Filenames: " << "\n";
	std::cout << "=========================== " << "\n"; 
	for( int i = 0; i < tumorPaths.size(); i++)
	{
		std::cout << tumorPaths[i] << "\n";
	}
}

void Caller::printLocations()
{
	std::unordered_map<std::string, Location>::iterator iter;
	for( iter = locationTable.begin(); iter != locationTable.end(); ++iter)
	{
		Location currentLocation = iter->second;
		currentLocation.printLocation();
	}
}
