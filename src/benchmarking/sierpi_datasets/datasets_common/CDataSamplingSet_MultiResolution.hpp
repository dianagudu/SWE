/*
 * Copyright (C) 2011 Technische Universitaet Muenchen
 * This file is part of the Sierpi project. For conditions of distribution and
 * use, please see the copyright notice in the file 'copyright.txt' at the root
 * directory of this package and the copyright notice at http://www5.in.tum.de/sierpi
 *
 *  Created on: July 19, 2012
 *      Author: Martin Schreiber <martin.schreiber@in.tum.de>
 */


#ifndef CDATASAMPLINGSET_MULTIRESOLUTION_HPP
#define CDATASAMPLINGSET_MULTIRESOLUTION_HPP

#include <cassert>
#include <netcdf.hh>
#include <stdlib.h>
#include "CDataSamplingSet_SingleLevel.hpp"


#if CONFIG_COMPILE_WITHOUT_SIERPI
#	ifndef nullptr
#		define nullptr		0
#	endif
#endif

/**
 * This class represents a multi-resolution data sampling set
 *
 */
template <typename T>
class CDataSamplingSet_MultiResolution
{
public:
	/**
	 * data levels with different resolutions
	 */
	CDataSamplingSet_SingleLevel<T> *dataLevels;

	/**
	 * number of data levels
	 */
	int number_of_data_levels;


	/**
	 * minimum x-coordinate
	 */
	T xmin;

	/**
	 * minimum y-coordinate
	 */
	T ymin;


	/**
	 * maximum x-coordinate
	 */
	T xmax;

	/**
	 * maximum y-coordinate
	 */
	T ymax;

	/**
	 * domain x length
	 */
	T domain_length_x;

	/**
	 * domain y length
	 */
	T domain_length_y;

	/**
	 * verbosity level
	 */
	int verbosity_level;


public:
	void setOutOfDomainValue(T i_value)
	{
		for (int i = 0; i < number_of_data_levels; i++)
		{
			dataLevels[i].out_of_domain_value = i_value;
		}
	}

	size_t getArraySizeX()
	{
		assert(number_of_data_levels > 0);
		return dataLevels[number_of_data_levels-1].getArraySizeX();
	}

	size_t getArraySizeY()
	{
		assert(number_of_data_levels > 0);
		return dataLevels[number_of_data_levels-1].getArraySizeY();
	}

	T getXMin()
	{
		return xmin;
	}

	T getXMax()
	{
		return xmax;
	}

	T getYMin()
	{
		return ymin;
	}

	T getYMax()
	{
		return ymax;
	}

	T getFloat2D(
		T i_x,
		T i_y,
		T i_level
	)
	{
		int l = (int)i_level;
		int nl = l+1;

		l = std::min(l, number_of_data_levels-1);
		nl = std::min(nl, number_of_data_levels-1);

		T val = dataLevels[l].getFloat2D_LinearInterpolation(i_x, i_y);
		T nval = dataLevels[nl].getFloat2D_LinearInterpolation(i_x, i_y);

		T ratio = (i_level-(T)l);
		return ratio*nval + ((T)1.0-ratio)*val;
	}


	T getFloat2D_BorderValue(
		T i_x,
		T i_y,
		T i_level
	)
	{
		int l = (int)i_level;
		int nl = l+1;

		l = std::min(l, number_of_data_levels-1);
		nl = std::min(nl, number_of_data_levels-1);

		T val = dataLevels[l].getFloat2D_BorderValue(i_x, i_y);
		T nval = dataLevels[nl].getFloat2D_BorderValue(i_x, i_y);

		T ratio = (i_level-(T)l);
		return ratio*nval + ((T)1.0-ratio)*val;
	}


	T getFloat2D_NearestNeighbor(
		T i_x,
		T i_y,
		T i_level
	)
	{
		int l = std::min(i_level, number_of_data_levels-1);

		return dataLevels[l].getFloat2D_NearestNeighbor(i_x, i_y);
	}


	/**
	 * constructor
	 */
	CDataSamplingSet_MultiResolution(
			int i_verbosity_level
	)	:
		verbosity_level(i_verbosity_level)
	{

		dataLevels = nullptr;
		number_of_data_levels = 0;
	}

	/**
	 * deconstructor
	 */
	~CDataSamplingSet_MultiResolution()
	{
		clear();
	}


	/**
	 * clear / free data
	 */
	void clear()
	{
		if (dataLevels)
		{
	        for (int i = number_of_data_levels-1; i >= 0; i--)
	        	dataLevels[i].clear();

			delete [] dataLevels;
			dataLevels = nullptr;
		}
	}



	/**
	 * setup multiresolution levels and return pointer to allocated array structure
	 *
	 * \return pointer to first (finest) level
	 */
	T *allocateMultiresolutionLevels(
		size_t i_array_size_x,		///< array elements in x-direction for finest resolution
		size_t i_array_size_y		///< array elements in y-direction for finest resolution
	)
	{
		/*
		 * determine number of multiresolution levels
		 */
		number_of_data_levels = 1;
		int max_size = std::max(i_array_size_x, i_array_size_y);

		while (max_size >> number_of_data_levels)
			number_of_data_levels++;

		if (1 << (number_of_data_levels-1) != max_size)
			number_of_data_levels++;

		dataLevels = new CDataSamplingSet_SingleLevel<T>[number_of_data_levels];

		/*
		 * setup first level
		 */
		int l = number_of_data_levels-1;
		dataLevels[l].setupArray(i_array_size_x, i_array_size_y);

		// return pointer to first level
		return dataLevels[l].getArrayPointer();
	}


	/**
	 * return the number of multi-resolution levels
	 */
	int getLevels()
	{
		return number_of_data_levels;
	}


	/**
	 * setup domain parameters for first level
	 */
	void setupDomainParameters(
			T i_xmin,
			T i_ymin,
			T i_xmax,
			T i_ymax
		)
	{
		xmin = i_xmin;
		ymin = i_ymin;

		xmax = i_xmax;
		ymax = i_ymax;

		domain_length_x = xmax - xmin;
		domain_length_y = ymax - ymin;

		if (dataLevels == nullptr)
		{
			std::cerr << "Parameters have to be set up AFTER allocateMultiresolutionLevels()" << std::endl;
			assert(false);
			exit(-1);
		}

		dataLevels[number_of_data_levels-1].setupDomainParameters(xmin, ymin, xmax, ymax);
	}



	/**
	 * setup the multiresolution map
	 */
	void setupMultiresolutionMap()
	{
		for (int i = number_of_data_levels-2; i >= 0; i--)
		{
			dataLevels[i].setupFromFinerLevel(dataLevels[i+1], verbosity_level);
			dataLevels[i].setupDomainParameters(dataLevels[i+1]);
			dataLevels[i].setOutOfDomainValue(0);
		}
	}



	/**
	 * output array values
	 */
	void debugPrintArray(int i_level)
	{
		assert(number_of_data_levels > i_level);

		dataLevels[i_level].debugPrintArray();
	}
};

#endif
