/*
 * Copyright (C) 2011 Technische Universitaet Muenchen
 * This file is part of the Sierpi project. For conditions of distribution and
 * use, please see the copyright notice in the file 'copyright.txt' at the root
 * directory of this package and the copyright notice at http://www5.in.tum.de/sierpi
 *
 *  Created on: July 24, 2012
 *      Author: Martin Schreiber <martin.schreiber@in.tum.de>
 */

#ifndef CDATASAMPLINGSET_SINGLELEVEL_HPP__
#define CDATASAMPLINGSET_SINGLELEVEL_HPP__

#include <cmath>


#if CONFIG_COMPILE_WITHOUT_SIERPI
#	ifndef nullptr
#		define nullptr		0
#	endif
#endif


/**
 * this class represents a single level in a multiresolution texture
 */
template <typename T>
class CDataSamplingSet_SingleLevel
{
private:
	T *data;

	size_t array_size_x;
	size_t array_size_y;

	T out_of_domain_value;

	/**
	 * domain parameters
	 */
	T domain_start_x;
	T domain_start_y;

	T domain_length_x;
	T domain_length_y;

	T domain_end_x;
	T domain_end_y;


	/**
	 * translation components to project from world-space to array-space
	 */
	T to_texel_coord_translate_x;
	T to_texel_coord_translate_y;
	T to_texel_coord_scale_x;
	T to_texel_coord_scale_y;



public:
	size_t getArraySizeX()
	{
		return array_size_x;
	}

	size_t getArraySizeY()
	{
		return array_size_y;
	}


	CDataSamplingSet_SingleLevel()	:
		data(nullptr),
		array_size_x(0),
		array_size_y(0)
	{
	}


	/**
	 * setup domain parameters (size)
	 */
	void setupDomainParameters(
			T i_domain_start_x,
			T i_domain_start_y,

			T i_domain_end_x,
			T i_domain_end_y
		)
	{
		domain_start_x = i_domain_start_x;
		domain_start_y = i_domain_start_y;

		domain_end_x = i_domain_end_x;
		domain_end_y = i_domain_end_y;

		domain_length_x = domain_end_x - domain_start_x;
		domain_length_y = domain_end_y - domain_start_y;

		p_setupSamplingHelperVariables();
	}



	/**
	 * setup the domain parameters using a parent's (finer) level
	 */
	void setupDomainParameters(
		CDataSamplingSet_SingleLevel &i_parent
	)	{
		domain_start_x = i_parent.domain_start_x;
		domain_start_y = i_parent.domain_start_y;

		domain_end_x = i_parent.domain_end_x;
		domain_end_y = i_parent.domain_end_y;

		domain_length_x = domain_end_x - domain_start_x;
		domain_length_y = domain_end_y - domain_start_y;

		p_setupSamplingHelperVariables();
	}


	/**
	 * setup helper variables to compute discrete sampling positions faster
	 *
	 *
	 * 3 4 5 6 7 8 9  <- world-space
	 * |   |   |   |
	 * |   |   |   |
	 * +---+---+---+
	 * | + | + | + |
	 * +---+---+---+
	 *   0   1   2    <- array-space
	 *
	 * dx = cell-size = domain_length_x / array_size_x
	 *
	 */
	void p_setupSamplingHelperVariables()
	{
		T cell_size_x = domain_length_x / (T)array_size_x;
		T cell_size_y = domain_length_y / (T)array_size_y;

		to_texel_coord_translate_x = (domain_start_x + cell_size_x*(T)0.5);
		to_texel_coord_translate_y = (domain_start_y + cell_size_y*(T)0.5);

		assert(array_size_x != 0);
		if (array_size_x > 1)
			to_texel_coord_scale_x = ((T)array_size_x-(T)1.0)/(domain_length_x - cell_size_x);
		else
			to_texel_coord_scale_x = (T)1.0;

		assert(array_size_y != 0);
		if (array_size_y > 1)
			to_texel_coord_scale_y = ((T)array_size_y-(T)1.0)/(domain_length_y - cell_size_y);
		else
			to_texel_coord_scale_y = (T)1.0;

		/*
		 * finally, the array index can be computed by rounding the following values:
		 *
		 * ( (x - translate_x) * scale_x )
		 * ( (y - translate_y) * scale_y )
		 */
	}


	/**
	 * setup array and allocate data storage
	 */
	void setupArray(
			size_t i_array_size_x,	///< array x-size
			size_t i_array_size_y	///< array y-size
	)	{
		array_size_x = i_array_size_x;
		array_size_y = i_array_size_y;

		data = new T[array_size_x*array_size_y];
	}


	/**
	 * return the pointer to the allocated array
	 */
	inline T* getArrayPointer()
	{
		assert(data != nullptr);
		return data;
	}


	/**
	 * free data
	 */
	void clear()
	{
		if (data)
		{
			delete [] data;
			data = nullptr;
		}
	}



	/**
	 * return a texel given at array index
	 */
	T getTexelFloat2D_BorderValue(
			int i_x,		///< x-array index
			int i_y			///< y-array index
	)
	{
		if (i_x < 0 || i_x >= array_size_x)
			return out_of_domain_value;

		if (i_y < 0 || i_y >= array_size_y)
			return out_of_domain_value;

		return data[i_x + i_y*array_size_x];
	}



	/**
	 * return a texel given at array index
	 */
	T getTexelFloat2D_Clamp(
			int i_x,
			int i_y
	)
	{
		if (i_x < 0)
			i_x = 0;

		if (i_x >= (int)array_size_x)
			i_x = array_size_x-1;

		if (i_y < 0)
			i_y = 0;

		if (i_y >= (int)array_size_y)
			i_y = array_size_y-1;

		return data[i_x + i_y*array_size_x];
	}



	/**
	 * return a sampling determined by nearest neighbor
	 */
	T getFloat2D_NearestNeighbor(
			T i_x,
			T i_y
	)	{
		T px = (i_x - to_texel_coord_translate_x)*to_texel_coord_scale_x;
		T py = (i_y - to_texel_coord_translate_y)*to_texel_coord_scale_y;

		// use int to also get negative sign
		int x = (int)(px + (T)0.5);
		int y = (int)(py + (T)0.5);

		return getTexelFloat2D_Clamp(x, y);
	}





	/**
	 * return an interpolated float value
	 */
	T getFloat2D_LinearInterpolation(
			T i_x,
			T i_y
	)	{
		T px = (i_x - to_texel_coord_translate_x)*to_texel_coord_scale_x;
		T py = (i_y - to_texel_coord_translate_y)*to_texel_coord_scale_y;

		/*
		 * now px and py are in texel space with integers aligned at texel cell centers
		 */

		/*
		 * get left & top most cell
		 */
		T x = std::floor(px);
		T y = std::floor(py);

		/*
		 * compute ratios
		 */
		T left_ratio = px - x;
		T right_ratio = (T)1.0 - left_ratio;

		T top_ratio = py - y;
		T bottom_ratio = (T)1.0 - top_ratio;

		int pixel_x = x;
		int pixel_y = y;

		/*
		 * interpolate using areas computed out of ratios
		 */
		return	(right_ratio*bottom_ratio)*getTexelFloat2D_Clamp(pixel_x, pixel_y)	+
				(right_ratio*top_ratio)*getTexelFloat2D_Clamp(pixel_x, pixel_y+1)	+
				(left_ratio*bottom_ratio)*getTexelFloat2D_Clamp(pixel_x+1, pixel_y)	+
				(left_ratio*top_ratio)*getTexelFloat2D_Clamp(pixel_x+1, pixel_y+1);
	}


	/**
	 * return an interpolated float value
	 */
	T getFloat2D_BorderValue(
			T i_x,
			T i_y
	)	{
		T px = (i_x - to_texel_coord_translate_x)*to_texel_coord_scale_x;
		T py = (i_y - to_texel_coord_translate_y)*to_texel_coord_scale_y;

		/*
		 * now px and py are in texel space with integers aligned at texel cell centers
		 */

		/*
		 * get left & top most cell
		 */
		T x = std::floor(px);
		T y = std::floor(py);

		/*
		 * compute ratios
		 */
		T left_ratio = px - x;
		T right_ratio = (T)1.0 - left_ratio;

		T top_ratio = py - y;
		T bottom_ratio = (T)1.0 - top_ratio;

		int pixel_x = x;
		int pixel_y = y;

		/*
		 * interpolate using areas computed out of ratios
		 */
		return	(right_ratio*bottom_ratio)*getTexelFloat2D_BorderValue(pixel_x, pixel_y)	+
				(right_ratio*top_ratio)*getTexelFloat2D_BorderValue(pixel_x, pixel_y+1)	+
				(left_ratio*bottom_ratio)*getTexelFloat2D_BorderValue(pixel_x+1, pixel_y)	+
				(left_ratio*top_ratio)*getTexelFloat2D_BorderValue(pixel_x+1, pixel_y+1);
	}



	void setOutOfDomainValue(
			T i_out_of_domain_value
	)
	{
		out_of_domain_value = i_out_of_domain_value;
	}


	/**
	 * print the currently stored array values for debugging purpose
	 */
	void debugPrintArray()
	{
		for (int y = 0; y < array_size_y; y++)
		{
			for (int x = 0; x < array_size_x; x++)
			{
				std::cout << data[x + array_size_x*y] << "\t";
			}
			std::cout << std::endl;
		}
	}


	/**
	 * setup current level from finer level
	 */
	void setupFromFinerLevel(
			CDataSamplingSet_SingleLevel &i_cSimpleNetCDF_DataLevel_parent,
			int i_verbosity_level
	)	{
		setOutOfDomainValue(i_cSimpleNetCDF_DataLevel_parent.out_of_domain_value);

		// use half array size...
		array_size_x = i_cSimpleNetCDF_DataLevel_parent.array_size_x / 2;
		array_size_y = i_cSimpleNetCDF_DataLevel_parent.array_size_y / 2;

		size_t non_padded_array_size_x = array_size_x;
		size_t non_padded_array_size_y = array_size_y;

		// ... with padding of one when above level had odd size
		if (i_cSimpleNetCDF_DataLevel_parent.array_size_x & 1)
			array_size_x++;
		if (i_cSimpleNetCDF_DataLevel_parent.array_size_y & 1)
			array_size_y++;


//		std::cout << "    PARENT: " << i_cSimpleNetCDF_DataLevel_parent.array_size_x << " " << i_cSimpleNetCDF_DataLevel_parent.array_size_y << std::endl;
//		std::cout << "     CHILD: " << array_size_x << " " << array_size_y << std::endl;

		assert(data == nullptr);

		try
		{
			data = new T[array_size_x*array_size_y];
		}
		catch (...)
		{
			std::cerr << "Error during allocation of Level for multi-resolution map" << std::endl;
			exit(-1);
		}

		domain_start_x = i_cSimpleNetCDF_DataLevel_parent.domain_start_x;
		domain_start_y = i_cSimpleNetCDF_DataLevel_parent.domain_start_y;


		/*
		 * TODO: optimize me! this setup is not really efficient!
		 */
		for (size_t y = 0; y < array_size_y; y++)
			for (size_t x = 0; x < array_size_x; x++)
				data[x + y*array_size_x] = 0;

		if (i_verbosity_level > 5)
			std::cout << "multi-res. map: (" << i_cSimpleNetCDF_DataLevel_parent.array_size_x << ", " << i_cSimpleNetCDF_DataLevel_parent.array_size_y << ") -> (" << array_size_x << ", " << array_size_y << ")" << std::endl;

		for (size_t y = 0; y < i_cSimpleNetCDF_DataLevel_parent.array_size_y; y++)
		{
			size_t local_y = y >> 1;		// ly = y/2

			for (size_t x = 0; x < i_cSimpleNetCDF_DataLevel_parent.array_size_x; x++)
			{
				size_t local_x = x >> 1;		// lx = x/2

				data[local_x + local_y*array_size_x] += i_cSimpleNetCDF_DataLevel_parent.getTexelFloat2D_Clamp(x, y);
			}
		}

		/*
		 * divide non-padded elements by 4
		 */
		for (size_t y = 0; y < non_padded_array_size_y; y++)
			for (size_t x = 0; x < non_padded_array_size_x; x++)
				data[x + y*array_size_x] *= (T)0.25;

		if (array_size_x != non_padded_array_size_x)
		{
			for (size_t y = 0; y < non_padded_array_size_y; y++)
				data[non_padded_array_size_x + y*array_size_x] *= (T)0.5;
		}

		if (array_size_y != non_padded_array_size_y)
		{
			for (size_t x = 0; x < non_padded_array_size_x; x++)
				data[x + non_padded_array_size_y*array_size_x] *= (T)0.5;
		}



		setupDomainParameters(
				i_cSimpleNetCDF_DataLevel_parent.domain_start_x,
				i_cSimpleNetCDF_DataLevel_parent.domain_start_y,
				i_cSimpleNetCDF_DataLevel_parent.domain_length_x,
				i_cSimpleNetCDF_DataLevel_parent.domain_length_y
			);
	}
};


#endif
