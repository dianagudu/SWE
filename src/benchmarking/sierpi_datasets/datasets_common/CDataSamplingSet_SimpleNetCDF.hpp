/*
 * Copyright (C) 2011 Technische Universitaet Muenchen
 * This file is part of the Sierpi project. For conditions of distribution and
 * use, please see the copyright notice in the file 'copyright.txt' at the root
 * directory of this package and the copyright notice at http://www5.in.tum.de/sierpi
 *
 *  Created on: July 19, 2012
 *      Author: Martin Schreiber <martin.schreiber@in.tum.de>
 */

#include <cassert>
#include <stdexcept>
#include <netcdf.hh>
#include "CDataSamplingSet_MultiResolution.hpp"



/**
 * Simple NetCDF datafile loader
 *
 * We assume pixel-node registration, see "B.2.2 Gridline and Pixel node registration" of GMT Docs
 */
template <typename T>
class CDataSamplingSet_SimpleNetCDF	:
	public CDataSamplingSet_MultiResolution<T>
{
		using CDataSamplingSet_MultiResolution<T>::number_of_data_levels;
		using CDataSamplingSet_MultiResolution<T>::dataLevels;

public:

	void setOutOfDomainValue(T i_value)
	{
		for (int i = 0; i < number_of_data_levels; i++)
		{
			dataLevels[i].out_of_domain_value = i_value;
		}
	}

	CDataSamplingSet_SimpleNetCDF(
			int i_verbosity_level
	)	:
		CDataSamplingSet_MultiResolution<T>(i_verbosity_level)
	{

	}

	size_t getXResolution()
	{
		return CDataSamplingSet_MultiResolution<T>::getXResolution();
	}

	size_t getYResolution()
	{
		return CDataSamplingSet_MultiResolution<T>::getYResolution();
	}

	T getXMin()
	{
		return CDataSamplingSet_MultiResolution<T>::getXMin();
	}

	T getXMax()
	{
		return CDataSamplingSet_MultiResolution<T>::getXMax();
	}

	T getYMin()
	{
		return CDataSamplingSet_MultiResolution<T>::getYMin();
	}

	T getYMax()
	{
		return CDataSamplingSet_MultiResolution<T>::getYMax();
	}

	T getFloat2D(
			T i_x,
			T i_y,
			int i_level
	)
	{
		assert(i_level >= 0);
		return CDataSamplingSet_MultiResolution<T>::getFloat2D(i_x, i_y, i_level);
	}


	T getFloat2D_BorderValue(
			T i_x,
			T i_y,
			int i_level
	)
	{
		assert(i_level >= 0);
		return CDataSamplingSet_MultiResolution<T>::getFloat2D_BorderValue(i_x, i_y, i_level);
	}


	~CDataSamplingSet_SimpleNetCDF()
	{
	}



	/**
	 * read range attribute from given variable
	 */
private:
	bool p_getRangeFromVariable(
			int ncid,
			int varid,
			T *o_min,
			T *o_max
	)
	{
		int retval;
		size_t len;

		if ((retval = nc_inq_attlen(ncid, varid, "actual_range", &len)) == -1)
		{
			throw(std::runtime_error(std::string("nc_inq_attlen 'actual_range' not found: ")+nc_strerror(retval)));
			return false;
		}

		T *range = new T[len];


		if (sizeof(T) == 4)
		{
			if ((retval = nc_get_att_float(ncid, varid, "actual_range", (float*)&range[0])) == -1)
			{
				throw(std::runtime_error(std::string("nc_get_att_float: ")+nc_strerror(retval)));
				return false;
			}
		}
		else
		{
			if ((retval = nc_get_att_double(ncid, varid, "actual_range", (double*)&range[0])) == -1)
			{
				throw(std::runtime_error(std::string("nc_get_att_double: ")+nc_strerror(retval)));
				return false;
			}
		}

		*o_min = range[0];
		*o_max = range[1];

		delete [] range;
		return true;
	}



	/**
	 * load the dataset from a netcdf file
	 *
	 * \return true if everything is ok
	 */
public:
	bool loadDataset(
			const char *i_dataset_filename
	)
	{
		CDataSamplingSet_MultiResolution<T>::clear();

		/*
		 * NETCDF LOADER
		 */
		int retval;
		int ncid;

		if ((retval = nc_open(i_dataset_filename, NC_NOWRITE, &ncid)) != NC_NOERR)
		{
			throw(std::runtime_error(std::string("nc_open ")+i_dataset_filename));
		}

		int ndims, nvars, ngatts, recdim;
		if ((retval = ncinquire(ncid, &ndims, &nvars, &ngatts, &recdim)) == -1)
		{
		    nc_close(ncid);
			throw(std::runtime_error(std::string("ncinquire ")+i_dataset_filename));
		}

		if (ndims != 2)
		{
		    nc_close(ncid);
			throw(std::runtime_error(std::string("Expected 2 dimensions")));
		}

		char char_buffer[MAX_NC_NAME];
		long size;

		if ((retval = ncdiminq(ncid, 0, &char_buffer[0], &size)) == -1)
		{
			nc_close(ncid);
			throw(std::runtime_error(std::string("ERROR ncdiminq: ")+nc_strerror(retval)));
		}
		size_t array_size_x = size;

		if ((retval = ncdiminq(ncid, 1, &char_buffer[0], &size)) == -1)
		{
			nc_close(ncid);
			throw(std::runtime_error(std::string("ERROR ncdiminq: ")+nc_strerror(retval)));
		}
		size_t array_size_y = size;


		int var_x_id;
		if ((retval = nc_inq_varid(ncid, "x", &var_x_id)) == -1)
		{
			nc_close(ncid);
			throw(std::runtime_error(std::string("ERROR var x not found: ")+nc_strerror(retval)));
		}

		T xmin, xmax;
		try		{
			p_getRangeFromVariable(ncid, var_x_id, &xmin, &xmax);
		}
		catch(std::runtime_error &r)		{
			nc_close(ncid);
			throw(r);
		}

		int var_y_id;
		if ((retval = nc_inq_varid(ncid, "y", &var_y_id)) == -1)
		{
			nc_close(ncid);
			throw(std::runtime_error(std::string("ERROR var y not found: ")+nc_strerror(retval)));
		}

		T ymin, ymax;
		try		{
			p_getRangeFromVariable(ncid, var_y_id, &ymin, &ymax);
		}
		catch(std::runtime_error &r)		{
			nc_close(ncid);
			throw(r);
		}

		int var_z_id;
		if ((retval = nc_inq_varid(ncid, "z", &var_z_id)) == -1)
		{
			nc_close(ncid);
			throw(std::runtime_error(std::string("ERROR var z not found: ")+nc_strerror(retval)));
		}


		/*
		 * allocate first level of multiresolution map
		 *
		 * the return value is the pointer to the first level data
		 */
		T *data = CDataSamplingSet_MultiResolution<T>::allocateMultiresolutionLevels(array_size_x, array_size_y);

		if (sizeof(T) == 4)
			retval =  nc_get_var_float(ncid, var_z_id, (float*)data);
		else
			retval =  nc_get_var_double(ncid, var_z_id, (double*)data);

		if (retval == -1)
		{
			nc_close(ncid);
			throw(std::runtime_error(std::string("ERROR while loading variable data: ")+nc_strerror(retval)));
		}

		nc_close(ncid);

		/*
		 * setup datasets
		 */
		CDataSamplingSet_MultiResolution<T>::setupDomainParameters(xmin, ymin, xmax, ymax);

		/*
		 * setup the multiresolution map
		 */
		CDataSamplingSet_MultiResolution<T>::setupMultiresolutionMap();

		return true;
	}
};
