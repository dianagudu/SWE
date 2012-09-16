/*
 * Copyright (C) 2011 Technische Universitaet Muenchen
 * This file is part of the Sierpi project. For conditions of distribution and
 * use, please see the copyright notice in the file 'copyright.txt' at the root
 * directory of this package and the copyright notice at http://www5.in.tum.de/sierpi
 *
 *  Created on: July 19, 2012
 *      Author: Martin Schreiber <martin.schreiber@in.tum.de>
 */

#ifndef CSIMPLE_NETCDF_HPP_
#define CSIMPLE_NETCDF_HPP_

#include <vector>
#include <string>

#ifndef CONFIG_COMPILE_WITHOUT_SIERPI
#	include "../subsimulation_generic/types/CTypes.hpp"
#endif

#include "CDataSet_Interface.hpp"


template <typename T>
class CDataSamplingSet_SimpleNetCDF;

/**
 * Dataset interface to the simple netCDF data handler
 *
 * the netCDF data is used in a mipmap
 */
class CSimpleNetCDF	:
	public CDataSet_Interface
{
	int verbosity_level;

	typedef CHyperbolicTypes::CSimulationTypes::CNodeData CNodeData;

public:
	void getOriginAndSize(
			T	*o_origin_x,
			T	*o_origin_y,
			T	*o_size_x,
			T	*o_size_y
	);

	CDataSamplingSet_SimpleNetCDF<T> *cSimpleNetCDF_TerrainData_private;

	CDataSamplingSet_SimpleNetCDF<T> *cSimpleNetCDF_SurfaceDisplacementData_private;

	/**
	 * get nodal data (surface height, momentums and bathymetry) for given coordinate
	 */
	void getNodalData(
			T i_x,
			T i_y,
			T i_level_of_detail,
			CNodeData *o_nodal_data
	);


	/**
	 * get displacement data for given coordinate
	 */
	T getDisplacementData(
			T i_x,		///< x-coordinate in model-space
			T i_y,		///< x-coordinate in model-space
			T i_level_of_detail		///< level of detail (0 = coarsest level)
	);



	/**
	 * get bathymetry data for given coordinate
	 */
	T getDatasetValue(
			int i_dataset_id,	/// id of dataset to return value
			T i_x,		///< x-coordinate in model-space
			T i_y,		///< x-coordinate in model-space
			T i_level_of_detail		///< level of detail (0 = coarsest level)
	);


	/**
	 * get boundary data at given position
	 */
	void getBoundaryData(
			T i_x,
			T i_y,
			T i_level_of_detail,
			T i_timestamp,
			CNodeData *o_nodal_data
	);


	/**
	 * get benchmark (analytical) nodal data
	 */
	bool getBenchmarkNodalData(
			T i_x,
			T i_y,
			T i_level_of_detail,
			T i_timestamp,
			CNodeData *o_nodal_data
	);

	/**
	 * return zero tolerance
	 */
	T getZeroTolerance();

	/**
	 * return dry tolerance
	 */
	T getDryTolerance();

	bool isDatasetLoaded();

	void outputVerboseInformation();

	bool loadDatasets(
			const std::vector<std::string> &i_datasets	/// vector with dataset filenames
	);

	void clear();

	CSimpleNetCDF(
			int i_verbosity_level
	);

	virtual ~CSimpleNetCDF();


private:

	/**
	 * set to true when setup was successful
	 */
	bool is_dataset_loaded;

//public:
	T domain_length_x;
	T domain_length_y;

	T domain_min_x;
	T domain_min_y;

	T domain_max_x;
	T domain_max_y;

	T domain_center_x;
	T domain_center_y;


	T displacements_size_x;
	T displacements_size_y;

	T displacements_min_x;
	T displacements_min_y;

	T displacements_max_x;
	T displacements_max_y;

	T displacements_center_x;
	T displacements_center_y;
};

#endif /* CASAGI_HPP_ */
