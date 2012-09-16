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
#include <cmath>
#include <stdexcept>

#include "CDataSamplingSet_SimpleNetCDF.hpp"
#include "CSimpleNetCDF.hpp"



CSimpleNetCDF::CSimpleNetCDF(
		int i_verbosity_level
)	:
	verbosity_level(i_verbosity_level),
	cSimpleNetCDF_TerrainData_private(nullptr),
	cSimpleNetCDF_SurfaceDisplacementData_private(nullptr),
	is_dataset_loaded(false)
{
}



bool CSimpleNetCDF::loadDatasets(
		const std::vector<std::string> &i_datasets
)
{
	clear();

	if (i_datasets.size() == 0)
		return false;

	if (i_datasets.size() < 2)
	{
		throw(std::runtime_error("Number of datasets != 2"));
	}


	is_dataset_loaded = false;

	cSimpleNetCDF_TerrainData_private = new CDataSamplingSet_SimpleNetCDF<T>(verbosity_level);
	if (cSimpleNetCDF_TerrainData_private->loadDataset(i_datasets[0].c_str()) == false)
	{
		delete cSimpleNetCDF_TerrainData_private;
		cSimpleNetCDF_TerrainData_private = nullptr;
		return false;
	}


	cSimpleNetCDF_SurfaceDisplacementData_private = new CDataSamplingSet_SimpleNetCDF<T>(verbosity_level);
	if (cSimpleNetCDF_SurfaceDisplacementData_private->loadDataset(i_datasets[1].c_str()) == false)
	{
		delete cSimpleNetCDF_TerrainData_private;
		cSimpleNetCDF_TerrainData_private = nullptr;
		return false;
	}

	domain_min_x = cSimpleNetCDF_TerrainData_private->getXMin();
	domain_max_x = cSimpleNetCDF_TerrainData_private->getXMax();

	domain_min_y = cSimpleNetCDF_TerrainData_private->getYMin();
	domain_max_y = cSimpleNetCDF_TerrainData_private->getYMax();

	domain_length_x = domain_max_x - domain_min_x;
	domain_length_y = domain_max_y - domain_min_y;

	domain_center_x = ((T)0.5)*(domain_min_x+domain_max_x);
	domain_center_y = ((T)0.5)*(domain_min_y+domain_max_y);


    if (domain_length_x < domain_length_y)
    {
    	domain_max_y = domain_min_y + domain_length_x;
    	domain_length_y = domain_length_x;
    }
    else
    {
    	domain_max_x = domain_min_x + domain_length_y;
    	domain_length_x = domain_length_y;
    }


	/**
	 * DISPLACEMENTS
	 */
	displacements_min_x = cSimpleNetCDF_SurfaceDisplacementData_private->getXMin();
	displacements_max_x = cSimpleNetCDF_SurfaceDisplacementData_private->getXMax();

	displacements_min_y = cSimpleNetCDF_SurfaceDisplacementData_private->getYMin();
	displacements_max_y = cSimpleNetCDF_SurfaceDisplacementData_private->getYMax();

	displacements_size_x = displacements_max_x - displacements_min_x;
	displacements_size_y = displacements_max_y - displacements_min_y;

	displacements_center_x = ((T)0.5)*(displacements_min_x+displacements_max_x);
	displacements_center_y = ((T)0.5)*(displacements_min_y+displacements_max_y);

	is_dataset_loaded = true;

	outputVerboseInformation();

	return true;
}


void CSimpleNetCDF::clear()
{
	if (cSimpleNetCDF_TerrainData_private)
	{
		delete cSimpleNetCDF_TerrainData_private;
		cSimpleNetCDF_TerrainData_private = nullptr;
	}

	if (cSimpleNetCDF_SurfaceDisplacementData_private)
	{
		delete cSimpleNetCDF_SurfaceDisplacementData_private;
		cSimpleNetCDF_SurfaceDisplacementData_private = nullptr;
	}

	is_dataset_loaded = false;
}


CSimpleNetCDF::~CSimpleNetCDF()
{
	delete cSimpleNetCDF_TerrainData_private;
	delete cSimpleNetCDF_SurfaceDisplacementData_private;
}



void CSimpleNetCDF::outputVerboseInformation()
{
	std::cout << "Simple NetCDF:" << std::endl;

	std::cout << " + Bathymetry data-window: (" << cSimpleNetCDF_TerrainData_private->getXMin() << ", " << cSimpleNetCDF_TerrainData_private->getYMin() << 	") x (" << cSimpleNetCDF_TerrainData_private->getXMax() << ", " << cSimpleNetCDF_TerrainData_private->getYMax() << ")" << std::endl;
	std::cout << " + Bathymetry data-size: (" << (cSimpleNetCDF_TerrainData_private->getXMax() - cSimpleNetCDF_TerrainData_private->getXMin()) << ", " 	<< (cSimpleNetCDF_TerrainData_private->getYMax() - cSimpleNetCDF_TerrainData_private->getYMin())<< ")" << std::endl;

	std::cout << " + Used bathymetry start: (" << domain_min_x << ", " << domain_min_y << ")" << std::endl;
	std::cout << " + Used bathymetry size: (" << domain_length_x << ", " << domain_length_y << ")" << std::endl;
	std::cout << " + Used bathymetry end: (" << domain_max_x << ", " << domain_max_y << ")" << std::endl;

	std::cout << " + Displacement window: (" << cSimpleNetCDF_SurfaceDisplacementData_private->getXMin() << ", " << cSimpleNetCDF_SurfaceDisplacementData_private->getYMin() << 	") x (" << cSimpleNetCDF_SurfaceDisplacementData_private->getXMax() << ", " << cSimpleNetCDF_SurfaceDisplacementData_private->getYMax() << ")" << std::endl;
	std::cout << " + Displacement size: (" << (cSimpleNetCDF_SurfaceDisplacementData_private->getXMax() - cSimpleNetCDF_SurfaceDisplacementData_private->getXMin()) << ", " 	<< (cSimpleNetCDF_SurfaceDisplacementData_private->getYMax() - cSimpleNetCDF_SurfaceDisplacementData_private->getYMin())<< ")" << std::endl;

	std::cout << " + Terrain raw-data array size: " << cSimpleNetCDF_TerrainData_private->getArraySizeX() << ", " << cSimpleNetCDF_TerrainData_private->getArraySizeY() << std::endl;

	size_t min_array_size = std::min(cSimpleNetCDF_TerrainData_private->getArraySizeX(), cSimpleNetCDF_TerrainData_private->getArraySizeY());

	T log2_array_size = std::log(static_cast<T>(min_array_size)) / std::log((T)2.0);
	std::cout << " + Max meaningful cartesian-grid refinement depth: " << log2_array_size << std::endl;
	std::cout << " + Max meaningful triangle-grid refinement depth: " << (log2_array_size*2) << std::endl;
	std::cout << " + Terrain raw-data cell size: (" << (cSimpleNetCDF_TerrainData_private->getXMax() - cSimpleNetCDF_TerrainData_private->getXMin())/(T)cSimpleNetCDF_TerrainData_private->getArraySizeX() << ", " << (cSimpleNetCDF_TerrainData_private->getYMax() - cSimpleNetCDF_TerrainData_private->getYMin())/(T)cSimpleNetCDF_TerrainData_private->getArraySizeY() << ")" << std::endl;
}



void CSimpleNetCDF::getOriginAndSize(
		T	*o_origin_x,
		T	*o_origin_y,
		T	*o_size_x,
		T	*o_size_y
)
{
	if (cSimpleNetCDF_TerrainData_private)
	{
		*o_origin_x = domain_min_x+domain_length_x*(T)0.5;
		*o_origin_y = domain_min_y+domain_length_y*(T)0.5;

		*o_size_x = domain_length_x;
		*o_size_y = domain_length_y;
		return;
	}

	throw(std::runtime_error("No bathymetry dataset loaded!"));
}





/**
 * get nodal data for given coordinate
 */
void CSimpleNetCDF::getNodalData(
		T i_x,
		T i_y,
		T i_level_of_detail,
		CHyperbolicTypes::CSimulationTypes::CNodeData *o_nodal_data
)
{
	if (cSimpleNetCDF_SurfaceDisplacementData_private)
	{
		o_nodal_data->b = cSimpleNetCDF_TerrainData_private->getFloat2D(i_x, i_y, i_level_of_detail);
		o_nodal_data->hu = 0;
		o_nodal_data->hv = 0;
		o_nodal_data->h = -o_nodal_data->b + cSimpleNetCDF_SurfaceDisplacementData_private->getFloat2D(i_x, i_y, i_level_of_detail);
		return;
	}
	
	throw(std::runtime_error("No bathymetry data loaded"));
}


/**
 * get dataset value for given coordinate and dataset id
 */
CHyperbolicTypes::CSimulationTypes::T CSimpleNetCDF::getDatasetValue(
		int i_dataset_id,	///< id of dataset to get value from
		T i_x,		///< x-coordinate in model-space
		T i_y,		///< x-coordinate in model-space
		T i_level_of_detail		///< level of detail (0 = coarsest level)
)
{
	if (cSimpleNetCDF_TerrainData_private)
	{
		switch(i_dataset_id)
		{
		case 0:
			return cSimpleNetCDF_TerrainData_private->getFloat2D(i_x, i_y, i_level_of_detail);

		case 1:
			return cSimpleNetCDF_SurfaceDisplacementData_private->getFloat2D(i_x, i_y, i_level_of_detail);

		default:
			throw(std::runtime_error("Invalid dataset id"));
		}
	}

	throw(std::runtime_error("No bathymetry data loaded"));
}



/**
 * get nodal data for benchmark
 */
bool CSimpleNetCDF::getBenchmarkNodalData(
		T i_x,
		T i_y,
		T i_level_of_detail,
		T i_timestamp,
		CHyperbolicTypes::CSimulationTypes::CNodeData *o_nodal_data
)
{
	return false;
}


/**
 * get boundary data at given position
 */
void CSimpleNetCDF::getBoundaryData(
		T i_x,
		T i_y,
		T i_level_of_detail,
		T i_timestamp,
		CHyperbolicTypes::CSimulationTypes::CNodeData *o_nodal_data
)
{
	throw(std::runtime_error("No boundary implemented yet"));
}



bool CSimpleNetCDF::isDatasetLoaded()
{
	return is_dataset_loaded;
}
