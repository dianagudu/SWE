/*
 * Copyright (C) 2011 Technische Universitaet Muenchen
 * This file is part of the Sierpi project. For conditions of distribution and
 * use, please see the copyright notice in the file 'copyright.txt' at the root
 * directory of this package and the copyright notice at http://www5.in.tum.de/sierpi
 *
 *  Created on: Feb 20, 2012
 *      Author: Martin Schreiber <martin.schreiber@in.tum.de>
 */



#include "CAsagi.hpp"
#include <cassert>
#include <stdexcept>


CAsagi::CAsagi(int i_verbosity_level)	:
	verbosity_level(i_verbosity_level),
	bathymetry_grid(nullptr),
	displacements_grid(nullptr),
	is_dataset_loaded(false)
{
}



bool CAsagi::loadDatasets(
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


	/*
	 * BATHYMETRY
	 */
	bathymetry_grid = asagi::Grid::create(asagi::Grid::FLOAT);

	if (bathymetry_grid->open(i_datasets[0].c_str()) != asagi::Grid::SUCCESS)
	{
		std::cerr << "NO ASAGI :-( !" << std::endl;
		std::cerr << "Failed to open bathymetry file '" << i_datasets[0].c_str() << std::endl;
		exit(-1);
	}

	bathymetry_min_x = bathymetry_grid->getXMin();
	bathymetry_min_y = bathymetry_grid->getYMin();

	bathymetry_max_x = bathymetry_grid->getXMax();
	bathymetry_max_y = bathymetry_grid->getYMax();

	bathymetry_size_x = bathymetry_max_x - bathymetry_min_x;
	bathymetry_size_y = bathymetry_max_y - bathymetry_min_y;

	/*
	 * domain size reduced to quad
	 */
	double min_size = CMath::min(bathymetry_size_x, bathymetry_size_y);

	double padding = min_size*0.001;
	double new_size = min_size - 2.0*padding;	// apply padding on both sides

	bathymetry_size_x = new_size;
	bathymetry_size_y = new_size;

	bathymetry_min_x += padding;
	bathymetry_min_y += padding;

	bathymetry_max_x = bathymetry_min_x + bathymetry_size_x;
	bathymetry_max_y = bathymetry_min_y + bathymetry_size_y;

	bathymetry_center_x = ((T)0.5)*(bathymetry_min_x+bathymetry_max_x);
	bathymetry_center_y = ((T)0.5)*(bathymetry_min_y+bathymetry_max_y);

	/**
	 * DISPLACEMENTS
	 */
	displacements_grid = asagi::Grid::create(asagi::Grid::FLOAT);

	if (displacements_grid->open(i_datasets[1].c_str()) != asagi::Grid::SUCCESS)
	{
		std::cerr << "NO ASAGI :-( !" << std::endl;
		std::cerr << "Failed to open file '" << i_datasets[1].c_str() << std::endl;
		is_dataset_loaded = false;
		return false;
	}

	displacements_min_x = displacements_grid->getXMin();
	displacements_min_y = displacements_grid->getYMin();

	displacements_max_x = displacements_grid->getXMax();
	displacements_max_y = displacements_grid->getYMax();

	displacements_size_x = displacements_max_x - displacements_min_x;
	displacements_size_y = displacements_max_y - displacements_min_y;

	displacements_center_x = ((T)0.5)*(displacements_min_x+displacements_max_x);
	displacements_center_y = ((T)0.5)*(displacements_min_y+displacements_max_y);

	is_dataset_loaded = true;

	outputVerboseInformation();
	return true;
}



void CAsagi::clear()
{
	if (bathymetry_grid)
	{
		delete bathymetry_grid;
		bathymetry_grid = nullptr;
	}

	if (displacements_grid)
	{
		delete displacements_grid;
		displacements_grid = nullptr;
	}

	is_dataset_loaded = false;
}



CAsagi::~CAsagi()
{
	clear();
}



void CAsagi::outputVerboseInformation()
{
	std::cout << "ASAGI:" << std::endl;
	std::cout << " + Bathymetry data-window: (" << bathymetry_grid->getXMin() << ", " << bathymetry_grid->getYMin() << 	") x (" << bathymetry_grid->getXMax() << ", " << bathymetry_grid->getYMax() << ")" << std::endl;
	std::cout << " + Bathymetry data-size: (" << (bathymetry_grid->getXMax() - bathymetry_grid->getXMin()) << ", " 	<< (bathymetry_grid->getYMax() - bathymetry_grid->getYMin())<< ")" << std::endl;

	std::cout << " + Used bathymetry start: (" << bathymetry_min_x << ", " << bathymetry_min_y << ")" << std::endl;
	std::cout << " + Used bathymetry size: (" << bathymetry_size_x << ", " << bathymetry_size_y << ")" << std::endl;
	std::cout << " + Used bathymetry end: (" << bathymetry_max_x << ", " << bathymetry_max_y << ")" << std::endl;

	std::cout << " + Displacement window: (" << displacements_grid->getXMin() << ", " << displacements_grid->getYMin() << 	") x (" << displacements_grid->getXMax() << ", " << displacements_grid->getYMax() << ")" << std::endl;
	std::cout << " + Displacement size: (" << (displacements_grid->getXMax() - displacements_grid->getXMin()) << ", " 	<< (displacements_grid->getYMax() - displacements_grid->getYMin())<< ")" << std::endl;
}



void CAsagi::getOriginAndSize(
		T	*o_origin_x,
		T	*o_origin_y,
		T	*o_size_x,
		T	*o_size_y
)
{
	if (bathymetry_grid)
	{
		*o_origin_x = bathymetry_min_x+bathymetry_size_x*(T)0.5;
		*o_origin_y = bathymetry_min_y+bathymetry_size_y*(T)0.5;

		*o_size_x = bathymetry_size_x;
		*o_size_y = bathymetry_size_y;
		return;
	}

	throw(std::runtime_error("No bathymetry dataset loaded!"));
}



CHyperbolicTypes::CSimulationTypes::T CAsagi::getDisplacementData(
		T i_x,		///< x-coordinate in world-space
		T i_y,		///< x-coordinate in world-space
		T i_level_of_detail		///< level of detail (0 = coarsest level)
)
{
	if (i_x < displacements_min_x)
		return 0;

	if (i_x > displacements_max_x)
		return 0;

	if (i_y < displacements_min_y)
		return 0;

	if (i_y > displacements_max_y)
		return 0;

	i_level_of_detail = 0;

	if (displacements_grid)
		return displacements_grid->getFloat2D(i_x, i_y, i_level_of_detail);

	throw(std::runtime_error("No displacements dataset loaded"));

}





/**
 * get nodal data for given coordinate
 */
void CAsagi::getNodalData(
		T i_x,
		T i_y,
		T i_level_of_detail,
		CTsunamiSimulationDOFs *o_nodal_data
)
{
	i_level_of_detail = 0;

	if (bathymetry_grid)
	{
		o_nodal_data->b = bathymetry_grid->getFloat2D(i_x, i_y, i_level_of_detail);
		o_nodal_data->hu = 0;
		o_nodal_data->hv = 0;

		o_nodal_data->h = -o_nodal_data->b + getDisplacementData(i_x, i_y, i_level_of_detail);

		return;
	}

	throw(std::runtime_error("No bathymetry data loaded"));
}


/**
 * get nodal data for given coordinate
 */
CHyperbolicTypes::CSimulationTypes::T CAsagi::getDatasetValue(
		int i_dataset_id,
		T i_x,		///< x-coordinate in model-space
		T i_y,		///< x-coordinate in model-space
		T i_level_of_detail		///< level of detail (0 = coarsest level)
)
{
	i_level_of_detail = 0;

	switch(i_dataset_id)
	{
	case 0:
		if (bathymetry_grid)
			return bathymetry_grid->getFloat2D(i_x, i_y, i_level_of_detail) + getDisplacementData(i_x, i_y, i_level_of_detail);
		throw(std::runtime_error("No bathymetry data loaded"));

	case 1:
		return getDisplacementData(i_x, i_y, i_level_of_detail);

	default:
		throw(std::runtime_error("Invalid dataset"));
	}
}



/**
 * get nodal data for benchmark
 */
bool CAsagi::getBenchmarkNodalData(
		T i_x,
		T i_y,
		T i_level_of_detail,
		T i_timestamp,
		CTsunamiSimulationDOFs *o_nodal_data
)
{
	return false;
}



/**
 * get boundary data at given position
 */
void CAsagi::getBoundaryData(
		T i_x,
		T i_y,
		T i_level_of_detail,
		T i_timestamp,
		CTsunamiSimulationDOFs *o_nodal_data
)
{
	throw(std::runtime_error("No boundary implemented yet"));
}



bool CAsagi::isDatasetLoaded()
{
	return is_dataset_loaded;
}
