/*
 * Copyright (C) 2011 Technische Universitaet Muenchen
 * This file is part of the Sierpi project. For conditions of distribution and
 * use, please see the copyright notice in the file 'copyright.txt' at the root
 * directory of this package and the copyright notice at http://www5.in.tum.de/sierpi
 *
 *  Created on: Feb 20, 2012
 *      Author: Martin Schreiber <martin.schreiber@in.tum.de>
 */

#ifndef CASAGI_HPP_
#define CASAGI_HPP_


#if !CONFIG_ENABLE_MPI
	#define ASAGI_NOMPI
#endif

#include "asagi.h"
#include "../subsimulation_generic/types/CTypes.hpp"
#include "CDataSet_Interface.hpp"


/**
 * dataset interface to asagi
 */
class CAsagi	:
	public CDataSet_Interface
{
private:
	int verbosity_level;

public:
	void getOriginAndSize(
			T	*o_origin_x,
			T	*o_origin_y,
			T	*o_size_x,
			T	*o_size_y
	);


	/**
	 * get nodal data (surface height, momentums and bathymetry) for given coordinate
	 */
	void getNodalData(
			T i_x,
			T i_y,
			T i_level_of_detail,
			CHyperbolicTypes::CSimulationTypes::CNodeData *o_nodal_data
	);

	/**
	 * get dataset data for given dataset id
	 */
	T getDatasetValue(
			int i_dataset_id,
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
			CHyperbolicTypes::CSimulationTypes::CNodeData *o_nodal_data
	);


	/**
	 * get benchmark (analytical) nodal data
	 */
	bool getBenchmarkNodalData(
			T i_x,
			T i_y,
			T i_level_of_detail,
			T i_timestamp,
			CHyperbolicTypes::CSimulationTypes::CNodeData *o_nodal_data
	);


	/**
	 * return zero tolerance
	 */
	T getZeroTolerance();


	/**
	 * return dry tolerance
	 */
	T getDryTolerance();


	/**
	 * return true if dataset is loaded
	 */
	bool isDatasetLoaded();

	/**
	 * output some verbose information
	 */
	void outputVerboseInformation();

	bool loadDatasets(
			const std::vector<std::string> &i_datasets	/// vector with dataset filenames
	);

	void clear();

	CAsagi(int i_verbosity_level);

	virtual ~CAsagi();



private:

	/**
	 * get displacement data for given coordinate
	 */
	T getDisplacementData(
			T i_x,		///< x-coordinate in model-space
			T i_y,		///< x-coordinate in model-space
			T i_level_of_detail		///< level of detail (0 = coarsest level)
	);



	T bathymetry_size_x;
	T bathymetry_size_y;

	T bathymetry_min_x;
	T bathymetry_min_y;

	T bathymetry_max_x;
	T bathymetry_max_y;

	T bathymetry_center_x;
	T bathymetry_center_y;


	T displacements_size_x;
	T displacements_size_y;

	T displacements_min_x;
	T displacements_min_y;

	T displacements_max_x;
	T displacements_max_y;

	T displacements_center_x;
	T displacements_center_y;

	/**
	 * set to true when setup was successful
	 */
	bool is_dataset_loaded;

	/**
	 * singleton for bathymetry datasets
	 */
	asagi::Grid *bathymetry_grid;

	/**
	 * singleton for displacements
	 */
	asagi::Grid *displacements_grid;

};

#endif /* CASAGI_HPP_ */
