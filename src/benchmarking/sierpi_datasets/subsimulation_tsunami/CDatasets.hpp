/*
 * Copyright (C) 2011 Technische Universitaet Muenchen
 * This file is part of the Sierpi project. For conditions of distribution and
 * use, please see the copyright notice in the file 'copyright.txt' at the root
 * directory of this package and the copyright notice at http://www5.in.tum.de/sierpi
 *
 *  Created on: Mar 16, 2012
 *      Author: Martin Schreiber <martin.schreiber@in.tum.de>
 */

#ifndef CTSUNAMISIMULATION_DATASETS_HPP_
#define CTSUNAMISIMULATION_DATASETS_HPP_

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <cassert>

// only include when compiled with sierpi
#ifndef CONFIG_COMPILE_WITHOUT_SIERPI
#	include "../CParameters.hpp"
#	include "../datasets_common/CParameters_Datasets.hpp"
#endif


#if CONFIG_ENABLE_ASAGI
#	include "../datasets_common/CAsagi.hpp"
#endif

#if CONFIG_ENABLE_NETCDF
#	include "../datasets_common/CSimpleNetCDF.hpp"
#	include "tsunami_benchmarks/CSingleWaveOnSimpleBeach.hpp"
# 	include "tsunami_benchmarks/CSingleWaveOnCompositeBeach.hpp"
#endif

#include "../datasets_common/CDataSet_Interface.hpp"



class CDatasets	:
	public CHyperbolicTypes::CSimulationTypes,	///< import tsunami simulation types
	public CDataSet_Interface			///< this class has to provide the same dataset interface as each dataset
{
	typedef CHyperbolicTypes::CSimulationTypes::T T;

public:

	/**
	 * possible choices to setup water height
	 */
	enum
	{
		/*
		 * use simple netcdf
		 */
		SIMULATION_WATER_HEIGHT_SIMPLE_NETCDF = -2,

		/*
		 * use asagi to determine water surface height
		 */
		SIMULATION_WATER_HEIGHT_ASAGI = -1,

		/*
		 * set simulation water surface height to zero
		 */
		SIMULATION_WATER_HEIGHT_DEFAULT = 0,

		/*
		 * setup water surface height with cylinder parameters (simulation_dataset_cylinder*)
		 */
		SIMULATION_WATER_HEIGHT_CYLINDER = 1,
		SIMULATION_WATER_HEIGHT_LINEAR_SLOPE = 2,
		SIMULATION_WATER_HEIGHT_SQUARE = 3,
		SIMULATION_WATER_HEIGHT_SQUARE_DIST = 4,
		SIMULATION_WATER_HEIGHT_PARABOLA = 5,
		SIMULATION_WATER_HEIGHT_CYLINDER_RELATIVE_TO_BATHYMETRY = 6,
		SIMULATION_WATER_HEIGHT_RADIAL_DAM_BREAK_OUTER_AREA_MINUS_INF = 7,
		SIMULATION_WATER_HEIGHT_RADIAL_DAM_BREAK = 8,
		SIMULATION_WATER_HEIGHT_SINGLE_WAVE_ON_SINGLE_BEACH = 100,
		SIMULATION_WATER_HEIGHT_SINGLE_WAVE_ON_COMPOSITE_BEACH = 101,


		SIMULATION_INTERACTIVE_UPDATE = 9
	};


	/**
	 * possible choices to setup terrain distance
	 */
	enum
	{
		SIMULATION_TERRAIN_HEIGHT_SIMPLE_NETCDF = -2,
		SIMULATION_TERRAIN_HEIGHT_ASAGI = -1,

		SIMULATION_TERRAIN_HEIGHT_DEFAULT = 0,

		SIMULATION_TERRAIN_HEIGHT_LINEAR_SLOPE = 1,
		SIMULATION_TERRAIN_HEIGHT_PARABOLA = 2,
		SIMULATION_TERRAIN_HEIGHT_FANCY = 3,
		SIMULATION_TERRAIN_HEIGHT_SINGLE_WAVE_ON_SINGLE_BEACH = 100,
		SIMULATION_TERRAIN_HEIGHT_SINGLE_WAVE_ON_COMPOSITE_BEACH = 101
	};


	/**
	 * simulation parameters for tsunami simulation
	 */
	CParameters_Datasets &cParameters_Datasets;

	/**
	 * ASAGI dataset interface
	 */
#if CONFIG_ENABLE_ASAGI
	CAsagi cAsagi;
#endif

	/**
	 * simple netcdf interface
	 */
#if CONFIG_ENABLE_NETCDF
	CSimpleNetCDF cSimpleNetCDF;

	/**
	 * analytical benchmark: single wave on single beach
	 */
	CSingleWaveOnSimpleBeach<T> cSingleWaveOnSingleBeach;
	
	/**
	 * analytical benchmark: single wave on composite beach
	 */
	CSingleWaveOnCompositeBeach<T> cSingleWaveOnCompositeBeach;
#endif

	/**
	 * verbosity level
	 */
	int verbosity_level;



public:
	/**
	 * constructor
	 */
	CDatasets(
			CParameters_Datasets &i_cParameters_Datasets,
			int i_verbosity_level
	)	:
		cParameters_Datasets(i_cParameters_Datasets),
#if CONFIG_ENABLE_ASAGI
		cAsagi(i_verbosity_level),
#endif
#if CONFIG_ENABLE_NETCDF
		cSimpleNetCDF(i_verbosity_level),
#endif
		verbosity_level(i_verbosity_level)
	{
	}


	void loadDatasets()
	{
		loadDatasets(cParameters_Datasets.simulation_datasets);
	}

	/**
	 * \brief setup datasets
	 *
	 * This is important for the datafile based simulation scenarios
	 */
	bool loadDatasets(
			const std::vector<std::string> &i_datasets
	)
	{
		switch(cParameters_Datasets.simulation_dataset_0_id)
		{
#if CONFIG_ENABLE_ASAGI
		case SIMULATION_TERRAIN_HEIGHT_ASAGI:
			if (	i_datasets[0].length() > 0		&&
					i_datasets[1].length() > 0
			)
			{
				cAsagi.loadDatasets(i_datasets);


				if (!cAsagi.isDatasetLoaded())
					return false;

				cAsagi.getOriginAndSize(
						&cParameters_Datasets.simulation_dataset_default_domain_origin_x,
						&cParameters_Datasets.simulation_dataset_default_domain_origin_y,
						&cParameters_Datasets.simulation_dataset_default_domain_size_x,
						&cParameters_Datasets.simulation_dataset_default_domain_size_y
				);

				return true;
			}
			break;
#endif

#if CONFIG_ENABLE_NETCDF
		case SIMULATION_TERRAIN_HEIGHT_SIMPLE_NETCDF:
			if (	i_datasets[0].length() > 0		&&
					i_datasets[1].length() > 0
			)
			{
				cSimpleNetCDF.loadDatasets(i_datasets);

				if (!cSimpleNetCDF.isDatasetLoaded())
					return false;

				cSimpleNetCDF.getOriginAndSize(
						&cParameters_Datasets.simulation_dataset_default_domain_origin_x,
						&cParameters_Datasets.simulation_dataset_default_domain_origin_y,
						&cParameters_Datasets.simulation_dataset_default_domain_size_x,
						&cParameters_Datasets.simulation_dataset_default_domain_size_y
				);

				return true;
			}
			break;
#endif
		}

		return true;
	}



public:
	/**
	 * deconstructor
	 */
	virtual ~CDatasets()
	{
	}




public:
	/**
	 * terrain dimensions
	 */
	inline void getOriginAndSize(
			T *o_origin_x,	///< origins x-coordinate
			T *o_origin_y,	///< origins y-coordinate
			T *o_size_x,	///< size of terrain in x-dimension
			T *o_size_y		///< size of terrain in y-dimension
	)
	{
		switch(cParameters_Datasets.simulation_dataset_0_id)
		{
#if CONFIG_ENABLE_ASAGI
		case SIMULATION_TERRAIN_HEIGHT_ASAGI:
			cAsagi.getOriginAndSize(o_origin_x, o_origin_y, o_size_x, o_size_y);
			break;
#endif

#if CONFIG_ENABLE_NETCDF
		case SIMULATION_TERRAIN_HEIGHT_SIMPLE_NETCDF:
			cSimpleNetCDF.getOriginAndSize(o_origin_x, o_origin_y, o_size_x, o_size_y);
			break;

#endif

		default:
			*o_origin_x = 0;
			*o_origin_y = 0;

			*o_size_x = cParameters_Datasets.simulation_dataset_default_domain_size_x;
			*o_size_y = cParameters_Datasets.simulation_dataset_default_domain_size_y;
			break;
		}
	}



	/**
	 * store benchmarks nodal data for given timestamp to o_nodal_data
	 */
	inline bool getBenchmarkNodalData(
			T i_x,		///< x-coordinate in world-space
			T i_y,		///< y-coordinate in world-space
			T i_level_of_detail,	///< level of detail (0 = coarsest level)
			T i_timestamp,			///< timestamp for boundary data
			CHyperbolicTypes::CSimulationTypes::CNodeData *o_cNodeData	///< nodal data is written to this position
	)
	{
		switch(cParameters_Datasets.simulation_dataset_0_id)
		{
#if CONFIG_ENABLE_ASAGI
		case SIMULATION_TERRAIN_HEIGHT_ASAGI:
			return cAsagi.getBenchmarkNodalData(i_x, i_y, i_level_of_detail, i_timestamp, o_cNodeData);
#endif

#if CONFIG_ENABLE_NETCDF
		case SIMULATION_TERRAIN_HEIGHT_SIMPLE_NETCDF:
			return cSimpleNetCDF.getBenchmarkNodalData(i_x, i_y, i_level_of_detail, i_timestamp, o_cNodeData);


		case SIMULATION_WATER_HEIGHT_SINGLE_WAVE_ON_SINGLE_BEACH:
			return cSingleWaveOnSingleBeach.getBenchmarkNodalData(i_x, i_y, i_level_of_detail, i_timestamp, o_cNodeData);
			
		case SIMULATION_WATER_HEIGHT_SINGLE_WAVE_ON_COMPOSITE_BEACH:
			cSingleWaveOnCompositeBeach.getBenchmarkNodalData(i_x, i_y, i_level_of_detail, i_timestamp, o_cNodeData);
			break;
#endif
		}
		return false;
	}


	inline void getNodalData(
			T i_x,					///< x-coordinate in world-space
			T i_y,					///< y-coordinate in world-space
			T i_level_of_detail,	///< level of detail (0 = coarsest level)
			CHyperbolicTypes::CSimulationTypes::CNodeData *o_cNodeData	///< nodal data is written to this position
	)
	{
		if (cParameters_Datasets.simulation_dataset_1_id == SIMULATION_INTERACTIVE_UPDATE)
		{
			T rx = (i_x - cParameters_Datasets.simulation_dataset_breaking_dam_posx);
			T ry = (i_y - cParameters_Datasets.simulation_dataset_breaking_dam_posy);

			if (rx*rx + ry*ry >= cParameters_Datasets.simulation_dataset_breaking_dam_radius*cParameters_Datasets.simulation_dataset_breaking_dam_radius)
				return;

			T horizon = o_cNodeData->b+o_cNodeData->h;
			if (horizon < cParameters_Datasets.simulation_dataset_default_values[1])
				o_cNodeData->h = -o_cNodeData->b+cParameters_Datasets.simulation_dataset_default_values[1];

			return;
		}

		switch(cParameters_Datasets.simulation_dataset_0_id)
		{
#if CONFIG_ENABLE_NETCDF
		case SIMULATION_WATER_HEIGHT_SINGLE_WAVE_ON_SINGLE_BEACH :
			cSingleWaveOnSingleBeach.getNodalData(i_x,i_y,i_level_of_detail,o_cNodeData);
			break;

		case SIMULATION_WATER_HEIGHT_SINGLE_WAVE_ON_COMPOSITE_BEACH:
			cSingleWaveOnCompositeBeach.getNodalData(i_x,i_y,i_level_of_detail,o_cNodeData);
			break;
#endif

		default:
			o_cNodeData->b = getDatasetValue(0, i_x, i_y, i_level_of_detail);
			o_cNodeData->h = std::max<T>(getDatasetValue(1, i_x, i_y, i_level_of_detail) - o_cNodeData->b, 0);
			o_cNodeData->hu = 0;
			o_cNodeData->hv = 0;

			break;
		}
	}



public:
	/**
	 * terrain setup
	 */
	inline T getDatasetValue(
			int i_dataset_id,
			T i_x,
			T i_y,
			T i_level_of_detail
	)
	{
		if (i_dataset_id == 0)
		{
			/*
			 * bathymetry
			 */

			switch(cParameters_Datasets.simulation_dataset_0_id)
			{
#if CONFIG_ENABLE_ASAGI
				case SIMULATION_TERRAIN_HEIGHT_ASAGI:
					return cAsagi.getDatasetValue(0, i_x, i_y, i_level_of_detail);
#endif

#if CONFIG_ENABLE_NETCDF
				case SIMULATION_TERRAIN_HEIGHT_SIMPLE_NETCDF:
					return cSimpleNetCDF.getDatasetValue(0, i_x, i_y, i_level_of_detail);

				case SIMULATION_TERRAIN_HEIGHT_SINGLE_WAVE_ON_SINGLE_BEACH:
					return cSingleWaveOnSingleBeach.getDatasetValue(0, i_x, i_y, i_level_of_detail);
					
				case SIMULATION_TERRAIN_HEIGHT_SINGLE_WAVE_ON_COMPOSITE_BEACH:
					return cSingleWaveOnCompositeBeach.getDatasetValue(0, i_x,i_y,i_level_of_detail);
#endif

				case SIMULATION_TERRAIN_HEIGHT_DEFAULT:
					return -cParameters_Datasets.simulation_dataset_default_values[0];

				case SIMULATION_TERRAIN_HEIGHT_LINEAR_SLOPE:
					return (i_x-cParameters_Datasets.simulation_dataset_default_domain_origin_x)/(cParameters_Datasets.simulation_dataset_default_domain_size_x*0.5)*cParameters_Datasets.simulation_dataset_default_values[0];

				case SIMULATION_TERRAIN_HEIGHT_PARABOLA:
					return -(T)2.0*((T)1.0-std::max(i_x*i_x, i_y*i_y))*cParameters_Datasets.simulation_dataset_default_values[0] +
							cParameters_Datasets.simulation_dataset_default_values[0]*(T)0.3;

				case SIMULATION_TERRAIN_HEIGHT_FANCY:
				{
					T d = (T)-1.0;

					d = std::max(d, i_x*i_x*(T)0.9-(T)0.8);

					if (i_x*i_x + i_y*i_y < (T)0.05)
						d = (T)0.15;

					T nx = i_x-(T)1.0;
					d = std::max(d, nx*nx*nx*(T)0.3+nx*nx*(T)1.0+nx+(T)0.29);

					d -= 0.1;
					return d*cParameters_Datasets.simulation_dataset_default_values[0]*0.01;
				}

				default:
					return -cParameters_Datasets.simulation_dataset_default_values[0];
			}
		}


		if (i_dataset_id == 1)
		{
			/*
			 * displacements
			 */

			switch(cParameters_Datasets.simulation_dataset_1_id)
			{
#if CONFIG_ENABLE_ASAGI
				case SIMULATION_WATER_HEIGHT_ASAGI:		// setup with asagi
					return cAsagi.getDatasetValue(1, i_x, i_y, i_level_of_detail);
#endif

#if CONFIG_ENABLE_NETCDF
				case SIMULATION_WATER_HEIGHT_SIMPLE_NETCDF:
					return cSimpleNetCDF.getDatasetValue(1, i_x, i_y, i_level_of_detail);
#endif

				case SIMULATION_WATER_HEIGHT_DEFAULT:
					return cParameters_Datasets.simulation_dataset_default_values[1];

				case SIMULATION_WATER_HEIGHT_CYLINDER:
				{
					CHyperbolicTypes::CSimulationTypes::T rx = (i_x - cParameters_Datasets.simulation_dataset_breaking_dam_posx);
					CHyperbolicTypes::CSimulationTypes::T ry = (i_y - cParameters_Datasets.simulation_dataset_breaking_dam_posy);

					if (rx*rx + ry*ry < cParameters_Datasets.simulation_dataset_breaking_dam_radius*cParameters_Datasets.simulation_dataset_breaking_dam_radius)
						return cParameters_Datasets.simulation_dataset_breaking_dam_height;
					else
						return 0;
				}

				case SIMULATION_WATER_HEIGHT_LINEAR_SLOPE:
					return (i_x-cParameters_Datasets.simulation_dataset_default_domain_origin_x)/(cParameters_Datasets.simulation_dataset_default_domain_size_x*0.5)*cParameters_Datasets.simulation_dataset_default_values[0];

				case SIMULATION_WATER_HEIGHT_SQUARE:
				{
					if (i_x > -1.0+1.0/std::pow(2.0, 2.0)+0.0001)
						return 0;
					return cParameters_Datasets.simulation_dataset_breaking_dam_height;
				}

				case SIMULATION_WATER_HEIGHT_SQUARE_DIST:
				{
					if (i_x < -1.0+1.0/std::pow(2.0, 4.0)+0.0001 || i_x > -1.0+1.0/std::pow(2.0, 2.0)+0.0001)
						return 0;
					return cParameters_Datasets.simulation_dataset_breaking_dam_height;
				}

				case SIMULATION_WATER_HEIGHT_PARABOLA:
					return 1.0-std::max(i_x*i_x, i_y*i_y);


					/*
					 * setup a cylinder relative to bathymetry if bathymetry is higher than 0
					 */
				case SIMULATION_WATER_HEIGHT_CYLINDER_RELATIVE_TO_BATHYMETRY:
				{
					CHyperbolicTypes::CSimulationTypes::T rx = (i_x - cParameters_Datasets.simulation_dataset_breaking_dam_posx);
					CHyperbolicTypes::CSimulationTypes::T ry = (i_y - cParameters_Datasets.simulation_dataset_breaking_dam_posy);

					if (rx*rx + ry*ry > cParameters_Datasets.simulation_dataset_breaking_dam_radius*cParameters_Datasets.simulation_dataset_breaking_dam_radius)
						return -99999999.0;	// return huge negative value since cylinder is set-up by using max() function with currently existing value and return value of this method

					return cParameters_Datasets.simulation_dataset_breaking_dam_height;
				}

				case SIMULATION_WATER_HEIGHT_RADIAL_DAM_BREAK_OUTER_AREA_MINUS_INF:
				{
					CHyperbolicTypes::CSimulationTypes::T rx = (i_x - cParameters_Datasets.simulation_dataset_breaking_dam_posx);
					CHyperbolicTypes::CSimulationTypes::T ry = (i_y - cParameters_Datasets.simulation_dataset_breaking_dam_posy);

					if (rx*rx + ry*ry < cParameters_Datasets.simulation_dataset_breaking_dam_radius*cParameters_Datasets.simulation_dataset_breaking_dam_radius)
						return cParameters_Datasets.simulation_dataset_breaking_dam_height;
					else
						return -99999999.0;	// return huge negative value since cylinder is set-up by using max() function with currently existing value and return value of this method
				}


				case SIMULATION_WATER_HEIGHT_RADIAL_DAM_BREAK:
				{
					CHyperbolicTypes::CSimulationTypes::T rx = (i_x - cParameters_Datasets.simulation_dataset_breaking_dam_posx);
					CHyperbolicTypes::CSimulationTypes::T ry = (i_y - cParameters_Datasets.simulation_dataset_breaking_dam_posy);

					if (rx*rx + ry*ry < cParameters_Datasets.simulation_dataset_breaking_dam_radius*cParameters_Datasets.simulation_dataset_breaking_dam_radius)
						return cParameters_Datasets.simulation_dataset_breaking_dam_height;
					else
						return 0;
				}

#if CONFIG_ENABLE_NETCDF
				case SIMULATION_WATER_HEIGHT_SINGLE_WAVE_ON_SINGLE_BEACH:
					return cSingleWaveOnSingleBeach.getDatasetValue(1, i_x, i_y, i_level_of_detail);
				
				case SIMULATION_WATER_HEIGHT_SINGLE_WAVE_ON_COMPOSITE_BEACH:
					return cSingleWaveOnCompositeBeach.getDatasetValue(1, i_x,i_y,i_level_of_detail);
#endif

				default:
					return cParameters_Datasets.simulation_dataset_1_id;
			}
		}

		throw(std::runtime_error("unknown dataset id"));
		return 0;
	}


	/**
	 * get boundary data at given position
	 */
	void getBoundaryData(
			T i_x,					///< x-coordinate in world-space
			T i_y,					///< y-coordinate in world-space
			T i_level_of_detail,	///< level of detail (0 = coarsest level)
			T i_timestamp,			///< timestamp for boundary data
			CHyperbolicTypes::CSimulationTypes::CNodeData *o_nodal_data	///< nodal data is written to this position
	)
	{
#if CONFIG_ENABLE_NETCDF
		switch(cParameters_Datasets.simulation_dataset_1_id)
		{
			case SIMULATION_WATER_HEIGHT_SINGLE_WAVE_ON_COMPOSITE_BEACH:
				cSingleWaveOnCompositeBeach.getBoundaryData(i_x,i_y,i_timestamp,i_level_of_detail,o_nodal_data);
				break;
		}
#endif
		std::cerr << "getBoundaryData not implemented yet" << std::endl;
		assert(false);
		return;
	}


	/**
	 * return true if the dataset is loaded
	 */
	bool isDatasetLoaded()
	{
		switch(cParameters_Datasets.simulation_dataset_1_id)
		{
#if CONFIG_ENABLE_ASAGI
		case SIMULATION_TERRAIN_HEIGHT_ASAGI:		// setup with asagi
			return cAsagi.isDatasetLoaded();
#endif

#if CONFIG_ENABLE_NETCDF
		case SIMULATION_TERRAIN_HEIGHT_SIMPLE_NETCDF:
			return cSimpleNetCDF.isDatasetLoaded();
#endif

		default:
			return true;
		}
	}


	/**
	 * output verbose information
	 */
	void outputVerboseInformation()
	{
		switch(cParameters_Datasets.simulation_dataset_1_id)
		{
#if CONFIG_ENABLE_ASAGI
		case SIMULATION_WATER_HEIGHT_ASAGI:		// setup with asagi
			cAsagi.outputVerboseInformation();
			break;
#endif

#if CONFIG_ENABLE_NETCDF
		case SIMULATION_WATER_HEIGHT_SIMPLE_NETCDF:
			cSimpleNetCDF.outputVerboseInformation();
			break;
#endif


		case SIMULATION_WATER_HEIGHT_DEFAULT:
			std::cout << "Water setup: default displacement" << std::endl;
			break;

		case SIMULATION_WATER_HEIGHT_CYLINDER:
			std::cout << "Water setup: cylinder" << std::endl;
			break;

		case SIMULATION_WATER_HEIGHT_LINEAR_SLOPE:
			std::cout << "Water setup: linear slope" << std::endl;
			break;

		case SIMULATION_WATER_HEIGHT_SQUARE:
			std::cout << "Water setup: square" << std::endl;
			break;

		case SIMULATION_WATER_HEIGHT_SQUARE_DIST:
			std::cout << "Water setup: square dist" << std::endl;
			break;

		case SIMULATION_WATER_HEIGHT_PARABOLA:
			std::cout << "Water setup: parabola" << std::endl;
			break;

		case SIMULATION_WATER_HEIGHT_CYLINDER_RELATIVE_TO_BATHYMETRY:
			std::cout << "Water setup: cylinder relative to bathymetry" << std::endl;
			break;

		case SIMULATION_WATER_HEIGHT_RADIAL_DAM_BREAK_OUTER_AREA_MINUS_INF:
			std::cout << "Water setup: radial dam break outer area minus inf" << std::endl;
			break;

		case SIMULATION_WATER_HEIGHT_RADIAL_DAM_BREAK:
			std::cout << "Water setup: radial dam break" << std::endl;
			break;

		case SIMULATION_WATER_HEIGHT_SINGLE_WAVE_ON_SINGLE_BEACH:
		{
			std::cout << "Water setup: single wave on single beach" << std::endl;
#if CONFIG_ENABLE_NETCDF
			cSingleWaveOnSingleBeach.outputVerboseInformation();
#endif
			break;
		}

		case SIMULATION_WATER_HEIGHT_SINGLE_WAVE_ON_COMPOSITE_BEACH:
		{
			std::cout << "Water setup: single wave on Composite beach" << std::endl;
#if CONFIG_ENABLE_NETCDF
			cSingleWaveOnCompositeBeach.outputVerboseInformation();
#endif
			break;
		}

		default:
			std::cout << "Water setup: [unknown]" << std::endl;
			break;
		}


		switch(cParameters_Datasets.simulation_dataset_0_id)
		{
		case SIMULATION_TERRAIN_HEIGHT_ASAGI:
			// no output since output was already done for water information
			break;

		case SIMULATION_TERRAIN_HEIGHT_SIMPLE_NETCDF:
			// no output since output was already done for water information
			break;

		case SIMULATION_TERRAIN_HEIGHT_DEFAULT:
			std::cout << "Terrain setup: default" << std::endl;
			break;

		case SIMULATION_TERRAIN_HEIGHT_LINEAR_SLOPE:
			std::cout << "Terrain setup: linear slope" << std::endl;
			break;

		case SIMULATION_TERRAIN_HEIGHT_PARABOLA:
			std::cout << "Terrain setup: parabola" << std::endl;
			break;

		case SIMULATION_TERRAIN_HEIGHT_FANCY:
			std::cout << "Terrain setup: fancy" << std::endl;
			break;

		case SIMULATION_TERRAIN_HEIGHT_SINGLE_WAVE_ON_SINGLE_BEACH:
			std::cout << "Terrain setup: single wave on single beach" << std::endl;
			break;

		case SIMULATION_TERRAIN_HEIGHT_SINGLE_WAVE_ON_COMPOSITE_BEACH:
			std::cout << "Terrain setup: single wave on composite beach" << std::endl;
			break;
		}
	}
};



#endif /* CTSUNAMISIMULATION_DATASETS_HPP_ */
