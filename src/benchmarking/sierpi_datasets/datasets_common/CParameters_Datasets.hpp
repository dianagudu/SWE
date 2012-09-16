/*
 * CParameters_Datasets.hpp
 *
 *  Created on: Aug 30, 2012
 *      Author: schreibm
 */

#ifndef CPARAMETERS_DATASETS_HPP_
#define CPARAMETERS_DATASETS_HPP_

#include <vector>
#include <string>
#include <limits>

class CParameters_Datasets
{
	typedef CHyperbolicTypes::CSimulationTypes::T	T;

public:

	/**
	 * terrain scene id
	 */
	int simulation_dataset_0_id;

	/**
	 * water surface scene id
	 */
	int simulation_dataset_1_id;

	/**
	 * surface scene id 2
	 */
	int simulation_dataset_2_id;

	/**
	 * surface scene id 3
	 */
	int simulation_dataset_3_id;


	/*******************************************************
	 * DATASET RELATED PARAMETERS
	 */

	/**
	 * parameters for default setup of column
	 */
	T simulation_dataset_breaking_dam_posx;
	T simulation_dataset_breaking_dam_posy;
	T simulation_dataset_breaking_dam_radius;
	T simulation_dataset_breaking_dam_height;

	/**
	 * values to setup the DOFs
	 */

	T simulation_dataset_default_values[4];

	/**
	 * filenames for datasets
	 *
	 * For tsunami simulation:
	 *   0: terrain data
	 *   1: displacement data
	 */
	std::vector<std::string> simulation_datasets;




	/**
	 * domain center (x-coordinate)
	 */
	T simulation_dataset_default_domain_origin_x;

	/**
	 * domain center (y-coordinate)
	 */
	T simulation_dataset_default_domain_origin_y;



	/**
	 * domain size (x-axis)
	 */
	T simulation_dataset_default_domain_size_x;

	/**
	 * domain size (y-axis)
	 */
	T simulation_dataset_default_domain_size_y;


	CParameters_Datasets()	:
		simulation_dataset_0_id(0),
		simulation_dataset_1_id(1),
		simulation_dataset_2_id(0),
		simulation_dataset_3_id(1),


		simulation_dataset_breaking_dam_posx(-1.0),
		simulation_dataset_breaking_dam_posy(-1.0),
		simulation_dataset_breaking_dam_radius(-1.0),
		simulation_dataset_breaking_dam_height(1.0),

		simulation_dataset_default_domain_origin_x(std::numeric_limits<T>::infinity()),
		simulation_dataset_default_domain_origin_y(std::numeric_limits<T>::infinity()),
		simulation_dataset_default_domain_size_x(std::numeric_limits<T>::infinity()),
		simulation_dataset_default_domain_size_y(std::numeric_limits<T>::infinity())
	{
		simulation_dataset_default_values[0] = -std::numeric_limits<T>::infinity();
		simulation_dataset_default_values[1] = -std::numeric_limits<T>::infinity();
		simulation_dataset_default_values[2] = -std::numeric_limits<T>::infinity();
		simulation_dataset_default_values[3] = -std::numeric_limits<T>::infinity();
	}
};



#endif /* CPARAMETERS_DATASETS_HPP_ */
