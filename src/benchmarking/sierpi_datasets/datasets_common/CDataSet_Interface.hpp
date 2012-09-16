/*
 * CDataSet.hpp
 *
 *  Created on: Jul 3, 2012
 *      Author: Martin Schreiber <martin.schreiber@in.tum.de>
 */

#ifndef CDATASET_INTERFACE_HPP_
#define CDATASET_INTERFACE_HPP_

#include <vector>
#include <string>

// only include when compiled with sierpi
#ifndef CONFIG_COMPILE_WITHOUT_SIERPI
#	include "../subsimulation_generic/CConfig.hpp"
#	include "../subsimulation_generic/types/CTypes.hpp"
#else
#	include "CHyperbolicTypes.hpp"
#endif

/**
 * The following interfaces have to be provided by the different data sets
 */
class CDataSet_Interface
{
public:
	typedef CONFIG_DEFAULT_FLOATING_POINT_TYPE T;
	typedef CHyperbolicTypes::CSimulationTypes::CNodeData CNodeData;

	virtual ~CDataSet_Interface()	{}


	/**
	 * load the datasets specified via cSimulationParameters
	 */
	virtual bool loadDatasets(
			const std::vector<std::string> &i_datasets
	) = 0;

	/**
	 * return true if the dataset was successfully loaded
	 */
	virtual bool isDatasetLoaded() = 0;


	/**
	 * load terrain origin coordinate and size
	 */
	virtual void getOriginAndSize(
			T	*o_origin_x,	///< origin of domain in world-space
			T	*o_origin_y,	///< origin of domain in world-space
			T	*o_size_x,		///< size of domain in world-space
			T	*o_size_y		///< size of domain in world-space
	) = 0;


	/**
	 * get nodal data (surface height, momentums and bathymetry) for given coordinate for setup (t=0)
	 */
	virtual void getNodalData(
			T i_x,		///< x-coordinate in world-space
			T i_y,		///< x-coordinate in world-space
			T i_level_of_detail,		///< level of detail (0 = coarsest level)
			CNodeData *o_nodal_data		///< nodal data is written to this position
	) = 0;


	/**
	 * get boundary data at given position
	 */
	virtual void getBoundaryData(
			T i_x,		///< x-coordinate in world-space
			T i_y,		///< y-coordinate in world-space
			T i_level_of_detail,	///< level of detail (0 = coarsest level)
			T i_timestamp,			///< timestamp for boundary data
			CNodeData *o_nodal_data	///< nodal data is written to this position
	) = 0;



	/**
	 * get benchmark (analytical) nodal data
	 */
	virtual bool getBenchmarkNodalData(
			T i_x,		///< x-coordinate in world-space
			T i_y,		///< y-coordinate in world-space
			T i_level_of_detail,	///< level of detail (0 = coarsest level)
			T i_timestamp,			///< timestamp for boundary data
			CNodeData *o_nodal_data	///< nodal data is written to this position
	) = 0;



	/**
	 * get dataset data for dataset i_dataset_id
	 *
	 * \return dataset value
	 */
	virtual T getDatasetValue(
			int i_dataset_id,	/// id of dataset
			T i_x,				///< x-coordinate in world-space
			T i_y,				///< x-coordinate in world-space
			T i_level_of_detail	///< level of detail (0 = coarsest level)
	) = 0;



	/**
	 * output some verbose information about this dataset
	 */
	virtual void outputVerboseInformation() = 0;
};


#endif /* CDATASET_HPP_ */
