/*
 * Copyright (C) 2011 Technische Universitaet Muenchen
 * This file is part of the Sierpi project. For conditions of distribution and
 * use, please see the copyright notice in the file 'copyright.txt' at the root
 * directory of this package and the copyright notice at http://www5.in.tum.de/sierpi
 *
 *  Created on: Aug 27, 2012
 *      Author: Alexander Breuer <breuera@in.tum.de>,
 *              Martin Schreiber <martin.schreiber@in.tum.de>,
 *              Aditya Ghantasala <shine.aditya@gmail.com>
 */



#ifndef CSINGLEWAVEONCOMPOUNDBEACH_HPP_
#define CSINGLEWAVEONCOMPOUNDBEACH_HPP_

#include <cmath>
#include "netcdfcpp.h"

#include "helper_routines/Variable.hpp"

#ifndef CONFIG_COMPILE_WITHOUT_SIERPI
#	include "../CConfig.hpp"
#	include "../types/CTypes.hpp"
#endif

template <typename T>
class CSingleWaveOnCompositeBeach
{
	typedef CHyperbolicTypes::CSimulationTypes::CNodeData CNodeData;

	/*!
	 * Contains all the functions needed for implementing the
	 * Single Wave on Composite Beach Benchmark problem.
	 * Also this Benchmark problem has three cases depending on the
	 * length(L) of the beach before the slope begins.
	 *
	 * Case A: L =  2.40 meters
	 * Case B: L =  0.98 meters
	 * Case C: L =  0.64 meters
	 *
	 * In this Benchmark the above L is taken as initial mid point of the
	 * wave.
	 *
	 */

private:
	// #################### VAIRABLES #######################
	T d;
	///< Undisturbed water Height. By default it is set to 0.218 meters(Also the maximum depth of water)
	T L;
	///< The length of the beach before the slope of the beach.
	T H;
	///< Total Height of water... h+bathymetry.
	T h;
	///< elevation of water surface above the mean height d (h is read from the .nc file according to the case)
	T gravity;
	///< Gravity in this problem is 9.81 m/s^2;
	T beach_slope_1;
	///< Slope of the first part of the beach.
	T beach_slope_2;
	///< Slope of the second/middle part of the beach.
	T beach_slope_3;
	///< Slope of the last part/part next to the wall of the beach.
	std::string Case;
	///< This variable stores the case which is to be simulated.
	std::string filename;
	///< Name and path of the NetCDF file containing the input values
	long timePosition;
	///< Variable to store the corresponding index of time from NetCDF file
	NcFile* dataFile;
	///< NetCDF File variable
	T* hVect;
	///< Pointer for storing the values of height at Guage G4
	T* tVect;
	///<
	long tVectLength;

	/// benchmark data
	benchMark::CVerify<T> l_h;

	///< Class object for

	// #################### Functions/Methods #######################

	/**
	 * Method to calculate the time index in the NetCDF file
	 * @param i_time 	- Global Simulation Time
	 * @param i_timeTol	- Time tolerance (OPTIONAL --> Set to 0.005 by Default)
	 */
	bool CalculateTimePosition(T i_time, T i_timeTol = T(0.005)){

		bool l_timeStatus = false;

		//searching for the corresponding time index
		for(long i = 0; i < tVectLength; i++) {
			//Calculating the limits for the time
			T timeLlim = (tVect[i]-i_timeTol*tVect[i]);
			T timeHlim = (tVect[i]+i_timeTol*tVect[i]);

			if(timeLlim < i_time && i_time < timeHlim) {
				timePosition = i;
				l_timeStatus = true;
				break;
			}
		}

		return l_timeStatus;
	}




public:


	/**
	 * Constructor for the Single wave on Composite beach Benchmark Class
	 * @param gravity is the acceleration due to gravity. (optional -> set to 9.81)
	 */
	CSingleWaveOnCompositeBeach(const T i_d = T(0.218), const T i_gravity = T(9.81) ){

		setup();

	}
	/**
	 * @param gravity is the acceleration due to gravity. (optional -> set to 9.81)
	 */
	void setup(const T i_d = T(0.218), const T i_gravity = T(9.81)){

		//Variable Initialization
		d = i_d;
		timePosition = 0;
		tVectLength = 0;
		tVect = NULL;
		std::string l_Case = "A";
		//Depending on Case set the L value.
		if(l_Case.compare("A") == 0)
			L = 2.40;
		if(l_Case.compare("B") == 0)
			L = 0.98;
		if(l_Case.compare("C") == 0)
			L = 0.64;
		filename = "data/tsunami_benchmarks/SingleWaveonCompositeBeach_A.nc";
		l_h.setup(filename,"time","x","y","h");
		///< Setup for Verify Class object


		//########## NetCDF Stuff ##########
		//Reading the .nc file in
		dataFile = new NcFile(filename.c_str(),NcFile::ReadOnly);
		if(!dataFile->is_valid())
			std::cout<<"Could not open the file containing the input guage Readings"<<std::endl;

		//Obtaining Time vector length from the NetCDF file
		tVectLength = dataFile->get_dim("time")->size();

		//Nc Variable Declaration
		NcVar* l_G4Var = dataFile->get_var("G4");
		NcVar* l_tVar = dataFile->get_var("time");

		//Memory Allocation
		tVect = new T[dataFile->get_dim("time")->size()];
		hVect = new T[dataFile->get_dim("time")->size()];

		//Reading the vectors of time and G4
		l_G4Var->get(hVect, tVectLength );
		l_tVar->get(tVect, tVectLength );

	}

	/**
	 * Returns the Bathymetry at a given co ordinate
	 *
	 * @param i_x 				- X coordinate
	 * @param i_y 				- Y coordinate
	 * @param i_level_of_detail	- Grid Fineness
	 * @return 					- Returns the Bathymetry at the given X coordinate
	 *
	 */
	T p_getBathymetryData(T i_x, T i_y, T i_level_of_detail){

		T l_bathymetry=0.;
		//If the location is in the first part of the slope
		if(i_x<0 && i_x>-4.36){
			l_bathymetry = -(d+i_x/53);
			return l_bathymetry;
		}
		//If the location is in the second part/middle part of the slope
		if(i_x<0 && i_x<-4.36 && i_x>-7.29){
			l_bathymetry = 4.36/53 - 4.36/150 - i_x/150 -d;
			return l_bathymetry;
		}
		//If the location is in the last part/next to the wall of the slope
		if(i_x<0 && i_x<-7.29 && i_x>=-8.19){
			l_bathymetry = 4.36/53 + 2.93/150 - 7.29/13 -i_x/13 -d;
			return l_bathymetry;
		}
		//If the location is not in the slope region
		if(i_x>=0){
			l_bathymetry = -d;
			return l_bathymetry;
		}
		if(i_x<=-8.19){
			l_bathymetry = 0.1;
		}

		return l_bathymetry;
	}

	/**
	 * Returns the Water height at any given time
	 *
	 * @param i_time 				- Global Simulation Time
	 * @param i_level_of_detail 	- Grid Fineness
	 *
	 */
	T p_getWaterHeightatGhost(T i_time, T i_level_of_detail){

		T l_waterHeight;
		bool timeStatus = CalculateTimePosition(i_time);
		if(timeStatus == true)
			l_waterHeight = hVect[timePosition] - p_getBathymetryData(L,L,i_level_of_detail);//getBathymetry(L) because we are considering only at the left boundary
		else
			l_waterHeight = -p_getBathymetryData(L,L,i_level_of_detail);

		return l_waterHeight;
	}

	/**
	 * Method to get the undisturbed water Height
	 * @return 	- Undisturbed Water Height
	 */
	T getDepth(){
		return d;
	}


	/**
	 * Returns the water height at any given location at initial time.
	 *
	 * @param i_x 					- x Coordinate
	 * @param i_y 					- y Coordinate (in 1D case its Zero)
	 * @param i_level_of_detail 	- Grid Fineness
	 */
	T p_getWaterSufaceData(T i_x, T i_y, T i_level_of_detail){

		if(i_x < -8.19){
			return +0.00;
		}
		else
			return -p_getBathymetryData(i_x,i_y,i_level_of_detail);
	}


	/**
	 * Returns the wave velocity at any given location at initial time.
	 *
	 * @param i_x 					- x coordinate
	 * @param i_y 					- y Coordinate (in 1D case its Zero)
	 * @param i_level_of_detail 	- Grid Fineness
	 */
	T p_getVelocity(T i_x,T i_y, T i_level_of_detail){
		T l_velocity;

		l_velocity = 2*(sqrt(gravity*p_getWaterSufaceData(i_x,i_y,i_level_of_detail)) - sqrt(gravity*d));

		return l_velocity;
	}



	/**
	 * Returns the momentum at any given location at initial time.
	 *
	 * @param i_x					- x coordinate
	 * @param i_y 					- y Coordinate (in 1D case its Zero)
	 * @param i_level_of_detail 	- Grid Fineness
	 *
	 */
	T p_getMomentum(T i_x, T i_y, T i_level_of_detail){

		T l_momentum;
		//CalculateTimePosition(i_time);
		l_momentum = p_getVelocity(i_x,i_y,i_level_of_detail)*hVect[timePosition];

		if(p_getWaterSufaceData(i_x,i_y,i_level_of_detail) <  SIMULATION_TSUNAMI_ZERO_THRESHOLD)
			return (T)0.;

		return l_momentum;
	}


	/**
	 * This function sets the value for the ghost cell according to the INFLOW boundary condition
	 *          @param  	i_ghostVelocity 	- A reference to the Momentum value on the ghost cell
	 * 			@param		time 				- Simulation global time
	 * 			@param		i_ghostVelocity 	- A reference to the velocity on the ghost cell
	 * 			@param		i_ghostHeight  		- A reference to the total Height on the ghost cell
	 *
	 */

	bool getBoundaryData(T i_x, T i_y, T i_level_of_detail, T i_time, CNodeData *o_nodal_data){

		if(CalculateTimePosition(i_time) == true){
			T l_ghostVelocity = 2*(std::sqrt(gravity*p_getWaterHeightatGhost(i_time,i_level_of_detail)) - std::sqrt(gravity*d));
			o_nodal_data->h = p_getWaterHeightatGhost(i_time,i_level_of_detail);
			o_nodal_data->hu = l_ghostVelocity * o_nodal_data->h;
			o_nodal_data->hv = 0; //TODO hard coded ... Check if its ok
			o_nodal_data->b = p_getBathymetryData(i_x,i_y,i_level_of_detail);
			return true;
		}
		else{
			return false;
		}

	}
	/**
	 * Method provides the initial data at given coordinates.
	 * @param i_x				- X Coordinate
	 * @param i_y 				- Y Coordinate
	 * @param i_level_of_detail - Fineness of Grid
	 * @param o_nodal_data 		- Pointer to the object containing all the data about the node.
	 */
	void getNodalData(T i_x,T i_y,T i_level_of_detail, CNodeData *o_nodal_data){

		o_nodal_data->h =  p_getWaterSufaceData(i_x,i_y,i_level_of_detail);
		o_nodal_data->b =  p_getBathymetryData(i_x,i_y,i_level_of_detail);
		o_nodal_data->hu = p_getMomentum(i_x,i_y,i_level_of_detail);
		o_nodal_data->hv = 0; //TODO Hard coded ... Check if its OK.

	}

	/**
	 * Method provides the initial surface elevation data at given coordinates.
	 *
	 * @param i_x				- X Coordinate
	 * @param i_y 				- Y Coordinate
	 * @param i_level_of_detail - Fineness of Grid
	 * @param o_nodal_data 		- Pointer to the object containing all the data about the node.
	 */
	T p_getWaterSurfaceElevation(T i_x, T i_y, T i_level_of_detail){
		return p_getWaterSufaceData(i_x,i_y,i_level_of_detail);
	}


	/**
	 * Method puts out the description of the current benchmark
	 */
	void outputVerboseInformation(){

	}


	/**
	 * return dataset value
	 *
	 * id 0: bathymetry
	 * id 1: displacements
	 */
	CHyperbolicTypes::CSimulationTypes::T getDatasetValue(
			int i_dataset_id,		///< id of dataset to get value from
			T i_x,					///< x-coordinate in model-space
			T i_y,					///< x-coordinate in model-space
			T i_level_of_detail		///< level of detail (0 = coarsest level)
	)
	{
		if (i_dataset_id == 0)
			return p_getBathymetryData(i_x, i_y, i_level_of_detail);

		if (i_dataset_id == 1)
			return p_getWaterSurfaceElevation(i_x, i_y, i_level_of_detail);

		throw(std::runtime_error("invalid dataset id"));
		return 1;
	}


	///< BenchMarking stuff ...


	/**
	 * Method to provide the Benchmark Data at a given location and time
	 *
	 * @param i_x				- X coordinate
	 * @param i_y				- Y coordinate
	 * @param i_level_of_detail	- Grid Fineness
	 * @param i_time			- Simulation Global time
	 * @param o_nodal_data 		- pointer for the nodal data object.
	 */
	void getBenchmarkNodalData(T i_x, T i_y, T i_level_of_detail, T i_time, CNodeData *o_nodal_data){

		T l_dummy=0;
		l_h.getBenchMarkNodalData(i_time,i_x,i_y,i_level_of_detail,&l_dummy);
		//std::cout<<"Calculated Benchmark Value is :: "<<l_dummy<<"   at :: "<< i_x << "  and time :: "<<i_time<<std::endl;
		o_nodal_data->h = l_dummy;

	}

	//Destructor
	~CSingleWaveOnCompositeBeach(){

	}


};

#endif /* CSINGLEWAVEONCOMPOUNDBEACH_HPP_ */
