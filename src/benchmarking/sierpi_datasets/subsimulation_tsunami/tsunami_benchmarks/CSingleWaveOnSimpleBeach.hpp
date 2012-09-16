/*
 * Copyright (C) 2011 Technische Universitaet Muenchen
 * This file is part of the Sierpi project. For conditions of distribution and
 * use, please see the copyright notice in the file 'copyright.txt' at the root
 * directory of this package and the copyright notice at http://www5.in.tum.de/sierpi
 *
 *  Created on: Apr 17, 2012
 *
 *      Author: Alexander Breuer <breuera@in.tum.de>,
 *              Martin Schreiber <martin.schreiber@in.tum.de>,
 *              Aditya Ghantasala <shine.aditya@gmail.com>
 */


#ifndef CBENCHMARK_SINGLE_WAVE_ON_SINGLE_BEACH_HPP
#define CBENCHMARK_SINGLE_WAVE_ON_SIMPLE_BEACH_HPP

#include <cmath>
#include <stdexcept>

#include "helper_routines/Variable.hpp"
#include "../datasets_separated/CHyperbolicTypes.hpp"

#ifndef CONFIG_COMPILE_WITHOUT_SIERPI
#	include "../CConfig.hpp"
#	include "../types/CTypes.hpp"
#endif


template <typename T>
class CSingleWaveOnSimpleBeach
{
	typedef CHyperbolicTypes::CSimulationTypes::CNodeData CNodeData;

	/**
	 * see
	 * "Standards, Criteria, and Procedures for NOAA Evaluation of Tsunami Numerical Models"
	 * page 26ff. for more details about this benchmark.
	 *
	 * Sketch:
	 *  _
	 * / \
	 *  |  *
	 * y|   *         ***
	 *  |    *       * | *
	 *  |     *     ** H **
	 *  |      ******  | *************************************
	 *          *                      |
	 *           *                     |
	 *            *                    d
	 *             *                   |
	 *              *                  |
	 *               *****************************************
	 *         |-X_0-|--------------L----------------|
	 *         |-----------------X_s-----------------|
	 *
	 *
	 *         ---------------------------->
	 *                                    x
	 */

	//private:
	//dimensional variables
	/// initial max wave height (surface elevation)
	T H;

	/// water height from the sea floor up to the surface in the constant bathymetry area.
	T d;

	/// slope of the "simple" beach
	T beach_slope;

	/// midpoint of initial wave
	T Xsd;

	/// used gravity constant
	T gravity;

	//non-dimensional variables

	/// ratio: max. surface elevation / water height up to the surface
	T Hd;

	/// wave midpoint
	T Xs;

	/// start position of beach
	T X0;

	T tau;

	T gamma;	// gamma value used for surface wave elevation

	//    T simulation_domain_length;		// length of simulation domain
	//    T dimensionless_scale_factor;	// scale factor for dimensionless space-parameters
	//    T simulation_domain_translate;	// displacement of simulation domain relative in real-world-space

	benchMark::CVerify<T> l_h;

public:
	/**
	 * Constructor of the "single wave on a simple beach"-benchmark
	 */
	CSingleWaveOnSimpleBeach(
			const T i_H = 0.019,	///< maximum height of the initial displacement.
			const T i_d = 1.,		///< water height from the sea floor up to the surface in the constant bathymetry area.
			const T i_gravity = 1.	///< gravity constant.
	) {
		setup(	i_H,
				i_d,
				i_gravity
		);
	}

	/**
	 * Method for setting up the parameters for the class object
	 */
	void setup(
			const T i_H,		///< maximum height of the initial displacement.
			const T i_d,		///< water height from the sea floor up to the surface in the constant bathymetry area.
			const T i_gravity	///< gravity constant.
	)	{
		//set problem specific variables
		//	  origin=70.;
		//	  d = 1.;
		//	  H = 0.019*d;
		//	  gamma = sqrt(.75*H/d);
		//	  L = acosh(sqrt(20.))/gamma;
		//	  xZero = d*19.85;
		//	  xOne = xZero + L;
		//	  g=1.;
		//	  tau=sqrt(d/g);

		// set member variables to input
		H = i_H;
		d = i_d;
		gravity = i_gravity;

		// start of beach
		X0 = 19.85 * d;

		// dimensionless H
		Hd = H/d;

		gamma = std::sqrt(0.25*3.0*Hd);

		// distance: start of the beach - midpoint of the initial displacement
		T l_L = std::log(sqrt(20.)+std::sqrt(pow(sqrt(20.0),2)-1))/gamma; //(C++11)

		// midpoint of wave
		Xs = X0 + l_L;

		//		// midpoint of wave
		//		Xs = 40.0;

		//		// height of wave
		//		H = 0.0185;

		// angle of beach
		// d = 1 (dimensionless)
		beach_slope = d/X0;

		tau = std::sqrt(d/gravity);


		//		simulation_domain_length = i_simulation_domain_length;
		//		dimensionless_scale_factor = i_simulation_domain_length/(Xs*0.8);
		//		simulation_domain_translate = simulation_domain_length*0.8;
		//
		//		Xs *= dimensionless_scale_factor;
		//		X0 *= dimensionless_scale_factor;
		//
		//		H *= dimensionless_scale_factor;


		//Setup of the benchmark Variable
//		l_h.setup("/home/aditya/ncfiles/SingleWaveOnaSimpleBeach_BenchMark_Values_1000.nc","time","x","y","h");
		l_h.setup("/home/diana/workspace_c++/sierpi_datasets/subsimulation_tsunami/tsunami_benchmarks/data/tsunami_benchmarks/SingleWaveOnaSimpleBeach_BenchMark_Values_1000.nc","time","x","y","h");
	


	}


	/**
	 * Get the bathymetry at a specific coordinate.
	 *
	 * \return bathymetry at given coordinates.
	 */
	T p_getBathymetryData(
			T i_x,					///< X coordinate
			T i_y,					///< Y coordinate
			T i_level_of_detail		///< Grid Fineness
	)	{
		T x = i_x;//*simulation_domain_length + simulation_domain_translate;

		//constant bathymetry after X_0
		if (x > X0)
			return -d;//-dimensionless_scale_factor;
#if 0
		std::cout << (X0-x) << std::endl;
		std::cout << (X0-x)*beach_slope << std::endl;
		std::cout << std::endl;
#endif

		//"simple" beach else
		return -d+(X0-x)*beach_slope;//-dimensionless_scale_factor;
	}

	/**
	 * Get the water height (height relative to the sea floor) at a specific coordinate.
	 *
	 * \return 					- Water height.
	 */
	T p_getWaterSufaceData(
			const T i_x,				///< X coordinate
			const T i_y,				///< Y Coordinate
			const T i_level_of_detail	///< Fineness of Grid
	) {
		//        if(cells[i].getCellData().b < (T)0)
		//          cells[i].getCellData().h -= cells[i].getCellData().b;
		T l_waterHeight;

		if( p_getBathymetryData(i_x,i_y,i_level_of_detail) < (T)0 )
			l_waterHeight = p_getWaterSurfaceElevation(i_x,i_y,i_level_of_detail) - p_getBathymetryData(i_x,i_y,i_level_of_detail);
		else
			l_waterHeight =  p_getWaterSurfaceElevation(i_x,i_y,i_level_of_detail);

		return l_waterHeight;
	}


	/**
	 * Get the surface elevation (height relative to the water surface) at a specific coordinate.
	 *
	 * \return
	 */
	T p_getWaterSurfaceElevation(
			T i_x,					///< x-coordinate
			T i_y,					///< y-coordinate
			T i_level_of_detail		///< LOD
	) {
		//      2*cosh(gamma*(x-xOne)/d)/(cosh(2* gamma*(x-xOne)/d)+1);
		//        initialWaveForm *= H*initialWaveForm;

		T x = i_x;//*simulation_domain_length + simulation_domain_translate;

		T surfaceElevation = 2*std::cosh(gamma*(x-Xs)/d)/(std::cosh(2* gamma*(x-Xs)/d)+1);
		surfaceElevation *= Hd * surfaceElevation;

		//      // parameter for sech(x)
		//      T sech_x = gamma * ((x - Xs)/d);///dimensionless_scale_factor);
		//
		//      T sech_tanh_x = CMath::tanh(sech_x);
		//
		//      T sech = 1.0 - sech_tanh_x*sech_tanh_x;
		//
		//      T h = H * sech*sech;
		//
		//      return h;

		return surfaceElevation;
	}

#if 0
	/**
	 * Get the value of tau (gamma/d).
	 * \return tau.
	 */
	T p_getTau() {
		return tau;
	}


	/**
	 * Get the water height from the sea floor up to the surface in the constant bathymetry area.
	 * \return d.
	 */
	T p_getD() {
		return d;
	}


	/**
	 * Get the dimensionless height.
	 *
	 * \return dimensionless height.
	 */
	T p_getHd() {
		return Hd;
	}
#endif

	/**
	 * Get the initial momentum according to the formula for the velocity in the SWOSB-benchmark:
	 *   u(x,0) = - \sqrt{g / d} \eta(x,0)
	 *
	 * \return 					- Returns the initial Momentum at a given X coordinate.
	 */
	T p_getMomentum(
			T i_x,				///< X Coordinate
			T i_y,				///< Y Coordinate
			T i_level_of_detail	///< LOD
	) {
		T l_momentum = p_getWaterSurfaceElevation(i_x,i_y,i_level_of_detail );
		l_momentum *= p_getWaterSufaceData( i_x, i_y, i_level_of_detail );
		l_momentum *= -std::sqrt(gravity / d);


		if(p_getWaterSufaceData(i_x,i_y,i_level_of_detail) < SIMULATION_TSUNAMI_ZERO_THRESHOLD)
			return (T)0.;

		return l_momentum;
	}


	/**
	 * Prints the specifications of the Benchmarks
	 *
	 */

	void outputVerboseInformation(){
		std::cout << std::endl;
		std::cout << "H: " << H << std::endl;
		std::cout << "d: " << d << std::endl;
		std::cout << "H/d: " << Hd << std::endl;
		std::cout << "g: " << gravity << std::endl;
		std::cout << "beach_slope: " << beach_slope << std::endl;
		std::cout << "gamma: " << gamma << std::endl;
		std::cout << "X0 (start of beach slope): " << X0 << std::endl;
		std::cout << "Xs (midpoint of initial wave): " << Xs << std::endl;
		std::cout << "tau: " << tau << std::endl;

		std::cout << std::endl;
	}

	/**
	 * Method assigns Bathymetry, Height, and Momentum for a coordinates
	 *
	 * @param o_nodal_data 		-
	 */
	void getNodalData(
			T i_x,		///< x-coordinate in world-space
			T i_y,		///< x-coordinate in world-space
			T i_level_of_detail,		///< level of detail (0 = coarsest level)
			CNodeData *o_nodal_data	///< pointer for the nodal data object.
	){

		o_nodal_data->b = p_getBathymetryData(i_x,i_y,i_level_of_detail);
		o_nodal_data->h = p_getWaterSufaceData(i_x,i_y,i_level_of_detail);
		o_nodal_data->hu = p_getMomentum(i_x,i_y,i_level_of_detail);
		o_nodal_data->hv = 0;
	}

	/**
	 * Obtains the boundary Cell Values
	 */
	void getBoundaryData(
			T i_x,		///< x-coordinate in world-space
			T i_y,		///< x-coordinate in world-space
			T i_level_of_detail,		///< level of detail (0 = coarsest level)
			CNodeData *o_nodal_data	///< pointer for the nodal data object.
	){
		//Nothing to do here as out let boundary conditions are already implemented in the frame work.
		assert(false);
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

		throw std::runtime_error("invalid dataset id");
		return 1;
	}


	///< BenchMarking stuff ...



	/**
	 * Method to provide the Benchmark Data at a given location and time
	 *
	 * return true to signal that valid benchmark data is available.
	 */
	bool getBenchmarkNodalData(
			T i_x,					///< X coordinate
			T i_y,					///< Y coordinate
			T i_level_of_detail,	///< Grid Fineness
			T i_time,				///< Simulation Global time
			CNodeData *o_nodal_data	///< pointer for the nodal data object.
	){
		T l_dummy=0;
		l_h.getBenchMarkNodalData(i_time, i_x, i_y, i_level_of_detail, &l_dummy);
		o_nodal_data->h = l_dummy;

		return true;
	}



};

#endif
