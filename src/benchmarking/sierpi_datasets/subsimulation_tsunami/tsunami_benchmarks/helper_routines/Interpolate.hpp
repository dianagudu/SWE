/*
 * Copyright (C) 2012 Technische Universitaet Muenchen
 * This file is part of the Sierpi project. For conditions of distribution and
 * use, please see the copyright notice in the file 'copyright.txt' at the root
 * directory of this package and the copyright notice at http://www5.in.tum.de/sierpi
 *
 *  Created on: August 16, 2012
 *      Author: Aditya ghantasala <shine.aditya@gmail.com>
 */

#ifndef INTERPOLATE_HPP_
#define INTERPOLATE_HPP_

#include <stdio.h>
#include <iostream>
#include <string>
#include "netcdfcpp.h"
#include <cmath>
#include "ReadNetcdfFile.hpp"




namespace benchMark {
template <typename T>
class CInterpolate;
}

/**
 * Class Interpolate inherits ReadData to access all the vectors.
 * Class Interpolate contains all the necessary methods to perform an interpolation
 * of the required Benchmark Value from the available values in the given NetCDF Value.
 */

template <typename T>
class benchMark::CInterpolate{



private:


	T* timeVector;
	T* xCordianteVector;
	T* yCordinateVector;
	T* heightMatrix;
	//Variavble for storing the attributes of NetCDF file
	benchMark::CReadNetcdfFile<T> ncFile;


	/**
	 * Method to give the position in a 1D array of Height martix (in row major storage)
	 * 			(This method is for 2D problems)
	 * @param xPosition		- X index calculated.
	 * @param yPosition		- Y index calculated.
	 * @param timePosition	- time index Calculated.
	 */
	long getPosition(long xPosition, long yPosition, long timePosition){


		return 0;
	}

	/**
	 * Method assumes a Row Major saving of a 2D matrix and give the position accordingly
	 *											(For 1D simulation problem)
	 * @param xPosition		- X index calculated.
	 * @param timePosition	- time index Calculated.
	 */
	long getPosition(long i_xPosition,long i_timePosition){

		return i_xPosition*ncFile.getTimeVectorLength() + i_timePosition;

	}



public:

	/**
	 * Empty Constructor for the class Interpolate

	Interpolate(){

	}*/
	/**
	 *
	 *
	 * @param i_ncFile				- The object of class CReadNetcdfFile which read in the NetCDF file.
	 * @param i_timeVector			- Pointer to the time vector which is read in.
	 * @param i_xCordinateVector	- Pointer to the X Coordinate vector which is read in.
	 * @param i_yCordinateVector	- Pointer to the Y Coordinate vector which is read in.
	 * @param i_heightMatrix		- Pointer to the Height Matrix vector which is read in.
	 */
	void setup(benchMark::CReadNetcdfFile<T> i_ncFile,T* i_timeVector, T* i_xCordinateVector, T* i_yCordinateVector, T* i_heightMatrix){

		ncFile = i_ncFile;
		timeVector = i_timeVector;
		xCordianteVector = i_xCordinateVector;
		yCordinateVector = i_yCordinateVector;
		heightMatrix = i_heightMatrix;

	}

	/**
	 * Method Performs Interpolation for the benchmark value from the values in the given NetCDF file
	 *
	 * @param i_simulationTime	- Global simulation time of simulation.
	 * @param i_xCordinate		- X Coordinate where the interpolation is to be done.
	 * @param i_yCordinate		- Y Coordinate where the interpolation is to be done.
	 * @param i_xPosition		- X Index corresponding to above x coordinate.
	 * @param i_yPosition		- Y Index corresponding to above Y coordinate.
	 * @param i_timePosition	- Time Index corresponding to above time.
	 * @return 					Returns the interpolated Value at above x and y coordinates.
	 */
	T DoIntepolation(T i_simulationTime,T i_xCordinate,T i_yCordinate, long i_xPosition,long i_yPosition,long i_timePosition){

		if(i_yCordinate == 0){
			//This is the 1D Case... We do a Bi-Linear Interpolation
			T l_interpolatedValue;

			T l_u00 = heightMatrix[getPosition(i_xPosition,i_timePosition)];
			T l_u01 = heightMatrix[getPosition(i_xPosition+1,i_timePosition)];
			T l_u10 = heightMatrix[getPosition(i_xPosition,i_timePosition+1)];
			T l_u11 = heightMatrix[getPosition(i_xPosition+1,i_timePosition+1)];

			T x = i_xCordinate;
			T y = i_simulationTime;
			T x1= xCordianteVector[i_xPosition];
			T x2= xCordianteVector[i_xPosition+1];
			T y1= yCordinateVector[i_yPosition+1];
			T y2= yCordinateVector[i_yPosition+2];

			//Formula for Bilinear Interpolation

			l_interpolatedValue = (l_u00*(x2-x)*(y2-y))/((x2-x1)*(y2-y1)) +
					(l_u10*(x-x1)*(y2-y))/((x2-x1)*(y2-y1)) +
					(l_u01*(x2-x)*(y-y1))/((x2-x1)*(y2-y1)) +
					(l_u11*(x-x1)*(y-y1))/((x2-x1)*(y2-y1));

			return l_interpolatedValue;




		}
		else{
			//This is the 2D Case ... We do a Tri-linear Interpolation
			return 0.001;

		}




	}



};



#endif /* INTERPOLATE_HPP_ */
