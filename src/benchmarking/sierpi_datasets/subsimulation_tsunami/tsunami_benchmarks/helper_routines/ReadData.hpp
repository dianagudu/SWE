/*
 * Copyright (C) 2012 Technische Universitaet Muenchen
 * This file is part of the Sierpi project. For conditions of distribution and
 * use, please see the copyright notice in the file 'copyright.txt' at the root
 * directory of this package and the copyright notice at http://www5.in.tum.de/sierpi
 *
 *  Created on: May 6, 2012
 *      Author: Aditya ghantasala <shine.aditya@gmail.com>
 */



#ifndef READDATA_HPP_
#define READDATA_HPP_

#include <stdio.h>
#include <iostream>
#include <string>
#include <cmath>
#include <netcdfcpp.h>

#include "CheckData.hpp"
#include "ReadNetcdfFile.hpp"
#include "Interpolate.hpp"


namespace benchMark {
template <typename T>
class CReadData;
}


template <typename T>
class benchMark::CReadData {

protected:
	T time;
	T xCordinate;
	long l_timePosition;
	T* timeVector;
	T* xCordianteVector;
	T* yCordinateVector;
	T* heightMatrix;

	// Position Variables
	long xPosition;
	long yPosition;
	long timePosition;
	T timeDifference;

	//Variable for Checking the data
	benchMark::CCheckData<T> checkData;

	//Variable for storing the attributes of NetCDF file
	benchMark::CReadNetcdfFile<T> ncFile;

	//Interpolation class object
	benchMark::CInterpolate<T> interpolate;

	/**
	 * Method assumes a Row Major saving of a 2D matrix and give the position accordingly
	 *
	 * @param xPosition		- X index as calculated by other Method.
	 * @param timePosition	- Time Index as calculated by other Method.
	 */
	long getPosition(long i_xPosition,long i_timePosition){

		return i_xPosition + i_timePosition*ncFile.getXVectorLength();

	}




public:
	/**
	 *Empty Constructor for the Class ReadData

	ReadData(){
	}*/
	/**
	 * Method for setting up the variables for the objects of class ReadData
	 *
	 * @param i_ncFile 				- The Object of CReadNetcdfFile which read in the netcdf File.
	 * @param i_timeVector			- Pointer to the Time Vector which is read from NetCDF file.
	 * @param i_xCordinateVector 	- Pointer to the X Coordinate Vector which is read from NetCDF file.
	 * @param i_yCordinateVector	- Pointer to the Y Coordinate Vector which is read from NetCDF file.
	 * @param i_heightMatrix		- Pointer to the Height matrix which is read in by the NetCDF file Reader.
	 */
	void setup(benchMark::CReadNetcdfFile<T> i_ncFile,T* i_timeVector, T* i_xCordinateVector, T* i_yCordinateVector, T* i_heightMatrix){


		ncFile = i_ncFile;
		timeVector = i_timeVector;
		xCordianteVector = i_xCordinateVector;
		yCordinateVector = i_yCordinateVector;
		heightMatrix = i_heightMatrix;
		checkData.setup(i_ncFile, i_timeVector,i_xCordinateVector,i_yCordinateVector);
		interpolate.setup(i_ncFile,i_timeVector,i_xCordinateVector,i_yCordinateVector,i_heightMatrix);

	}


	/**
	 * Method to obtain the Benchmark value at the given time and coordinates
	 *
	 * @param i_simulationTime		- Global simulation time where the benchmark value is required.
	 * @param i_xCordinate			- X coordinate where simulation time is required.
	 * @param i_yCordinate			- Y Coordinate where simulation time is required.
	 * @return  TRUE if it finds data in the NetCDF file or else FALSE
	 */
	bool getBenchMarkData(T i_simulationTime,T i_xCordinate,T i_yCordinate, T &o_benchMarkValue){

		timeDifference = 10;
		bool l_dataStatus = checkData.GetPositions(i_simulationTime,i_xCordinate,i_yCordinate,xPosition,yPosition,timePosition,timeDifference);

		//Making sure that all the positions are positive
		assert(xPosition>=0);
		assert(yPosition>=0);
		assert(timePosition>=0);

		if(l_dataStatus == true){
			//o_benchMarkValue = interpolate.DoIntepolation(i_simulationTime,i_xCordinate, i_yCordinate, xPosition, yPosition, timePosition);
			o_benchMarkValue = heightMatrix[getPosition(xPosition,timePosition)];

			return true;
		}
		else
			return false;

	}

	/**
	 * Method to obtain the time index.
	 */
	long getTimePosition(){
		return timePosition;
	}
	/**
	 * Method to obtain the x coordinate Index
	 */
	long getXPosition(){
		return xPosition;
	}
	/**
	 * Method to obtain the y coordinate Index
	 */
	long getYPosition(){
		return yPosition;
	}

	~CReadData() {
		// TODO Auto-generated destructor stub
	}

};

#endif /* READDATA_HPP_ */
