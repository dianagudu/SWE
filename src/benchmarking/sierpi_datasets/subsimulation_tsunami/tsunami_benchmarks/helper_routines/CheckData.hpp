/*
 * Copyright (C) 2012 Technische Universitaet Muenchen
 * This file is part of the Sierpi project. For conditions of distribution and
 * use, please see the copyright notice in the file 'copyright.txt' at the root
 * directory of this package and the copyright notice at http://www5.in.tum.de/sierpi
 *
 *  Created on: August 22, 2012
 *      Author: Aditya ghantasala <shine.aditya@gmail.com>
 */

#ifndef CHECKDATA_HPP_
#define CHECKDATA_HPP_
#include <stdio.h>
#include <iostream>
#include <string>
#include "netcdfcpp.h"
#include <cmath>
#include "ReadNetcdfFile.hpp"
#define xTolerence 0.2
#define timeTolerence 0.005

namespace benchMark {
template <typename T>
class CCheckData;
}

template <typename T>
class benchMark::CCheckData{

private:

	T time;
	T xCordinate;
	long l_timePosition;
	T* timeVector;
	T* xCordianteVector;
	T* yCordinateVector;
	benchMark::CReadNetcdfFile<T> ncFile;


	/**
	 * Method Calculates the corresponding indices of the given X and Y  coordinates in
	 * NetCDF File
	 * @param i_xPosition	- The X Coordinate
	 * @param i_yPosition	- The Y Coordinate
	 * @return				- Returns True if a corresponding index is found in NetCDF file or else a False
	 */
	bool CalculatePositions(long & i_xPosition, long & i_yPosition, T xTol = T(xTolerence)){

		bool l_position_status = false;
		long xPosit = 0;

		double l_ncX;


		for(long i = 0; i < ncFile.getXVectorLength(); i++) {
			//To calculate the xTol% of the x co ordinate from .nc file
			if(xCordianteVector[i]<0)
				l_ncX = xCordianteVector[i]*-1;
			else
				l_ncX = xCordianteVector[i];

			T l_xHlimit = (xCordianteVector[i]-xTol*l_ncX);
			T l_xLlimit = (xCordianteVector[i]+xTol*l_ncX);

			//Finding the x position
			if( l_xHlimit < xCordinate && xCordinate < l_xLlimit ) {
				xPosit = i;
				//std::cout<<"xPosition found in CheckData"<<xPosit<<std::endl;
				l_position_status = true;
				break;
			}
			l_ncX = 0.0;
		}

		//If it is a 1D file then we set the y position always to zero.
		if(ncFile.getYVectorLength() == 1){
			i_yPosition = 0;
		}
		else{
			i_yPosition = i_xPosition;
		}
		i_xPosition = xPosit;

		return l_position_status;

	}

	/**
	 * Calculates the time index in NetCDF file
	 * @param i_timeTol - Optional input of Time tolerance.
	 * @return 	- TRUE if it finds a corresponding index or else FALSE
	 */
	bool CalculateTimePosition(T i_timeTol = T(timeTolerence)){

		bool l_timeStatus = false;

		//searching for the corresponding time index
		for(long i = 0; i < ncFile.getTimeVectorLength(); i++) {
			//Calculating the limits for the time
			T timeLlim = (timeVector[i]-i_timeTol*timeVector[i]);
			T timeHlim = (timeVector[i]+i_timeTol*timeVector[i]);

			if(timeLlim < time && time < timeHlim) {
				l_timePosition = i;
				l_timeStatus = true;
				break;
			}
		}


		return l_timeStatus;
	}


	/**
	 * Method for Calculation of absolute of Time Difference
	 * @param i_timeDiff
	 */
	void CalculateTimeDifference(T& i_timeDiff){

		//Calculating Time Diff.
		i_timeDiff = time-timeVector[l_timePosition];

		if(i_timeDiff<0)
			i_timeDiff = -1*i_timeDiff;

	}





public:
	/**
	 * Empty constructor for the Class CheckData

	CheckData(){
	}*/
	/**
	 * Method to setup all the parameters for the checkData Class
	 * @param i_ncFile				- Object of class CReadNetcdfFile with the properties of netcdf File
	 * @param i_timeVector			- Pointer for the time Vector read in from NetCDF file.
	 * @param i_xCordinateVector	- Pointer for the X Coordinate Vector read in from NetCDF File.
	 * @param i_yCordinateVector	- Pointer for the Y Coordinate Vector read in from NetCDF File.
	 */
	void setup(benchMark::CReadNetcdfFile<T> i_ncFile, T* i_timeVector, T* i_xCordinateVector, T* i_yCordinateVector){

		ncFile = i_ncFile;
		timeVector = i_timeVector;
		xCordianteVector = i_xCordinateVector;
		yCordinateVector = i_yCordinateVector;

	}


	/**
	 *
	 * @param i_time			- Global Simulation Time.
	 * @param i_xCordinate		- X Coordinate for which a corresponding index in the NetCDF
	 * @param i_xPosition		- A Reference to the X Index variable.
	 * @param i_yPosition		- A Reference to the Y Index Variable.
	 * @param i_timePosition	- A Reference to the time Index Variable.
	 * @param i_timeDiff		- A Reference to the previous Time difference which is used to check the given time.
	 * @return					- Returns a TRUE if all the positions are found or else a FALSE.
	 */
	bool GetPositions(T i_time, T i_xCordinate,T i_yCordinate,long& i_xPosition,long& i_yPosition, long& i_timePosition,T& i_timeDiff){

		time = i_time;
		xCordinate = i_xCordinate;

		T l_timeDiff;
		bool l_positionStatus = false;

		//Calculating time Position
		//Stores the time position in the class variable i_timePosition
		bool l_timeStatus = CalculateTimePosition();

		if(l_timeStatus == false){
			return false;
		}else{

			if(i_timePosition == l_timePosition){

				i_timePosition = l_timePosition;
				CalculateTimeDifference(l_timeDiff);

				if(l_timeDiff<=i_timeDiff){
					l_positionStatus = CalculatePositions(i_xPosition,i_yPosition);

					//If the given x coordinate is in between the coordinates in the Netcdf file and still its given a
					//false position status we do a re calculation with higher tolerance
					if(i_xCordinate > xCordianteVector[0] && i_xCordinate < xCordianteVector[ncFile.getXVectorLength()-1] && l_positionStatus == false)
					{
						l_positionStatus = CalculatePositions(i_xPosition,i_yPosition, 1);
						if(i_xPosition == 0)
							i_xPosition = i_xPosition+92;
					}

					i_timeDiff = l_timeDiff;
				}

				return l_positionStatus;

			}
			else{

				CalculateTimeDifference(l_timeDiff);
				i_timePosition = l_timePosition;
				i_timeDiff = l_timeDiff;
				l_positionStatus = CalculatePositions(i_xPosition,i_yPosition);

				return l_positionStatus;

			}
		}
	}


};

#endif /* CHECKDATA_HPP_ */
