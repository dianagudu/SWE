/*
 * Copyright (C) 2012 Technische Universitaet Muenchen
 * This file is part of the Sierpi project. For conditions of distribution and
 * use, please see the copyright notice in the file 'copyright.txt' at the root
 * directory of this package and the copyright notice at http://www5.in.tum.de/sierpi
 *
 *  Created on: August 15, 2012
 *      Author: Aditya ghantasala <shine.aditya@gmail.com>
 */

#ifndef READNETCDFFILE_HPP_
#define READNETCDFFILE_HPP_

#include <stdio.h>
#include <string>
#include "netcdfcpp.h"
#include <cmath>
#include <cassert>

namespace benchMark {
template <typename T>
class CReadNetcdfFile;
}

template <typename T>
class benchMark::CReadNetcdfFile{

private:

	long timeVectLength;
	///< length of the time vector in the given NetCDF file.
	long xVectLength;
	///< length of the x coordinate vector in given NetCDF file.
	long yVectLength;
	///< length of the y coordinate vector in given NetCDF file.
	std::string fileName;
	///< Name of the NetCDF File given.
	std::string timeVariableName;
	///< Name of time variable as specified in NetCDF file.
	std::string xVariableName;
	///< Name of the X coordinate as specified in NetCDF file.
	std::string yVariableName;
	///< Name of the Y coordinate as specified in NetCDF file.
	std::string heightVariableName;
	///< Name of height variable as specified in NetCDF file.
	NcFile* dataFile;
	///< NetCDF file variable for the given file.
	NcVar* timeVariable;
	///< NetCDF time variable.
	NcVar* heightVariable;
	///< NetCDF height Variable.
	NcVar* xVariable;
	///< NetCDF X variable;
	NcVar* yVariable;
	///< NetCDF Y variable;



public:
	/**
	 * Empty Constructor for the Class

	ReadNetcdfFile(){
	} */


	/**
	 * Method for setting up the variables and other values for the class
	 *
	 * @param i_fileName  				- Name of the file along with the path which has to be read in to memory.
	 * @param i_timeVariableName		- Name of the time variable/dimensions as mentioned in the NetCDF file.
	 * @param i_xVariableName			- Name of the X variable/dimensions as mentioned in the NetCDF file.
	 * @param i_yVariableName			- Name of the Y variable/dimensions as mentioned in the NetCDF file.
	 * @param i_heightVaribaleName		- Name of the Height variable/dimensions as mentioned in the NetCDF file.
	 */
	void setup(std::string i_fileName, std::string i_timeVariableName, std::string i_xVariableName,
			std::string i_yVariableName, std::string i_heightVaribaleName){


		//Variable Initialization
		fileName = i_fileName;
		timeVariableName = i_timeVariableName;
		xVariableName = i_xVariableName;
		yVariableName = i_yVariableName;
		heightVariableName = i_heightVaribaleName;

		timeVectLength = 0;
		xVectLength = 0;
		yVectLength = 0;

		//Reading the Netcdf File in to the pointer
		dataFile = new NcFile(fileName.c_str(),NcFile::ReadOnly);
		assert(dataFile->is_valid());

		//Obtaining the lengths of all the vectors
		timeVectLength = dataFile->get_dim(timeVariableName.c_str())->size();
		xVectLength = dataFile->get_dim(xVariableName.c_str())->size();
		yVectLength = dataFile->get_dim(yVariableName.c_str())->size();


	}

	/**
	 * Method to read all the values from NetCDF file into the pointers given.
	 *
	 * @param o_timeVector		- Pointer to the Array in to which the time vector in the NetCDF file is to be read.
	 * @param o_xVector			- Pointer to the Array in to which the X vector in the NetCDF file is to be read.
	 * @param o_yVector			- Pointer to the Array in to which the Y vector in the NetCDF file is to be read.
	 * @param o_heightMatrix	- Pointer to the Array in to which the Height Matrix in the NetCDF file is to be read.
	 */
	void readFile(T* o_timeVector, T* o_xVector, T* o_yVector, T* o_heightMatrix){

		//Reading all the variables
		if(!(timeVariable = dataFile->get_var(timeVariableName.c_str())))
			std::cout<<"Variable "<<timeVariableName <<" is not in the given NetCDF File"<<std::endl;
		if(!(xVariable = dataFile->get_var(xVariableName.c_str())))
			std::cout<<"Variable "<<xVariableName <<" is not in the given NetCDF File"<<std::endl;
		//yVariable = dataFile->get_var(yVariableName.c_str());
		if(!(heightVariable = dataFile->get_var(heightVariableName.c_str())))
			std::cout<<"Variable "<<heightVariableName <<" is not in the given NetCDF File"<<std::endl;

		//Reading the vectors into the memory allocated above
		if(!(timeVariable->get(o_timeVector, timeVectLength)))
			std::cout<<"Vector of  "<<timeVariableName <<"is not read from NetCDF File"<<std::endl;
		if(!(xVariable->get(o_xVector,xVectLength)))
			std::cout<<"Vector of  "<<xVariableName<<"is not read from NetCDF File"<<std::endl;
		//yVariable->get(o_yVector,yVectLength);
		if(!(heightVariable->get(o_heightMatrix, yVectLength,timeVectLength,xVectLength)))
			std::cout<<"Vector of  "<<heightVariableName<<"is not read from NetCDF File"<<std::endl;

	}

	/**
	 * Method to obtain the Time vector length in the NetCDF file given.
	 * @return 	- Returns the Length of time Vector in given NetCDF File.
	 */
	long getTimeVectorLength(){
		return timeVectLength;
	}

	/**
	 * Method to obtain the X Coordinate vector length in the NetCDF file given.
	 * @return 	- Returns the length of the X coordinate Vector Length.
	 */
	long getXVectorLength(){
		return xVectLength;
	}

	/**
	 * Method to obtain the Y Coordinate vector length in the NetCDF file given.
	 * @return 	- Returns the length of the Y coordinate Vector Length.
	 */
	long getYVectorLength(){
		return yVectLength;
	}



};


#endif /* READNETCDFFILE_HPP_ */
