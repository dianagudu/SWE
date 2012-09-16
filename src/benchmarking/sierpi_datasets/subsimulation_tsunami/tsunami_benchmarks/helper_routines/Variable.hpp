/*
 * Copyright (C) 2012 Technische Universitaet Muenchen
 * This file is part of the Sierpi project. For conditions of distribution and
 * use, please see the copyright notice in the file 'copyright.txt' at the root
 * directory of this package and the copyright notice at http://www5.in.tum.de/sierpi
 *
 *  Created on: May 12, 2012
 *      Author: Aditya ghantasala <shine.aditya@gmail.com>
 */

#ifndef VARIABLE_HPP_
#define VARIABLE_HPP_

#include <cmath>
#include <cassert>

#include <iostream>
#include <string.h>
#include <netcdfcpp.h>

#include "FindNorm.hpp"
#include "ReadData.hpp"
#include "CheckData.hpp"
#include "ReadNetcdfFile.hpp"



namespace benchMark {
template <typename T>
class CVerify;

}


template <typename T>
class benchMark::CVerify{

private:

	//Stores the norm every time the bench mark is performed
	T* norm;
	int timeStepNumber;
	T* benchMarkTimeVector;

	// Variables from solver
	T xCordinate;
	T yCordinate;


	T benchMarkValue;
	std::string varName;

	// Variables used to store the quantities
	// from previous time a benchmark is performed


	int *counter;

	// Variable of the class which reads the NetCDF file
	benchMark::CReadNetcdfFile<T> ncFile;

	//Pointers to store the values read from the NetCDF file
	T* timeVector;
	T* xVector;
	T* yVector;
	T* heightMatrix;

	benchMark::CReadData<T> readData;



public:

	/**
	 * Empty Constructor for the Class Variable

	Verify(){

	}*/

	/**
	 * Method for setting up various parameters for the object of class Variable
	 *
	 * @param i_fileName			- Name of the NetCDF file containing the Benchmark variable
	 * @param i_timeVarableName		- Name of the Time variable/Dimension given in the NetCDF file.
	 * @param i_xVariableName		- Name of the X variable/Dimension given in the NetCDF file.
	 * @param i_yVariableName		- Name of the Y variable/Dimension given in the NetCDF file.
	 * @param i_heightVariableName	- Name of the Height variable/Dimension given in the NetCDF file.
	 */
	void setup(std::string i_fileName, std::string i_timeVarableName, std::string i_xVariableName, std::string i_yVariableName, std::string i_heightVariableName){

		varName = i_heightVariableName;


		counter = 0;

		//Initializing the positions


		//Reading the NetCDF File given
		ncFile.setup(i_fileName,i_timeVarableName,i_xVariableName,i_yVariableName,i_heightVariableName); //This also reads the file parameters into memory

		//Allocating memory to the Vector Pointers
		timeVector = new T[ncFile.getTimeVectorLength()];
		xVector = new T[ncFile.getXVectorLength()];
		yVector = new  T[ncFile.getYVectorLength()];
		heightMatrix = new T[ncFile.getTimeVectorLength()*ncFile.getXVectorLength()*ncFile.getYVectorLength()];

		ncFile.readFile(timeVector,xVector,yVector,heightMatrix);

		readData.setup(ncFile,timeVector,xVector,yVector,heightMatrix);
		//############ Allocating Memory to Norm and benchMarkTimeVector #############

		benchMarkTimeVector = new T[ncFile.getTimeVectorLength()];
		counter = new int[ncFile.getTimeVectorLength()];

		//Allocating the memory for norm vector
		norm = new T[ncFile.getTimeVectorLength()];

		//Making all the norm variables zero
		for(int i=0;i<ncFile.getTimeVectorLength();i++){
			norm[i] = 0;
		}
		//Making the counter Zero
		for(int i=0;i<ncFile.getTimeVectorLength();i++){
			counter[i] = 1;
		}

	}


	/**
	 * Method to execute the Benchmarking routines and store the norms of the calculated differences in the norm variables.
	 *
	 * @param i_simulationTime	- Global Simulation time of the Simulation.
	 * @param i_xCordinate		- X Coordinate where a benchmark should be carried out.
	 * @param i_yCordinate		- Y Coordinate where a benchmark should be carried out.
	 * @param i_solverVarValue	- Value of height calculated from the Solver which is then compared to the benchmark value.
	 */

	void runBenchmark(const T i_simulationTime, const T i_xCordinate,const T i_yCordinate, const T i_solverVarValue){

		if(ncFile.getYVectorLength() == 1)
			yCordinate = 0;
		else
			yCordinate = i_yCordinate;

		bool l_data_status = readData.getBenchMarkData(i_simulationTime,i_xCordinate,yCordinate,benchMarkValue);

		// ############ If data is present continuing with reading data ############
		if(l_data_status == true){
			long l_timePositon=readData.getTimePosition();
			// ############ Calculating the norm of the difference of the data from solver and .nc file ############
			benchMark::CFindNorm<T> l_findnorm(i_solverVarValue,benchMarkValue,l_timePositon);
			l_findnorm.EucleadianNorm(norm,counter);

		}
	}

	/**
	 * Method to read the Benchmark Value from the NetCDF File
	 *
	 * @param i_simulationTime	- Global Simulation time of the Simulation.
	 * @param i_xCordinate		- X Coordinate where a benchmark should be carried out.
	 * @param i_yCordinate		- Y Coordinate where a benchmark should be carried out.
	 * @param o_benchMarkValue	- Pointer to a variable where the read benchmark Variable should be stored.
	 */
	void getBenchMarkNodalData(const T i_simulationTime, const T i_xCordinate,const T i_yCordinate, const T i_level_of_detail, T* o_benchMarkValue){
		if(i_simulationTime > timeVector[0] && i_simulationTime<5+timeVector[ncFile.getTimeVectorLength()-1]){
			if(ncFile.getYVectorLength() == 1)
				yCordinate = 0;
			else
				yCordinate = i_yCordinate;

			bool data_status = readData.getBenchMarkData(i_simulationTime,i_xCordinate,yCordinate,benchMarkValue);
			if(data_status == true)
				*o_benchMarkValue = benchMarkValue;
		}
		else{
			*o_benchMarkValue = 0;
		}

	}

	/**
	 * Method to complete the Benchmark and output various norms and other parameters
	 * Calculated as a part of the Benchmark.
	 */
	void completeBenchmark(){

		int i;

		//		tools::Logger logger(0, "BENCHMARKING RESULTS");
		int l_length = ncFile.getTimeVectorLength();

		//Checking if the all counter elements are more than zero
		for(i=0;i<l_length;i++){
			assert(counter[i]>0);
		}

		//		logger.printNorms(norm,timeVector,l_length,counter,"Profiles");

	}


	//Destructor
	~CVerify(){

	}


};



#endif /* VARIABLE_HPP_ */
