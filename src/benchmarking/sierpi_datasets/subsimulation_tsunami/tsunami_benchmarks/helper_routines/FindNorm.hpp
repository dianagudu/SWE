/*
 * Copyright (C) 2012 Technische Universitaet Muenchen
 * This file is part of the Sierpi project. For conditions of distribution and
 * use, please see the copyright notice in the file 'copyright.txt' at the root
 * directory of this package and the copyright notice at http://www5.in.tum.de/sierpi
 *
 *  Created on: May 12, 2012
 *      Author: Aditya ghantasala <shine.aditya@gmail.com>
 */
#ifndef FINDNORM_HPP_
#define FINDNORM_HPP_

#include <stdio.h>
#include <iostream>
#include <string>
#include <netcdfcpp.h>
#include <cmath>

namespace benchMark {
  template <typename T>
  class CFindNorm;
}

template <typename T>
class benchMark::CFindNorm {

private:
	T l_solverValue;
	T l_benchMarkValue;
	int l_timeStepNumber;


public:

	//Consturctor
	/**
	 * Constructor for the class.
	 *
	 * @param i_solverValue		- Value calculated by the Solver which is to be compared with the benchmark Value
	 * @param i_timeStepNumber		- The time Index corresponding to the present simulation time.
	 */
	CFindNorm(const T i_solverValue, const T i_benchMarkValue,const int i_timeStepNumber){

		l_solverValue = i_solverValue;
		l_benchMarkValue = i_benchMarkValue;
		l_timeStepNumber = i_timeStepNumber;
	}

	//######################################
	//########   Various Norms   ###########
	//######################################

	/**
	 * Method to calculate the Eucleadean Norm.
	 *
	 * @param i_norm 	- Pointer to the norm Vector where the calculated norm is stored(corresponding to the time index).
	 * @param i_counter	- Pointer to the counter vector which keeps track of how many times a norm is calculated to normalize the calculated norm.
	 */
	void EucleadianNorm(T* i_norm, int* i_counter){

		i_norm[l_timeStepNumber] = i_norm[l_timeStepNumber] + std::pow((l_solverValue-l_benchMarkValue),2);
		i_counter[l_timeStepNumber] = i_counter[l_timeStepNumber]+1;

	}

	void InfNorm(T& i_norm){


	}
	//Distuctor
	~CFindNorm(){

	}

};

#endif /* NORM_HPP_ */
