/**
 * @file
 * This file is part of SWE.
 *
 * @author Diana Gudu
 *
 * @section LICENSE
 *
 * SWE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SWE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SWE.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * @section DESCRIPTION
 *
 * Class for collecting data in a given point in the domain at different times
 *
 */

#ifndef BENCHMARKING_DATA_RECEIVER_HPP_
#define BENCHMARKING_DATA_RECEIVER_HPP_

#include <string>
#include <iostream>
#include <ctime>
#include <cmath>
#include <utility>
#include <map>

#include "BenchmarkingDataStructures.hpp"


/**
 * Class for collecting benchmarking data
 * at different points and times.
 * The data are written to output files on the disk.
 */
class BenchmarkingDataReceiver {
public:
	/**
	 * The constructor
	 * Builds an object that collects benchmarking data
	 *
	 * @param i_basename	base name for output files
	 */
	BenchmarkingDataReceiver(const std::string i_basename): baseName(i_basename) {}

	/**
	 * Adds a given time to the vector of times when
	 * spatial data must be written
	 *
	 * @param i_time time in simulation
	 */
	void addSpatialData(const float i_time) {
		SpatialDataReceiver receiver(baseName, i_time);
		spaceReceivers.insert(pair<float, SpatialDataReceiver> (i_time, receiver));
	}

	/**
	 * Adds a given position in the 2D domain to
	 * the vector of positions where time dependent
	 * data will be written
	 *
	 * @param i_xPos position in the x-direction
	 * @param i_yPos position in the y-direction
	 */
	void addTimeSeries(const float i_xPos,
					   const float i_yPos) {
		TimeSeriesDataReceiver receiver(baseName, i_xPos, i_yPos);
		timeReceivers.insert(pair<pair<float, float>,
								  TimeSeriesDataReceiver> (
								  std::make_pair(i_xPos, i_yPos),
								  receiver)
							 );
	}

	/**
	 * Write data in required positions at a given time. If the
	 * time is close to a required time, write spatial data as well.
	 *
	 * Call this method at each time-step in the simulation for
	 * each SWE_Block to write the water height h
	 *
	 * @param i_time time in simulation
	 * @param i_dt time-step
	 * @param i_offsetX offset in x-direction
	 * @param i_offsetY offset in y-direction
	 * @param i_dX mesh size in x-direction
	 * @param i_dY mesh size in y-direction
	 * @param i_nX number of cells in x-direction
	 * @param i_nY number of cells in y-direction
	 * @param i_nghosts number of ghost cells at each boundary edge
	 * @param i_h 2D array of water height values in each cell, including ghost cells
	 * @param i_b 2D array of bathymetry values in each cell, including ghost cells
	 */
	void writeData(const float i_time,
				   const float i_dt,
				   const float i_offsetX,
				   const float i_offsetY,
				   const float i_dX,
				   const float i_dY,
				   const int i_nX,
				   const int i_nY,
				   const int i_nghosts,
				   const Float2D& i_h,
				   const Float2D& i_b) {
		for ( std::map<pair<float, float>,TimeSeriesDataReceiver>::iterator it = timeReceivers.begin(); it != timeReceivers.end(); ++it ) {
			float l_xPos = it->first.first;
			float l_yPos = it->first.second;
			// condition for required position inside the domain given as input
			if (l_xPos >= i_offsetX && l_xPos < i_offsetX + i_nX*i_dX &&
				l_yPos >= i_offsetY && l_yPos < i_offsetY + i_nY*i_dY) {
				int l_i = (l_xPos - i_offsetX) / i_dX + i_nghosts;
				int l_j = (l_yPos - i_offsetY) / i_dY + i_nghosts;
				it->second.writeData(i_time, i_h[l_i][l_j], i_b[l_i][l_j]);
			}
		}

		for ( std::map<float,SpatialDataReceiver>::iterator it = spaceReceivers.begin(); it != spaceReceivers.end(); ++it ) {
			// if this time is required for benchmarking
			if (abs(i_time - it->first) < i_dt)
				it->second.writeData(i_offsetX, i_offsetY, i_dX, i_dY, i_nX, i_nY, i_nghosts, i_h, i_b);
		}
	}

private:
	//! base name for output files
	std::string baseName;

	//! spatial data receivers for each position
	std::map<float, SpatialDataReceiver> spaceReceivers;

	//! time series receivers for each time
	std::map<pair<float, float>, TimeSeriesDataReceiver> timeReceivers;
};

#endif //BENCHMARKING_DATA_RECEIVER_HPP_
