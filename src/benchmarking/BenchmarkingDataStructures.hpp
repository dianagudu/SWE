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
 * Collection of data structures used for benchmarking and validation
 *
 */

#ifndef BENCHMARKING_DATA_STRUCTURES_HPP_
#define BENCHMARKING_DATA_STRUCTURES_HPP_

#include <string>
#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>

#include "../tools/help.hh"

/**
 * Class that holds data at a specific time
 * in a specific position in the domain
 */
class PointData {
public:
	/**
	 * Constructor
	 *
	 * Sets all the fields in the class to some given values
	 * @param i_t	time of measurement
	 * @param i_x	position of measurement point in x-direction
	 * @param i_y	position of measurement point in y-direction
	 * @param i_h	value of water height at given position and time
	 * @param i_b	value of bathymetry at given position and time
	 * @param i_hu	value of momentum in x-direction at given position and time
	 * @param i_hv	value of momentum in y-direction at given position and time
	 */
	PointData(const float i_t,
			  const float i_x,
			  const float i_y,
			  const float i_h,
			  const float i_b,
			  const float i_hu,
			  const float i_hv):
			  t(i_t),
			  x(i_x), y(i_y),
			  h(i_h), b(i_b),
			  hu(i_hu), hv(i_hv) {};

	/**
	 * overriding the < operator to compare two PointData objects
	 * with respect to time: older objects are "smaller"
	 *
	 * @param i_other	the PointData object the current object is compared against
	 * @return	0 if the time of measurement is the same
	 * 		   -1 if this object has a smaller time (older)
	 * 			1 if this object has a bigger time (more recent)
	 */
	bool operator<(const PointData& i_other) const {
		return (t < i_other.t);
	}

private:

	//! time of measurement
	float t;

	//! measurement position in 2D domain
	float x;
	float y;

	//! simulated values at time t in (x,y)
	float h;
	float b;
	float hu;
	float hv;
};

/**
 *
 */
class TimeSeriesDataReceiver {
public:
	/**
	 * The constructor
	 * Builds an object that collects data over time in
	 * one single point in the 2D domain
	 *
	 * @param i_xPos position of point in x-direction
	 * @param i_yPos	position of point in y-direction
	 */
	TimeSeriesDataReceiver(const std::string i_baseName,
						   const float i_xPos,
				   	   	   const float i_yPos):
				   	   	   xPos(i_xPos),
				   	   	   yPos(i_yPos) {
		// create output file
		// name: i_baseName + "_dynamics_" + xPos + "_" + yPos
		std::ostringstream oss;
		oss << i_baseName << "_dynamics_" << xPos << "_" << yPos;
		outFileName = oss.str();

		// write header in the output file
		ofstream l_outFile(outFileName.c_str());
//		if (l_outFile.is_open())
//			l_outFile << "t\t" << "x=" << xPos << ",yPos=" << yPos << endl;
		l_outFile.close();
	};

	/**
	 * Writes a measurement for a given time
	 *
	 * @param i_time moment in time when the data were measured
	 * @param i_value simulated value for given position and time
	 */
	void writeData(const float i_time, const float i_value) {
		// append value to the output file
		ofstream l_outFile (outFileName.c_str(), ios::app);
		if (l_outFile.is_open()) {
			l_outFile << i_time << " " << i_value << endl;
		}
		l_outFile.close();
	}

private:
	//! fixed position in the 2D domain
	float xPos;
	float yPos;

	//! output file name
	std::string outFileName;

	//! data at different times in fixed position
	// to delete
	std::vector<PointData> timeSeries;
};

/**
 *
 */
class SpatialDataReceiver {
public:
	/**
	 * The constructor
	 * Builds an object that collects data
	 * at a given time over a 2D domain
	 *
	 * @param i_time time of current data collection
	 */
	SpatialDataReceiver(const std::string i_baseName,
					    const float i_time):
					    time(i_time) {
		// create output file
		// name: i_baseName + "_profile_" + time
		std::ostringstream oss;
		oss << i_baseName << "_profile_" << time;
		outFileName = oss.str();
		// write header
		ofstream l_outFile(outFileName.c_str());
//		if (l_outFile.is_open())
//			l_outFile << "x" << "\t" << "y" << "\t" << "t=" << time << endl;
		l_outFile.close();
	};

	/**
	 *
	 */
	void writeData(const float i_offsetX,
				   const float i_offsetY,
				   const float i_dX,
				   const float i_dY,
				   const int i_nX,
				   const int i_nY,
				   const int i_nghosts,
				   const Float2D& i_h,
				   const Float2D& i_b) {
		ofstream l_outFile(outFileName.c_str(), ios::app);
		if (l_outFile.is_open()) {
			for (int i=0; i<i_nX; i++)
				for (int j=0; j<i_nY; j++)
					l_outFile << (i_offsetX + (i+0.5f)*i_dX) << " "
							  << (i_offsetY + (j+0.5f)*i_dY) << " "
							  << (i_h[i+i_nghosts][j+i_nghosts] +
								  i_b[i+i_nghosts][j+i_nghosts])
							  << endl;
		}
		l_outFile.close();
	}

private:
	//! fixed time for snapshot
	float time;

	//! output file name
	std::string outFileName;

	//! data in different positions at fixed time
	// to delete
	std::vector<PointData> spatialData;
};

#endif //BENCHMARKING_DATA_STRUCTURES_HPP_
