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
	 * @param i_h simulated water height for given position and time
	 * @param i_b bathymetry for given position and time
	 */
	void writeData(const float i_time, const float i_h, const float i_b) {
		// append value to the output file
		ofstream l_outFile (outFileName.c_str(), ios::app);
		if (l_outFile.is_open())
			l_outFile << i_time << " " << i_h + i_b << endl;
		l_outFile.close();
	}

private:
	//! fixed position in the 2D domain
	float xPos;
	float yPos;

	//! output file name
	std::string outFileName;
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
							  << i_h[i+i_nghosts][j+i_nghosts]
							   + i_b[i+i_nghosts][j+i_nghosts] << endl;

		}
		l_outFile.close();
	}

private:
	//! fixed time for snapshot
	float time;

	//! output file name
	std::string outFileName;
};

#endif //BENCHMARKING_DATA_STRUCTURES_HPP_
