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
 * TODO
 */

#ifndef SWE_BLOCKMANAGER_HH_
#define SWE_BLOCKMANAGER_HH_

#include <queue>
#include "SWE_BlockAMR.hh"

#ifdef BENCHMARKING
#include "benchmarking/BenchmarkingDataReceiver.hpp"
#endif

class CompareSWE_BlockAMR {
public:
    bool operator()(SWE_BlockAMR* b1, SWE_BlockAMR* b2);
};

class SWE_BlockManager {
public:
	SWE_BlockManager(SWE_BlockAMR*** i_blocks,
						 const int i_blockX=1,
						 const int i_blockY=1,
						 const InterpolationType i_interpolationScheme = APPROX_TIME_SPACE);
	float simulate_gts(const float i_dt_c);
	float simulate(const float i_dt);
	float simulate_level(const float i_dt_c,
						 const int i_level,
						 priority_queue<SWE_BlockAMR*, vector<SWE_BlockAMR*>, CompareSWE_BlockAMR> i_pq);
	void initBenchmarkingDataReceiver(const std::string i_baseName);

#ifdef BENCHMARKING
	void addSpatialData(const float i_time);
	void addTimeSeries(const float i_xPos, const float i_yPos);
#endif

protected:
	//! number of blocks in each dimension
	int blockX;
	int blockY;
	//! interpolation scheme used for ghost layers
	InterpolationType interpolationScheme;
	//! blocks sorted by refinement levels
	priority_queue<SWE_BlockAMR*, vector<SWE_BlockAMR*>, CompareSWE_BlockAMR> blocks;
	//! global time
	float time;
	//! benchmarking data receiver
#ifdef BENCHMARKING
	BenchmarkingDataReceiver* receiver;
#endif
};

#endif /* SWE_BLOCKMANAGER_HH_ */
