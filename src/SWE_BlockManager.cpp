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

#include <cstdlib>
#include <queue>
#include <vector>

#include "SWE_BlockManager.hh"

#ifdef BENCHMARKING
#include "benchmarking/BenchmarkingDataReceiver.hpp"
#endif

/**
 * Compares two blocks
 */
bool CompareSWE_BlockAMR::operator()(SWE_BlockAMR* b1, SWE_BlockAMR* b2) {
   return (b1->getRefinementLevel() > b2->getRefinementLevel());
}

SWE_BlockManager::SWE_BlockManager(SWE_BlockAMR*** i_blocks,
										   const int i_blockX,
										   const int i_blockY,
										   const InterpolationType i_interpolationScheme):
										   blockX(i_blockX),
										   blockY(i_blockY),
										   interpolationScheme(i_interpolationScheme) {
	for (int i=0; i<blockX; i++)
		for (int j=0; j<blockY; j++)
			blocks.push( i_blocks[i][j]);
	time = 0.;
}

#ifdef BENCHMARKING
void SWE_BlockManager::initBenchmarkingDataReceiver(const std::string i_baseName) {
	receiver = new BenchmarkingDataReceiver(i_baseName);
}

void SWE_BlockManager::addSpatialData(const float i_time) {
	receiver->addSpatialData(i_time);
}

void SWE_BlockManager::addTimeSeries(const float i_xPos,
											 const float i_yPos) {
	receiver->addTimeSeries(i_xPos, i_yPos);
}

#endif

// for global time-stepping, interpolation strategy doesn't matter
float SWE_BlockManager::simulate_gts(const float i_dt_c) {
	float l_t = 0;

	priority_queue<SWE_BlockAMR*, vector<SWE_BlockAMR*>, CompareSWE_BlockAMR> l_pq(blocks);
	vector<SWE_BlockAMR*> l_blocks;
	while (!l_pq.empty()) {
		l_blocks.push_back(l_pq.top());
		l_pq.pop();
	}

	while (l_t < i_dt_c) {
		float l_dt = (float) 600.;

		// update copy layers
		for(vector<SWE_BlockAMR*>::iterator it = l_blocks.begin(); it != l_blocks.end(); ++it) {
			for (int edge=0; edge<4; edge++) {
				if ((*it)->getNeighbour(edge) != NULL)
					(*it)->getNeighbour(edge)->synchCopyLayerBeforeRead(GTS, SWE_BlockAMR::getOppositeEdge(edge), l_t, l_t);
				(*it)->setGhostLayerEdge(edge);
			}
		}

		// set values in ghost cells
//		for(vector<SWE_BlockAMR*>::iterator it = l_blocks.begin(); it != l_blocks.end(); ++it) {
//			(*it)->setGhostLayer();
//		}

		// execute Euler time step:
		for(vector<SWE_BlockAMR*>::iterator it = l_blocks.begin(); it != l_blocks.end(); ++it) {
			(*it)->computeNumericalFluxes();
			float l_stepMax = (*it)->getMaxTimestep();
			l_dt = (l_stepMax < l_dt) ? l_stepMax : l_dt;
		}

		if (l_dt > i_dt_c - l_t)
			l_dt = i_dt_c - l_t;

		// update unknowns
		for(vector<SWE_BlockAMR*>::iterator it = l_blocks.begin(); it != l_blocks.end(); ++it) {
//			(*it)->synchBeforeRead();
			(*it)->updateUnknowns(l_dt);
//			(*it)->synchAfterWrite();
		}

		l_t += l_dt;
		time += l_dt;

#ifdef BENCHMARKING
		for(vector<SWE_BlockAMR*>::iterator it = l_blocks.begin(); it != l_blocks.end(); ++it) {
			receiver->writeData(time,
								l_dt,
								(*it)->getOffsetX(),
								(*it)->getOffsetY(),
								(*it)->getDx(),
								(*it)->getDy(),
								(*it)->getNx(),
								(*it)->getNy(),
								(*it)->getNghosts(),
								(*it)->getWaterHeight(),
								(*it)->getBathymetry());
		}
#endif
	}
	return l_t;
}

float SWE_BlockManager::simulate(const float dt_c) {
	priority_queue<SWE_BlockAMR*, vector<SWE_BlockAMR*>, CompareSWE_BlockAMR> l_pq(blocks);
	return simulate_level(dt_c, l_pq.top()->getRefinementLevel(), l_pq);
}

float SWE_BlockManager::simulate_level(const float i_dt_c,
											   const int i_level,
											   const priority_queue<SWE_BlockAMR*, vector<SWE_BlockAMR*>, CompareSWE_BlockAMR>& i_pq) {
	float l_t = 0;
	// build a temporary vector with the blocks on level l
	vector<SWE_BlockAMR*> l_blocks;
	while (!i_pq.empty() && i_pq.top()->getRefinementLevel() == i_level) {
		i_pq.top()->resetComputationalDomainMax();
		l_blocks.push_back(i_pq.top());
		i_pq.pop();
	}
	int l_num_ts = 0;
	// time-stepping
	while (l_t < i_dt_c) {
		float l_dt = (float) 600.;
		// update copy and ghost layers by using the two-phase update scheme:
		// 1. update left-right boundary cells
		// 2. update bottom-top boundary cells

		for (int edge=0; edge<4; edge++) {
			// update copy layer
			for(vector<SWE_BlockAMR*>::iterator it = l_blocks.begin(); it != l_blocks.end(); ++it) {
				if ((*it)->getNeighbour(edge) != NULL) {
					(*it)->getNeighbour(edge)->synchCopyLayerBeforeRead(LTS, SWE_BlockAMR::getOppositeEdge(edge), l_t, i_dt_c);
				}
			}
			// set ghost layer
			for(vector<SWE_BlockAMR*>::iterator it = l_blocks.begin(); it != l_blocks.end(); ++it) {
				(*it)->setGhostLayerEdge(edge);
			}
		}

		// execute Euler time step:
		for(vector<SWE_BlockAMR*>::iterator it = l_blocks.begin(); it != l_blocks.end(); ++it) {
			(*it)->computeNumericalFluxes();
			float l_stepMax = (*it)->getMaxTimestep();
			l_dt = (l_stepMax < l_dt) ? l_stepMax : l_dt;
		}

		if (l_dt > i_dt_c - l_t)
			l_dt = i_dt_c - l_t;

		// update unknowns
		for(vector<SWE_BlockAMR*>::iterator it = l_blocks.begin(); it != l_blocks.end(); ++it) {
			if (interpolationScheme != SPACE) (*it)->synchBeforeRead();
			(*it)->updateUnknowns(l_dt);
			if (interpolationScheme != SPACE) (*it)->synchAfterWrite();
			// for all blocks that are not the coarsest,
			if (i_level != blocks.top()->getRefinementLevel() && interpolationScheme == SPACE)
				(*it)->decreaseComputationalDomain();
		}

		// simulate blocks on the next refinement level
		if (!i_pq.empty())
			simulate_level(l_dt, i_pq.top()->getRefinementLevel(), i_pq);
#ifdef BENCHMARKING
		else {
			time += l_dt;
			priority_queue<SWE_BlockAMR*, vector<SWE_BlockAMR*>, CompareSWE_BlockAMR> l_tmp(blocks);
			while (!l_tmp.empty()) {
				receiver->writeData(time,
									l_dt,
									l_tmp.top()->getOffsetX(),
									l_tmp.top()->getOffsetY(),
									l_tmp.top()->getDx(),
									l_tmp.top()->getDy(),
									l_tmp.top()->getNx(),
									l_tmp.top()->getNy(),
									l_tmp.top()->getNghosts(),
									l_tmp.top()->getWaterHeight(),
									l_tmp.top()->getBathymetry());
				l_tmp.pop();
			}
		}
#endif

		l_t += l_dt;
		l_num_ts++;
	}
	return l_t;
}
