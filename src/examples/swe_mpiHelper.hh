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
 * helper functions for ghost cell exchange with MPI
 */

#ifndef __SWE_MPI_HELPER_HH
#define __SWE_MPI_HELPER_HH

#include <mpi.h>
#include "../SWE_BlockGhost.hh"
#include "../SWE_WavePropagationAMR.hh"

float simulateLTSTimeSpace(SWE_WavePropagationAMR* i_wavePropagationBlock,
				  const float i_t, const float i_dtCoarse,

				  const int i_leftNeighborRank,
				  SWE_BlockGhost* o_leftInflow,
				  SWE_BlockGhost* i_leftOutflow,
				  const int i_leftNeighbourRefinementLevel,

				  const int i_rightNeighborRank,
				  SWE_BlockGhost* o_rightInflow,
				  SWE_BlockGhost* i_rightOutflow,
				  const int i_rightNeighbourRefinementLevel,

				  const int i_bottomNeighborRank,
				  SWE_BlockGhost* o_bottomInflow,
				  SWE_BlockGhost* i_bottomOutflow,
				  const int i_bottomNeighbourRefinementLevel,

				  const int i_topNeighborRank,
				  SWE_BlockGhost* o_topInflow,
				  SWE_BlockGhost* i_topOutflow,
				  const int i_topNeighbourRefinementLevel,

				  const int i_refinementLevel,
				  const MPI_Datatype i_mpiCols,
				  const MPI_Datatype i_mpiRows,
				  MPI_Datatype i_mpiSendType[4],
				  MPI_Comm i_refinementLevelComm);

// Exchanges the left and right ghost layers.
void exchangeLeftRightGhostLayers(SWE_WavePropagationAMR* i_wavePropagationBlock,
								  const float i_t, const float i_dtCoarse,
								  const TimeSteppingType i_timeSteppingStrategy,
								  const int i_leftNeighborRank,
						  		  SWE_BlockGhost* o_leftInflow,
								  SWE_BlockGhost* i_leftOutflow,
								  const int i_leftNeighbourRefinementLevel,
								  const int i_rightNeighborRank,
								  SWE_BlockGhost* o_rightInflow,
								  SWE_BlockGhost* i_rightOutflow,
								  const int i_rightNeighbourRefinementLevel,
								  const int i_refinementLevel,
								  const MPI_Datatype i_mpiCols,
								  MPI_Datatype i_mpiSendType[4]);

// Exchanges the bottom and top ghost layers.
void exchangeBottomTopGhostLayers(SWE_WavePropagationAMR* i_wavePropagationBlock,
								  const float i_t, const float i_dtCoarse,
								  const TimeSteppingType i_timeSteppingStrategy,
								  const int i_bottomNeighborRank,
								  SWE_BlockGhost* o_bottomNeighborInflow,
								  SWE_BlockGhost* i_bottomNeighborOutflow,
								  const int i_bottomNeighbourRefinementLevel,
								  const int i_topNeighborRank,
								  SWE_BlockGhost* o_topNeighborInflow,
								  SWE_BlockGhost* i_topNeighborOutflow,
								  const int i_topNeighbourRefinementLevel,
								  const int i_refinementLevel,
								  const MPI_Datatype i_mpiRows,
								  MPI_Datatype i_mpiSendType[4]);

// send ghost layers
void sendGhostLayers(SWE_WavePropagationAMR* i_wavePropagationBlock,
					  const float i_t,
					  const float i_dtCoarse,
					  const TimeSteppingType i_timeSteppingStrategy,
					  const int i_neighborRank,
					  SWE_BlockGhost* i_outflow,
					  MPI_Datatype i_mpiSendType,
					  BoundaryEdge i_edge);

// send ghost layers using Isend
void isendGhostLayers(SWE_WavePropagationAMR* i_wavePropagationBlock,
					  const float i_t,
					  const float i_dtCoarse,
					  const TimeSteppingType i_timeSteppingStrategy,
					  const int i_neighborRank,
					  SWE_BlockGhost* i_outflow,
					  MPI_Datatype i_mpiSendType,
					  BoundaryEdge i_edge,
					  MPI_Request* i_request);

// receive ghost layers
void receiveGhostLayers(SWE_WavePropagationAMR* i_wavePropagationBlock,
						  const int i_neighborRank,
						  SWE_BlockGhost* o_inflow,
						  const MPI_Datatype i_mpiRecvType);


#endif // __SWE_MPI_HELPER_HH
