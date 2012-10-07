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

#include "swe_mpiHelper.hh"

/**
 * Simulation using local time-stepping and (approximate-)time-space interpolation
 * of the ghost cells
 */
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
				  MPI_Comm i_refinementLevelComm) {

	// fixed time-step for coarsest grid
	float l_dt = 0.12;

	// allocate space for send requests
	int l_numRequestsLeft = (i_leftNeighbourRefinementLevel > i_refinementLevel ? i_leftNeighbourRefinementLevel / i_refinementLevel : 0);
	int l_numRequestsRight = (i_rightNeighbourRefinementLevel > i_refinementLevel ? i_rightNeighbourRefinementLevel / i_refinementLevel : 0);
	int l_numRequestsBottom = (i_bottomNeighbourRefinementLevel > i_refinementLevel ? i_bottomNeighbourRefinementLevel / i_refinementLevel : 0);
	int l_numRequestsTop = (i_topNeighbourRefinementLevel > i_refinementLevel ? i_topNeighbourRefinementLevel / i_refinementLevel : 0);

	for (int l_ts = 0; l_ts < i_refinementLevel; l_ts++) {
		// coarse grids - receive ghost layers
		// fine grids - send ghost layers if needed
		// exchange ghost layers for same resolution grids
		if (i_leftNeighbourRefinementLevel >= i_refinementLevel)
			receiveGhostLayers(i_wavePropagationBlock, i_leftNeighborRank, o_leftInflow, i_mpiCols);
		if (i_rightNeighbourRefinementLevel <= i_refinementLevel && i_rightNeighbourRefinementLevel != 0)
			if (i_rightNeighbourRefinementLevel == i_refinementLevel || l_ts % ( i_refinementLevel / i_rightNeighbourRefinementLevel ) == 0)
				sendGhostLayers(i_wavePropagationBlock, 0, 0, LTS, i_rightNeighborRank, i_rightOutflow, i_mpiSendType[BND_RIGHT], BND_RIGHT);
		if (i_rightNeighbourRefinementLevel >= i_refinementLevel)
			receiveGhostLayers(i_wavePropagationBlock, i_rightNeighborRank, o_rightInflow, i_mpiCols);
		if (i_leftNeighbourRefinementLevel <= i_refinementLevel && i_leftNeighbourRefinementLevel != 0)
			if (i_leftNeighbourRefinementLevel == i_refinementLevel || l_ts % ( i_refinementLevel / i_leftNeighbourRefinementLevel ) == 0)
				sendGhostLayers(i_wavePropagationBlock, 0, 0, LTS, i_leftNeighborRank, i_leftOutflow, i_mpiSendType[BND_LEFT], BND_LEFT);


		if (i_bottomNeighbourRefinementLevel >= i_refinementLevel)
			receiveGhostLayers(i_wavePropagationBlock, i_bottomNeighborRank, o_bottomInflow, i_mpiRows);
		if (i_topNeighbourRefinementLevel <= i_refinementLevel && i_topNeighbourRefinementLevel != 0)
			if (i_topNeighbourRefinementLevel == i_refinementLevel || l_ts % ( i_refinementLevel / i_topNeighbourRefinementLevel ) == 0)
				sendGhostLayers(i_wavePropagationBlock, 0, 0, LTS, i_topNeighborRank, i_topOutflow, i_mpiSendType[BND_TOP], BND_TOP);
		if (i_topNeighbourRefinementLevel >= i_refinementLevel)
			receiveGhostLayers(i_wavePropagationBlock, i_topNeighborRank, o_topInflow, i_mpiRows);
		if (i_bottomNeighbourRefinementLevel <= i_refinementLevel && i_bottomNeighbourRefinementLevel != 0)
			if (i_bottomNeighbourRefinementLevel == i_refinementLevel || l_ts % ( i_refinementLevel / i_bottomNeighbourRefinementLevel ) == 0)
				sendGhostLayers(i_wavePropagationBlock, 0, 0, LTS, i_bottomNeighborRank, i_bottomOutflow, i_mpiSendType[BND_BOTTOM], BND_BOTTOM);

		// receive ghost layers from coarser grids, if any
		// left neighbour
		if (i_leftNeighbourRefinementLevel < i_refinementLevel)
			receiveGhostLayers(i_wavePropagationBlock, i_leftNeighborRank, o_leftInflow, i_mpiCols);
		// right neighbour
		if (i_rightNeighbourRefinementLevel < i_refinementLevel)
			receiveGhostLayers(i_wavePropagationBlock, i_rightNeighborRank, o_rightInflow, i_mpiCols);
		// bottom neighbour
		if (i_bottomNeighbourRefinementLevel < i_refinementLevel)
			receiveGhostLayers(i_wavePropagationBlock, i_bottomNeighborRank, o_bottomInflow, i_mpiRows);
		// top neighbour
		if (i_topNeighbourRefinementLevel < i_refinementLevel)
			receiveGhostLayers(i_wavePropagationBlock, i_topNeighborRank, o_topInflow, i_mpiRows);


		// set values in ghost cells
		i_wavePropagationBlock->setGhostLayer();

		// compute numerical flux on each edge
		i_wavePropagationBlock->computeNumericalFluxes();

		// keep the copy layer at the beginning of time-step before updating the unknowns
		i_wavePropagationBlock->synchBeforeRead();

		// update the cell values
		i_wavePropagationBlock->updateUnknowns(l_dt/i_refinementLevel);

		// keep the copy layer at the end of the time-step for time-interpolation
		i_wavePropagationBlock->synchAfterWrite();

		// send all ghost layers to finer neighbours (neighbour refinement level number of sends)
		// left neighbour
		if (i_leftNeighbourRefinementLevel > i_refinementLevel)
			for (int l_fineTs = 0; l_fineTs < l_numRequestsLeft; l_fineTs++)
				sendGhostLayers(i_wavePropagationBlock,
									l_dt * l_fineTs * i_refinementLevel / i_leftNeighbourRefinementLevel,
									l_dt,
									LTS,
									i_leftNeighborRank,
									i_leftOutflow,
									i_mpiSendType[BND_LEFT],
									BND_LEFT);
		// right neighbour
		if (i_rightNeighbourRefinementLevel > i_refinementLevel)
			for (int l_fineTs = 0; l_fineTs < l_numRequestsRight; l_fineTs++)
				sendGhostLayers(i_wavePropagationBlock,
								 l_dt * l_fineTs * i_refinementLevel / i_rightNeighbourRefinementLevel,
								 l_dt,
								 LTS,
								 i_rightNeighborRank,
								 i_rightOutflow,
								 i_mpiSendType[BND_RIGHT],
								 BND_RIGHT);
		// bottom neighbour
		if (i_bottomNeighbourRefinementLevel > i_refinementLevel)
			for (int l_fineTs = 0; l_fineTs < l_numRequestsBottom; l_fineTs++)
				sendGhostLayers(i_wavePropagationBlock,
								  l_dt * l_fineTs * i_refinementLevel / i_bottomNeighbourRefinementLevel,
								  l_dt,
								  LTS,
								  i_bottomNeighborRank,
								  i_bottomOutflow,
								  i_mpiSendType[BND_BOTTOM],
								  BND_BOTTOM);
		// top neighbour
		if (i_topNeighbourRefinementLevel > i_refinementLevel)
			for (int l_fineTs = 0; l_fineTs < l_numRequestsTop; l_fineTs++)
				sendGhostLayers(i_wavePropagationBlock,
								   l_dt * l_fineTs * i_refinementLevel / i_topNeighbourRefinementLevel,
								   l_dt,
								   LTS,
								   i_topNeighborRank,
								   i_topOutflow,
								   i_mpiSendType[BND_TOP],
								   BND_TOP);
	}
	return l_dt;
}

/**
 * TODO
 * LTS with (approx-)time-space interpolation:
 * - use non-blocking communication when sending fine ghost layers to fine grids
 *
 * Check speedup:
 * 2x2 blocks of 100x100, one of them refined by 2
 * LTS time-space interpolation
 * serial: 49s
 * parallel 4procs (2cores): 53s
 *
 * 2x1 blocks of 100x100, right one refined by 2
 *
 * LTS approximate time-space interpolation
 * serial: 27s
 * parallel 2procs: 20s -> 1.35
 *
 * LTS time-space interpolation
 * serial: 40s
 * parallel 2procs: 22s -> 1.81
 *
 * LTS space interpolation
 * serial: 40s
 * parallel 2procs with allreduce: 43s -> 0.93
 * parallel 2procs with fixed ts:  19s -> 2.1
 *
 * 2x1 blocks of 200x200, right one refined by 2i_refinementLevel

 * LTS space interpolation
 * serial: 356s
 * parallel 2procs with allreduce: 352s -> 1.01
 * parallel 2procs with fixed ts: 86s -> 4.13
 */

/**
 * The updated ghost layers are needed to refine a copy layer
 * Therefore, the coarser grids first have to receive the ghost layers,
 * update (refine) the copy layers and then send them to the finer grids.
 * The finer grids first update (coarsen) the copy layers, send them to the
 * coarser grids and then receive the ghost layers.
 * This scheme is implemented with blocking operations: MPI_Send, MPI_Recv
 */

/**
 * Exchanges the left and right ghost layers with MPI's SendReceive.
 *
 * @param i_wavePropagationBlock pointer to the SWE block on this MPI process.
 * @param i_t time when the exchange is made.
 * @param i_leftNeighborRank MPI rank of the  left neighbor.
 * @param o_leftInflow ghost layer, where the left neighbor writes into.
 * @param i_leftOutflow layer where the left neighbor reads from.
 * @param i_leftNeighbourRefinementLevel refinement level of the left neighbour
 * @param i_rightNeighborRank MPI rank of the right neighbor.
 * @param o_rightInflow ghost layer, where the right neighbor writes into.
 * @param i_rightOutflow layer, where the right neighbor reads form.
 * @param i_rightNeighbourRefinementLevel refinement level of the right neighbour
 * @param i_mpiCols MPI data type for the vertical ghost layers.
 * @param i_mpiSendType data types for the all the copy layers of the current block.
 */
void exchangeLeftRightGhostLayers( SWE_WavePropagationAMR* i_wavePropagationBlock, const float i_t, const float i_dtCoarse, const TimeSteppingType i_timeSteppingStrategy,
								   const int i_leftNeighborRank,  SWE_BlockGhost* o_leftInflow,  SWE_BlockGhost* i_leftOutflow, const int i_leftNeighbourRefinementLevel,
                                   const int i_rightNeighborRank, SWE_BlockGhost* o_rightInflow, SWE_BlockGhost* i_rightOutflow, const int i_rightNeighbourRefinementLevel,
                                   const int i_refinementLevel, const MPI_Datatype i_mpiCols,  MPI_Datatype i_mpiSendType[4]) {

	  // if left neighbour is coarser, send and then receive
	  if (i_leftNeighbourRefinementLevel < i_refinementLevel) {
		  sendGhostLayers(i_wavePropagationBlock, i_t, i_dtCoarse, i_timeSteppingStrategy, i_leftNeighborRank, i_leftOutflow, i_mpiSendType[BND_LEFT], BND_LEFT);
		  receiveGhostLayers(i_wavePropagationBlock, i_leftNeighborRank, o_leftInflow, i_mpiCols);
	  }

	  // if right neighbour is finer, receive and then send
	  if (i_rightNeighbourRefinementLevel > i_refinementLevel) {
		  receiveGhostLayers(i_wavePropagationBlock, i_rightNeighborRank, o_rightInflow, i_mpiCols);
		  sendGhostLayers(i_wavePropagationBlock, i_t, i_dtCoarse, i_timeSteppingStrategy, i_rightNeighborRank, i_rightOutflow, i_mpiSendType[BND_RIGHT], BND_RIGHT);
	  }

	  // if right neighbour is coarser, send and then receive
	  if (i_rightNeighbourRefinementLevel <= i_refinementLevel) {
		  sendGhostLayers(i_wavePropagationBlock, i_t, i_dtCoarse, i_timeSteppingStrategy, i_rightNeighborRank, i_rightOutflow, i_mpiSendType[BND_RIGHT], BND_RIGHT);
		  receiveGhostLayers(i_wavePropagationBlock, i_rightNeighborRank, o_rightInflow, i_mpiCols);
	  }

	  // if left neighbour is finer, receive and then send
	  if (i_leftNeighbourRefinementLevel >= i_refinementLevel) {
		  receiveGhostLayers(i_wavePropagationBlock, i_leftNeighborRank, o_leftInflow, i_mpiCols);
		  sendGhostLayers(i_wavePropagationBlock, i_t, i_dtCoarse, i_timeSteppingStrategy, i_leftNeighborRank, i_leftOutflow, i_mpiSendType[BND_LEFT], BND_LEFT);
	  }
}

/**
 * Exchanges the bottom and top ghost layers with MPI blocking routines Send and Recv.
 *
 * @param i_wavePropagationBlock pointer to the SWE block on this MPI process.
 * @param i_t time when the exchange is made.
 * @param i_bottomNeighborRank MPI rank of the bottom neighbor.
 * @param o_bottomInflow ghost layer, where the bottom neighbor writes into.
 * @param i_bottomOutflow host layer, where the bottom neighbor reads from.
 * @param i_bottomNeighbourRefinementLevel refinement level for the bottom neighbour
 * @param i_topNeighborRank MPI rank of the top neighbor.
 * @param o_topInflow ghost layer, where the top neighbor writes into.
 * @param i_topOutflow ghost layer, where the top neighbor reads from.
 * @param i_topNeighbourRefinementLevel refinement level for the top neighbour
 * @param i_mpiRow MPI data type for the horizontal ghost layers.
 * @param i_mpiSendType data types for the all the copy layers of the current block.
 */
void exchangeBottomTopGhostLayers( SWE_WavePropagationAMR* i_wavePropagationBlock, const float i_t, const float i_dtCoarse, const TimeSteppingType i_timeSteppingStrategy,
		  	  	  	  	  	  	   const int i_bottomNeighborRank, SWE_BlockGhost* o_bottomInflow, SWE_BlockGhost* i_bottomOutflow, const int i_bottomNeighbourRefinementLevel,
                                   const int i_topNeighborRank,    SWE_BlockGhost* o_topInflow,    SWE_BlockGhost* i_topOutflow, const int i_topNeighbourRefinementLevel,
                                   const int i_refinementLevel, const MPI_Datatype i_mpiRows,   MPI_Datatype i_mpiSendType[4]) {
	  // if bottom neighbour is coarser, send and then receive
	  if (i_bottomNeighbourRefinementLevel <= i_refinementLevel) {
		  sendGhostLayers(i_wavePropagationBlock, i_t, i_dtCoarse, i_timeSteppingStrategy, i_bottomNeighborRank, i_bottomOutflow, i_mpiSendType[BND_BOTTOM], BND_BOTTOM);
		  receiveGhostLayers(i_wavePropagationBlock, i_bottomNeighborRank, o_bottomInflow, i_mpiRows);
	  }

	  // if top neighbour is finer, receive and then send
	  if (i_topNeighbourRefinementLevel > i_refinementLevel) {
		  receiveGhostLayers(i_wavePropagationBlock, i_topNeighborRank, o_topInflow, i_mpiRows);
		  sendGhostLayers(i_wavePropagationBlock, i_t, i_dtCoarse, i_timeSteppingStrategy, i_topNeighborRank, i_topOutflow, i_mpiSendType[BND_TOP], BND_TOP);
	  }

	  // if top neighbour is coarser, send and then receive
	  if (i_topNeighbourRefinementLevel <= i_refinementLevel) {
		  sendGhostLayers(i_wavePropagationBlock, i_t, i_dtCoarse, i_timeSteppingStrategy, i_topNeighborRank, i_topOutflow, i_mpiSendType[BND_TOP], BND_TOP);
		  receiveGhostLayers(i_wavePropagationBlock, i_topNeighborRank, o_topInflow, i_mpiRows);
	  }

	  // if bottom neighbour is finer, receive and then send
	  if (i_bottomNeighbourRefinementLevel > i_refinementLevel) {
		  receiveGhostLayers(i_wavePropagationBlock, i_bottomNeighborRank, o_bottomInflow, i_mpiRows);
		  sendGhostLayers(i_wavePropagationBlock, i_t, i_dtCoarse, i_timeSteppingStrategy, i_bottomNeighborRank, i_bottomOutflow, i_mpiSendType[BND_BOTTOM], BND_BOTTOM);
	  }
}


// TODO: add doxy comments for functions

// send ghost layers
void sendGhostLayers(SWE_WavePropagationAMR* i_wavePropagationBlock,
					  const float i_t,
					  const float i_dtCoarse,
					  const TimeSteppingType i_timeSteppingStrategy,
					  const int i_neighborRank,
					  SWE_BlockGhost* i_outflow,
					  MPI_Datatype i_mpiSendType,
					  BoundaryEdge i_edge) {
	i_wavePropagationBlock->synchCopyLayerBeforeRead(i_timeSteppingStrategy, i_edge, i_t, i_dtCoarse);
	MPI_Send( i_outflow->h.elemVector(),	1,	i_mpiSendType,	i_neighborRank,  1, MPI_COMM_WORLD );
	MPI_Send( i_outflow->hu.elemVector(),	1,	i_mpiSendType,	i_neighborRank,  2, MPI_COMM_WORLD );
	MPI_Send( i_outflow->hv.elemVector(),	1,	i_mpiSendType,	i_neighborRank,  3, MPI_COMM_WORLD );
}

// send ghost layers using Isend
void isendGhostLayers(SWE_WavePropagationAMR* i_wavePropagationBlock,
					  const float i_t,
					  const float i_dtCoarse,
					  const TimeSteppingType i_timeSteppingStrategy,
					  const int i_neighborRank,
					  SWE_BlockGhost* i_outflow,
					  MPI_Datatype i_mpiSendType,
					  BoundaryEdge i_edge,
					  MPI_Request* i_request) {
	i_wavePropagationBlock->synchCopyLayerBeforeRead(i_timeSteppingStrategy, i_edge, i_t, i_dtCoarse);
	MPI_Isend( i_outflow->h.elemVector(),	1,	i_mpiSendType,	i_neighborRank,  1, MPI_COMM_WORLD, i_request );
	MPI_Isend( i_outflow->hu.elemVector(),	1,	i_mpiSendType,	i_neighborRank,  2, MPI_COMM_WORLD, i_request );
	MPI_Isend( i_outflow->hv.elemVector(),	1,	i_mpiSendType,	i_neighborRank,  3, MPI_COMM_WORLD, i_request );
}

// receive ghost layers
void receiveGhostLayers(SWE_WavePropagationAMR* i_wavePropagationBlock,
						  const int i_neighborRank,
						  SWE_BlockGhost* o_inflow,
						  const MPI_Datatype i_mpiRecvType) {
	MPI_Status l_status;

	MPI_Recv( o_inflow->h.elemVector(),		1, i_mpiRecvType,	i_neighborRank,  1, MPI_COMM_WORLD, &l_status );
	MPI_Recv( o_inflow->hu.elemVector(),	1, i_mpiRecvType, 	i_neighborRank,  2, MPI_COMM_WORLD, &l_status );
	MPI_Recv( o_inflow->hv.elemVector(),	1, i_mpiRecvType, 	i_neighborRank,  3, MPI_COMM_WORLD, &l_status );
}

