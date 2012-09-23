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

#include "SWE_BlockAMR.hh"
#include "tools/help.hh"
#include <math.h>
#include <iostream>
#include <cassert>

/**
 * constructor
 * @param _offsetX
 * @param _offsetY
 * @param _nx
 * @param _ny
 * @param _dx
 * @param _dy
 * @param _rx
 * @param _ry
 * @param _nghosts
 * @param _interpolationStrategy
 */
SWE_BlockAMR::SWE_BlockAMR(float i_offsetX,
						   float i_offsetY,
						   int i_nx,
						   int i_ny,
						   float i_dx,
						   float i_dy,
						   int i_rx,
						   int i_ry,
						   int i_nghosts,
						   InterpolationType i_interpolationStrategy):
						   SWE_Block(i_offsetX, i_offsetY, i_nx, i_ny, i_dx, i_dy, i_nghosts),
						   rx(i_rx),
						   ry(i_ry),
						   interpolationStrategy(i_interpolationStrategy) {
	resetComputationalDomainMax();
	block_neighbour[BND_LEFT] = NULL;
	block_neighbour[BND_RIGHT] = NULL;
	block_neighbour[BND_BOTTOM] = NULL;
	block_neighbour[BND_TOP] = NULL;
}

void SWE_BlockAMR::initScenario(SWE_Scenario &i_scenario,
		const bool i_multipleBlocks) {
	SWE_Block::initScenario(i_scenario, i_multipleBlocks);
	scene = &i_scenario;
}

/**
 * Stores a pointer to the neighbouring block
 * TODO: change this!! will not work for MPI
 */
void SWE_BlockAMR::setBlockNeighbour(SWE_BlockAMR* i_neighbour, BoundaryEdge i_edge) {
	block_neighbour[i_edge] = i_neighbour;
}

/**
 * @param i_edge
 * @return
 */
SWE_BlockGhost* SWE_BlockAMR::registerCopyLayer(BoundaryEdge i_edge) {
	// for same resolution blocks, use proxy copy layer
	if (block_neighbour[i_edge]->getRefinementLevel() == getRefinementLevel()) {
		copyLayer[i_edge] = SWE_Block::registerCopyLayer(i_edge);
		proxyCopyLayer[i_edge] = copyLayer[i_edge];
		return copyLayer[i_edge];
	}

	// neighbour sizes
	int l_nx = block_neighbour[i_edge]->getNx();
	int l_ny = block_neighbour[i_edge]->getNy();
	int l_nghosts = block_neighbour[i_edge]->getNghosts();
	int l_xCopy, l_yCopy; //size of copy layer

	// get the size of the copy layer depending on the edge
	switch (i_edge) {
	case BND_LEFT:
	case BND_RIGHT:
		l_xCopy = l_nghosts;
		l_yCopy = l_ny + 2 * l_nghosts;
		break;
	case BND_BOTTOM:
	case BND_TOP:
		l_xCopy = l_nx + 2 * l_nghosts;
		l_yCopy = l_nghosts;
		break;
	}

	// in case the neighbour on the edge is finer than the current block
	// create the necessary blocks for higher order reconstruction
	if (block_neighbour[i_edge]->getRefinementLevel() > getRefinementLevel()) {
		if (interpolationStrategy != SPACE) {
			// for refining using time interpolation,
			// keep both copy layers at beginning and end of coarse time step
			startCopyLayer[i_edge] = new SWE_BlockGhost(l_xCopy, l_yCopy);
			endCopyLayer[i_edge] = new SWE_BlockGhost(l_xCopy, l_yCopy);
		}
		if (interpolationStrategy == APPROX_TIME_SPACE) {
			// for refining using approximate time interpolation
			// keep delta - the difference between the copyLayer coarse
			//			values at end and beginning of time-step extended
			//			by a piecewise-constant scheme to the fine dimensions
			delta[i_edge] = new SWE_BlockGhost(l_xCopy, l_yCopy);
		}
	}
	// initialise the proxy copy layer that will be coarsened/refined at each time-step
	proxyCopyLayer[i_edge] = getProxyCopyLayer(i_edge);

	// coarsened/refined copy layer from which the ghost values of the neighbouring block are copied
	copyLayer[i_edge] = new SWE_BlockGhost(l_xCopy, l_yCopy);
	return copyLayer[i_edge];
}

/**
 * @param i_edge
 * @return
 */
SWE_BlockGhost* SWE_BlockAMR::getProxyCopyLayer(BoundaryEdge i_edge) {
	// depending on the interpolation strategy of the ghost layer,
	// there are different sizes for the ghost layers,
	// therefore different proxy copy layers
	switch (interpolationStrategy) {
	case APPROX_TIME_SPACE:
		if (block_neighbour[i_edge]->getRefinementLevel() > getRefinementLevel())
			return getFineProxyCopyLayer_singleLayer(i_edge);
		else
			return getCoarseProxyCopyLayer_singleLayer(i_edge);
	case TIME_SPACE:
	case SPACE:
		if (block_neighbour[i_edge]->getRefinementLevel() > getRefinementLevel())
			return getFineProxyCopyLayer_multiLayer(i_edge);
		else
			return getCoarseProxyCopyLayer_multiLayer(i_edge);
	}
}

// ==========================================================================
//    methods to grab a proxy object of the copy layer needed for refining
// ==========================================================================
/**
 * for the APPROX_TIME_SPACE interpolation, there is only one ghost layer
 */
SWE_BlockGhost* SWE_BlockAMR::getFineProxyCopyLayer_singleLayer(BoundaryEdge i_edge) {
	int c, r, nc, nr;

	switch (i_edge) {
	case BND_LEFT: c = r = 0; nc = 3; nr = ny + 2; break;
	case BND_RIGHT: c = nx - 1; r = 0; nc = 3; nr = ny + 2; break;
	case BND_BOTTOM: c = r = 0; nc = nx + 2; nr = 3; break;
	case BND_TOP: c = 0; r = ny - 1; nc = nx + 2; nr = 3; break;
	}
	return new SWE_BlockGhost(*h.getBlockProxy(c, r, nc, nr), *b.getBlockProxy(c, r, nc, nr),
							  *hu.getBlockProxy(c, r, nc, nr), *hv.getBlockProxy(c, r, nc, nr));
}

SWE_BlockGhost* SWE_BlockAMR::getFineProxyCopyLayer_multiLayer(BoundaryEdge i_edge) {
	int l_nx = block_neighbour[i_edge]->getNx();
	int l_ny = block_neighbour[i_edge]->getNy();
	int l_nghosts = block_neighbour[i_edge]->getNghosts();
	int l_rx = block_neighbour[i_edge]->getRefinementLevel() / getRefinementLevel();
	int l_ry = block_neighbour[i_edge]->getRefinementLevel() / getRefinementLevel();
	int c, r, nc, nr;

	switch (i_edge) {
	case BND_LEFT:
		c = nghosts - 1;
		r = 0;
		nc = l_nghosts/l_rx + 2;
		nr = ny + 2*nghosts;
		break;
	case BND_RIGHT:
		c = nx+nghosts - l_nghosts/l_rx - 1;
		r = 0;
		nc = l_nghosts/l_rx + 2;
		nr = ny + 2*nghosts;
		break;
	case BND_BOTTOM:
		c = 0;
		r = nghosts - 1;
		nc = nx + 2*nghosts;
		nr = l_nghosts/l_ry + 2;
		break;
	case BND_TOP:
		c = 0;
		r = ny+nghosts - l_nghosts/l_ry - 1;
		nc = nx + 2*nghosts;
		nr = l_nghosts/l_ry + 2;
		break;
	}

	return new SWE_BlockGhost(*h.getBlockProxy(c, r, nc, nr), *b.getBlockProxy(c, r, nc, nr),
							  *hu.getBlockProxy(c, r, nc, nr), *hv.getBlockProxy(c, r, nc, nr));
}

// ==========================================================================
//    methods to grab a proxy object of the copy layer needed for coarsening
// ==========================================================================
/**
 * for the APPROX_TIME_SPACE interpolation, there is only one ghost layer
 */
SWE_BlockGhost* SWE_BlockAMR::getCoarseProxyCopyLayer_singleLayer(BoundaryEdge i_edge) {
	int c, r, nc, nr;
	int l_rx = getRefinementLevel() / block_neighbour[i_edge]->getRefinementLevel();
	int l_ry = getRefinementLevel() / block_neighbour[i_edge]->getRefinementLevel();

	switch (i_edge) {
	case BND_LEFT: c = r = 1; nc = l_rx; nr = ny; break;
	case BND_RIGHT: c = nx+1-l_rx; r = 1; nc = l_rx; nr = ny; break;
	case BND_BOTTOM: c = r = 1; nc = nx; nr = l_ry; break;
	case BND_TOP: c = 1; r = ny+1-l_ry; nc = nx; nr = l_ry; break;
	}

	return new SWE_BlockGhost(*h.getBlockProxy(c, r, nc, nr), *b.getBlockProxy(c, r, nc, nr),
							  *hu.getBlockProxy(c, r, nc, nr), *hv.getBlockProxy(c, r, nc, nr));
}

SWE_BlockGhost* SWE_BlockAMR::getCoarseProxyCopyLayer_multiLayer(BoundaryEdge i_edge) {
	int l_ny = block_neighbour[i_edge]->getNy();
	int l_nghosts = block_neighbour[i_edge]->getNghosts();
	int c, r, nc, nr;
	int l_rx = getRefinementLevel() / block_neighbour[i_edge]->getRefinementLevel();
	int l_ry = getRefinementLevel() / block_neighbour[i_edge]->getRefinementLevel();

	switch (i_edge) {
	case BND_LEFT:
		c = nghosts;
		r = 0;
		nc = l_nghosts*l_rx;
		nr = ny+2*nghosts;
		break;
	case BND_RIGHT:
		c = nx+nghosts-l_nghosts*l_rx;
		r = 0;
		nc = l_nghosts*l_rx;
		nr = ny + 2*nghosts;
		break;
	case BND_BOTTOM:
		c = 0;
		r = nghosts;
		nc = nx + 2*nghosts;
		nr = l_nghosts*l_ry;
		break;
	case BND_TOP:
		c = 0;
		r = ny+nghosts-l_nghosts*l_ry;
		nc = nx + 2*nghosts;
		nr = l_nghosts*l_ry;
		break;
	}

	return new SWE_BlockGhost(*h.getBlockProxy(c, r, nc, nr), *b.getBlockProxy(c, r, nc, nr),
							  *hu.getBlockProxy(c, r, nc, nr), *hv.getBlockProxy(c, r, nc, nr));
}


/**
 * @param i_timeStepping
 * @param i_edge
 * @param t
 * @param dt_coarse
 * @return
 */
void SWE_BlockAMR::synchCopyLayerBeforeRead(TimeSteppingType i_timeStepping,
											BoundaryEdge i_edge,
											float t,
											float dt_coarse) {
	if (boundary[i_edge] == CONNECT) {
		// do nothing for same resolution neighbour
		if (block_neighbour[i_edge]->getRefinementLevel() == getRefinementLevel())
			return;
		// neighbour is coarser
		if (block_neighbour[i_edge]->getRefinementLevel() < getRefinementLevel())
			setCopyLayerCoarse(i_edge);
		else {
			// neighbour is finer
			if (i_timeStepping == GTS)
				// for global time-stepping scheme
				// perform higher order reconstruction
				setCopyLayerFine(i_edge, copyLayer[i_edge]);
			else if (i_timeStepping == LTS) {
				// for local time-stepping scheme
				// do time interpolation ([APPROX_]TIME_SPACE) or
				// decrease the computational domain (SPACE) -> in block manager
				if (interpolationStrategy == SPACE)
					setCopyLayerFine(i_edge, copyLayer[i_edge]);
				else
					timeInterpolateCopyLayer(i_edge, t, dt_coarse);
			}
		}
	}
}

/**
 * @return
 */
void SWE_BlockAMR::synchBeforeRead() {
	for (int l_edge = 0; l_edge < 4; ++l_edge)
		if (boundary[l_edge] == CONNECT && block_neighbour[l_edge]->getRefinementLevel() > getRefinementLevel()) {
			setCopyLayerFine(l_edge, startCopyLayer[l_edge]);
			int l_rx = block_neighbour[l_edge]->getRefinementLevel()/getRefinementLevel();
			int l_ry = block_neighbour[l_edge]->getRefinementLevel()/getRefinementLevel();

			// APPROX_TIME_SPACE: set delta
			// the difference between coarse value at beginning and end of time-step
			if (interpolationStrategy == APPROX_TIME_SPACE) {
				switch (l_edge) {
				case BND_LEFT:
					for (int i=1; i<=ny; ++i)
						for (int l=(i-1)*l_ry+1; l<=i*l_ry; ++l) {
							delta[l_edge]->h[0][l] = -h[1][i];
							delta[l_edge]->hu[0][l] = -hu[1][i];
							delta[l_edge]->hv[0][l] = -hv[1][i];
						}
					break;
				case BND_RIGHT:
					for (int i=1; i<=ny; ++i)
						for (int l=(i-1)*l_ry+1; l<=i*l_ry; ++l) {
							delta[l_edge]->h[0][l] = -h[nx][i];
							delta[l_edge]->hu[0][l] = -hu[nx][i];
							delta[l_edge]->hv[0][l] = -hv[nx][i];
						}
					break;
				case BND_BOTTOM:
					for (int i=1; i<=nx; ++i)
						for (int k=(i-1)*l_rx+1; k<=i*l_rx; ++k) {
							delta[l_edge]->h[k][0] = -h[i][1];
							delta[l_edge]->hu[k][0] = -hu[i][1];
							delta[l_edge]->hv[k][0] = -hv[i][1];
						}
					break;
				case BND_TOP:
					for (int i=1; i<=nx; ++i)
						for (int k=(i-1)*l_rx+1; k<=i*l_rx; ++k) {
							delta[l_edge]->h[k][0] = -h[i][ny];
							delta[l_edge]->hu[k][0] = -hu[i][ny];
							delta[l_edge]->hv[k][0] = -hv[i][ny];
						}
					break;
				}
			}
		}
}

/**
 *
 */
void SWE_BlockAMR::synchAfterWrite() {
	for (int l_edge = 0; l_edge < 4; ++l_edge)
		if (boundary[l_edge] == CONNECT && block_neighbour[l_edge]->getRefinementLevel() > getRefinementLevel()) {
			int l_rx = block_neighbour[l_edge]->getRefinementLevel()/getRefinementLevel();
			int l_ry = block_neighbour[l_edge]->getRefinementLevel()/getRefinementLevel();

			// TIME_SPACE: refine copy layer at end of time-step
			if (interpolationStrategy == TIME_SPACE)
				setCopyLayerFine(l_edge, endCopyLayer[l_edge]);

			// APPROX_TIME_SPACE: set delta
			// the difference between coarse value at beginning and end of time-step
			if (interpolationStrategy == APPROX_TIME_SPACE) {
				switch (l_edge) {
				case BND_LEFT:
					for (int i=1; i<=ny; ++i)
						for (int l=(i-1)*l_ry+1; l<=i*l_ry; ++l) {
							delta[l_edge]->h[0][l] += h[1][i];
							delta[l_edge]->hu[0][l] += hu[1][i];
							delta[l_edge]->hv[0][l] += hv[1][i];
						}
					break;
				case BND_RIGHT:
					for (int i=1; i<=ny; ++i)
						for (int l=(i-1)*l_ry+1; l<=i*l_ry; ++l) {
							delta[l_edge]->h[0][l] += h[nx][i];
							delta[l_edge]->hu[0][l] += hu[nx][i];
							delta[l_edge]->hv[0][l] += hv[nx][i];
						}
					break;
				case BND_BOTTOM:
					for (int i=1; i<=nx; ++i)
						for (int k=(i-1)*l_rx+1; k<=i*l_rx; ++k) {
							delta[l_edge]->h[k][0] += h[i][1];
							delta[l_edge]->hu[k][0] += hu[i][1];
							delta[l_edge]->hv[k][0] += hv[i][1];
						}
					break;
				case BND_TOP:
					for (int i=1; i<=nx; ++i)
						for (int k=(i-1)*l_rx+1; k<=i*l_rx; ++k) {
							delta[l_edge]->h[k][0] += h[i][ny];
							delta[l_edge]->hu[k][0] += hu[i][ny];
							delta[l_edge]->hv[k][0] += hv[i][ny];
						}
					break;
				}
			}
		}
}

/**
 * @param i_edge
 * @param t
 * @param dt_coarseoarse
 */
void SWE_BlockAMR::timeInterpolateCopyLayer(BoundaryEdge i_edge, float t, float dt_coarse) {
	if (interpolationStrategy == APPROX_TIME_SPACE)
		for (int i = 0; i < copyLayer[i_edge]->nx; ++i)
			for (int j = 0; j < copyLayer[i_edge]->ny; ++j) {
				copyLayer[i_edge]->h[i][j] = startCopyLayer[i_edge]->h[i][j] + delta[i_edge]->h[i][j] * t/dt_coarse;
				copyLayer[i_edge]->hu[i][j] = startCopyLayer[i_edge]->hu[i][j] + delta[i_edge]->h[i][j] * t/dt_coarse;
				copyLayer[i_edge]->hv[i][j] = startCopyLayer[i_edge]->hv[i][j] + delta[i_edge]->h[i][j] * t/dt_coarse;
			}
	else if (interpolationStrategy == TIME_SPACE)
		for (int i = 0; i < copyLayer[i_edge]->nx; ++i)
			for (int j = 0; j < copyLayer[i_edge]->ny; ++j) {
				copyLayer[i_edge]->h[i][j] = startCopyLayer[i_edge]->h[i][j]
						+ (endCopyLayer[i_edge]->h[i][j] - startCopyLayer[i_edge]->h[i][j]) * t / dt_coarse;
				copyLayer[i_edge]->hu[i][j] = startCopyLayer[i_edge]->hu[i][j]
						+ (endCopyLayer[i_edge]->hu[i][j] - startCopyLayer[i_edge]->hu[i][j]) * t / dt_coarse;
				copyLayer[i_edge]->hv[i][j] = startCopyLayer[i_edge]->hv[i][j]
						+ (endCopyLayer[i_edge]->hv[i][j] - startCopyLayer[i_edge]->hv[i][j]) * t / dt_coarse;
			}
}

/**
 * @param i_edge
 */
void SWE_BlockAMR::setCopyLayerCoarse(BoundaryEdge i_edge) {
	int l_nx = copyLayer[i_edge]->nx;
	int l_ny = copyLayer[i_edge]->ny;
	int l_rx = getRefinementLevel() / block_neighbour[i_edge]->getRefinementLevel();
	int l_ry = getRefinementLevel() / block_neighbour[i_edge]->getRefinementLevel();

	SWE_BlockGhost* tmp = proxyCopyLayer[i_edge]->coarsen(l_rx, l_ry);

	// copy the result in the copy layer
	if (interpolationStrategy == APPROX_TIME_SPACE) { // single layer
		switch (i_edge) {
		case BND_LEFT:
		case BND_RIGHT:
			for (int i=1; i<l_ny-1; i++) {
				copyLayer[i_edge]->h[0][i] = tmp->h[0][i-1];
				copyLayer[i_edge]->hu[0][i] = tmp->hu[0][i-1];
				copyLayer[i_edge]->hv[0][i] = tmp->hv[0][i-1];
			}
			break;
		case BND_BOTTOM:
		case BND_TOP:
			for (int i=1; i<l_nx-1; i++) {
				copyLayer[i_edge]->h[i][0] = tmp->h[i-1][0];
				copyLayer[i_edge]->hu[i][0] = tmp->hu[i-1][0];
				copyLayer[i_edge]->hv[i][0] = tmp->hv[i-1][0];
			}
			break;
		}
	} else { // multi layer
		// the result is always the same size as the copy layer
		for (int i=0; i<l_nx; i++)
			for (int j=0; j<l_ny; j++) {
				copyLayer[i_edge]->h[i][j] = tmp->h[i][j];
				copyLayer[i_edge]->hu[i][j] = tmp->hu[i][j];
				copyLayer[i_edge]->hv[i][j] = tmp->hv[i][j];
			}
	}

}

/**
 * @param i_edge
 * @param layer
 */
void SWE_BlockAMR::setCopyLayerFine(BoundaryEdge i_edge, SWE_BlockGhost* layer) {
	int l_nx = layer->nx;
	int l_ny = layer->ny;
	int l_rx = block_neighbour[i_edge]->getRefinementLevel() / getRefinementLevel();
	int l_ry = block_neighbour[i_edge]->getRefinementLevel() / getRefinementLevel();
	int l_nghosts = block_neighbour[i_edge]->getNghosts();

//	SWE_BlockGhost* tmp = proxyCopyLayer[i_edge]->refine_constant(l_rx, l_ry);

	SWE_BlockGhost* tmp;
	switch (i_edge) {
	case BND_LEFT: tmp = proxyCopyLayer[i_edge]->refine(l_rx, l_ry, dx, dy, scene, offsetX, offsetY+dy); break;
	case BND_RIGHT: tmp = proxyCopyLayer[i_edge]->refine(l_rx, l_ry, dx, dy, scene, offsetX+nx*dx, offsetY+dy); break;
	case BND_BOTTOM: tmp = proxyCopyLayer[i_edge]->refine(l_rx, l_ry, dx, dy, scene, offsetX+dx, offsetY); break;
	case BND_TOP: tmp = proxyCopyLayer[i_edge]->refine(l_rx, l_ry, dx, dy, scene, offsetX+dx, offsetY+ny*dy); break;
	}

	// copy the result in the copy layer
	if (interpolationStrategy == APPROX_TIME_SPACE) { // single layer
		switch (i_edge) {
		case BND_LEFT:
			for (int i=1; i<l_ny-1; i++) {
				layer->h[0][i] = tmp->h[0][i-1];
				layer->hu[0][i] = tmp->hu[0][i-1];
				layer->hv[0][i] = tmp->hv[0][i-1];
			}
			break;
		case BND_RIGHT:
			for (int i=1; i<l_ny-1; i++) {
				layer->h[0][i] = tmp->h[l_rx-1][i-1];
				layer->hu[0][i] = tmp->hu[l_rx-1][i-1];
				layer->hv[0][i] = tmp->hv[l_rx-1][i-1];
			}
			break;
		case BND_BOTTOM:
			for (int i=1; i<l_nx-1; i++) {
				layer->h[i][0] = tmp->h[i-1][0];
				layer->hu[i][0] = tmp->hu[i-1][0];
				layer->hv[i][0] = tmp->hv[i-1][0];
			}
			break;
		case BND_TOP:
			for (int i=1; i<l_nx-1; i++) {
				layer->h[i][0] = tmp->h[i-1][l_ry-1];
				layer->hu[i][0] = tmp->hu[i-1][l_ry-1];
				layer->hv[i][0] = tmp->hv[i-1][l_ry-1];
			}
			break;
		}
	} else { // multi layer -> the result is always 2*l_nghosts elements smaller in one direction
		switch (i_edge) {
		case BND_LEFT:
		case BND_RIGHT:
			for (int i=0; i<l_nx; i++)
				for (int j=l_nghosts; j<l_ny-l_nghosts; j++) {
					layer->h[i][j] = tmp->h[i][j-l_nghosts];
					layer->hu[i][j] = tmp->hu[i][j-l_nghosts];
					layer->hv[i][j] = tmp->hv[i][j-l_nghosts];
				}
			break;
		case BND_BOTTOM:
		case BND_TOP:
			for (int i=l_nghosts; i<l_nx-l_nghosts; i++)
				for (int j=0; j<l_ny; j++) {
					layer->h[i][j] = tmp->h[i-l_nghosts][j];
					layer->hu[i][j] = tmp->hu[i-l_nghosts][j];
					layer->hv[i][j] = tmp->hv[i-l_nghosts][j];
				}
			break;
		}
	}
}

/**
 *
 */
void SWE_BlockAMR::decreaseComputationalDomain() {
	if (boundary[BND_LEFT] == CONNECT) nxint_s++;
	if (boundary[BND_BOTTOM] == CONNECT) nyint_s++;
	if (boundary[BND_RIGHT] == CONNECT) nxint_e--;
	if (boundary[BND_TOP] == CONNECT) nyint_e--;
}

/**
 *
 */
void SWE_BlockAMR::resetComputationalDomainMax() {
	if (nghosts == 1) {
		nxint_s = nyint_s = 1;
		nxint_e = nx + 2 * nghosts - 2;
		nyint_e = ny + 2 * nghosts - 2;
	} else {
		nxint_s = nyint_s = nghosts / 2;
		nxint_e = nx + 3 * nghosts / 2 - 1;
		nyint_e = ny + 3 * nghosts / 2 - 1;
	}
}

/**
 *
 */
void SWE_BlockAMR::resetComputationalDomainMin() {
	nxint_s = nyint_s = nghosts;
	nxint_e = nx + nghosts - 1;
	nyint_e = ny + nghosts - 1;
}

/**
 * set the values of ghost cells for one edge depending on the specifed boundary conditions;
 * if the ghost layer replicates the variables of a remote SWE_Block, the values are copied
 *
 * the function has a BoundaryEdge parameter to make possible the implementation of the
 * two-phase ghost layer exchange (to update the corner ghost cells properly)
 *
 * @param i_edge	the boundary edge where the ghost layers must be set
 */
void SWE_BlockAMR::setGhostLayerEdge(BoundaryEdge i_edge) {

#ifdef DBG
	cout << "Set simple boundary conditions " << endl << flush;
#endif
	// call to virtual function to set ghost layer values
	setBoundaryConditions();

	// for a CONNECT boundary, data will be copied from a neighbouring
	// SWE_Block (via a SWE_BlockGhost proxy object)
	// -> these copy operations cannot be executed in GPU/accelerator memory, e.g.
	//    setBoundaryConditions then has to take care that values are copied.

#ifdef DBG
	cout << "Set CONNECT boundary conditions in main memory " << endl << flush;
#endif

	switch (i_edge) {
	case BND_LEFT:
		// left boundary
		if (boundary[BND_LEFT] == CONNECT) {
			for (int i = 0; i < nghosts; i++) {
				for (int j = 0; j < ny + 2 * nghosts; j++) {
					h[i][j] = neighbour[BND_LEFT]->h[i][j];
					hu[i][j] = neighbour[BND_LEFT]->hu[i][j];
					hv[i][j] = neighbour[BND_LEFT]->hv[i][j];
				}
			}
		}

	case BND_RIGHT:
		// right boundary
		if (boundary[BND_RIGHT] == CONNECT) {
			for (int i = 0; i < nghosts; i++) {

				for (int j = 0; j < ny + 2 * nghosts; j++) {
					h[nx + nghosts + i][j] = neighbour[BND_RIGHT]->h[i][j];
					hu[nx + nghosts + i][j] = neighbour[BND_RIGHT]->hu[i][j];
					hv[nx + nghosts + i][j] = neighbour[BND_RIGHT]->hv[i][j];
				}
			}
		}

	case BND_BOTTOM:
		// bottom boundary
		if (boundary[BND_BOTTOM] == CONNECT) {
			for (int i = 0; i < nx + 2 * nghosts; i++) {
				for (int j = 0; j < nghosts; j++) {
					h[i][j] = neighbour[BND_BOTTOM]->h[i][j];
					hu[i][j] = neighbour[BND_BOTTOM]->hu[i][j];
					hv[i][j] = neighbour[BND_BOTTOM]->hv[i][j];
				}
			}
		}

	case BND_TOP:
		// top boundary
		if (boundary[BND_TOP] == CONNECT) {
			for (int i = 0; i < nx + 2 * nghosts; i++) {
				for (int j = 0; j < nghosts; j++) {
					h[i][ny + nghosts + j] = neighbour[BND_TOP]->h[i][j];
					hu[i][ny + nghosts + j] = neighbour[BND_TOP]->hu[i][j];
					hv[i][ny + nghosts + j] = neighbour[BND_TOP]->hv[i][j];
				}
			}
		}
	}

#ifdef DBG
	cout << "Synchronize ghost layers (for heterogeneous memory) " << endl << flush;
#endif
	// synchronize the ghost layers (for PASSIVE and CONNECT conditions)
	// with accelerator memory
	synchGhostLayerAfterWrite();
}

