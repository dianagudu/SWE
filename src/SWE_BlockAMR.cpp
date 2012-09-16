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
SWE_BlockAMR::SWE_BlockAMR(float _offsetX, float _offsetY, int _nx,
		int _ny, float _dx, float _dy, int _rx, int _ry, int _nghosts,
		InterpolationType _interpolationStrategy) :
		SWE_Block(_offsetX, _offsetY, _nx, _ny, _dx, _dy, _nghosts),
		rx(_rx), ry(_ry), interpolationStrategy(_interpolationStrategy) {
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
void SWE_BlockAMR::setBlockNeighbour(SWE_BlockAMR* neigh, BoundaryEdge edge) {
	block_neighbour[edge] = neigh;
}

/**
 * @param edge
 * @return
 */
SWE_BlockGhost* SWE_BlockAMR::registerCopyLayer(BoundaryEdge edge) {
	// for same resolution blocks, use proxy copy layer
	if (block_neighbour[edge]->getRefinementLevel() == getRefinementLevel()) {
		copyLayer[edge] = SWE_Block::registerCopyLayer(edge);
		proxyCopyLayer[edge] = copyLayer[edge];
		return copyLayer[edge];
	}

	// neighbour sizes
	int _nx = block_neighbour[edge]->getNx();
	int _ny = block_neighbour[edge]->getNy();
	int _nghosts = block_neighbour[edge]->getNghosts();
	int xCopy, yCopy; //size of copy layer

	// get the size of the copy layer depending on the edge
	switch (edge) {
	case BND_LEFT: case BND_RIGHT:
		xCopy = _nghosts; yCopy = _ny + 2 * _nghosts; break;
	case BND_BOTTOM: case BND_TOP:
		xCopy = _nx + 2 * _nghosts; yCopy = _nghosts; break;
	}

	// in case the neighbour on the edge is finer than the current block
	// create the necessary blocks for higher order reconstruction
	if (block_neighbour[edge]->getRefinementLevel() > getRefinementLevel()) {
		if (interpolationStrategy != SPACE) {
			// for refining using time interpolation,
			// keep both copy layers at beginning and end of coarse time step
			startCopyLayer[edge] = new SWE_BlockGhost(xCopy, yCopy);
			endCopyLayer[edge] = new SWE_BlockGhost(xCopy, yCopy);
		}
		if (interpolationStrategy == APPROX_TIME_SPACE) {
			// for refining using approximate time interpolation
			// keep delta - the difference between the copyLayer coarse
			//			values at end and beginning of time-step extended
			//			by a piecewise-constant scheme to the fine dimensions
			delta[edge] = new SWE_BlockGhost(xCopy, yCopy);
		}
	}
	// initialise the proxy copy layer that will be coarsened/refined at each time-step
	proxyCopyLayer[edge] = getProxyCopyLayer(edge);

	// coarsened/refined copy layer from which the ghost values of the neighbouring block are copied
	copyLayer[edge] = new SWE_BlockGhost(xCopy, yCopy);
	return copyLayer[edge];
}

/**
 * @param edge
 * @return
 */
SWE_BlockGhost* SWE_BlockAMR::getProxyCopyLayer(BoundaryEdge edge) {
	// depending on the interpolation strategy of the ghost layer,
	// there are different sizes for the ghost layers,
	// therefore different proxy copy layers
	switch (interpolationStrategy) {
	case APPROX_TIME_SPACE:
		if (block_neighbour[edge]->getRefinementLevel() > getRefinementLevel())
			return getFineProxyCopyLayer_singleLayer(edge);
		else
			return getCoarseProxyCopyLayer_singleLayer(edge);
	case TIME_SPACE:
	case SPACE:
		if (block_neighbour[edge]->getRefinementLevel() > getRefinementLevel())
			return getFineProxyCopyLayer_multiLayer(edge);
		else
			return getCoarseProxyCopyLayer_multiLayer(edge);
	}
}

// ==========================================================================
//    methods to grab a proxy object of the copy layer needed for refining
// ==========================================================================
/**
 * for the APPROX_TIME_SPACE interpolation, there is only one ghost layer
 */
SWE_BlockGhost* SWE_BlockAMR::getFineProxyCopyLayer_singleLayer(BoundaryEdge edge) {
	int c, r, nc, nr;

	switch (edge) {
	case BND_LEFT: c = r = 0; nc = 3; nr = ny + 2; break;
	case BND_RIGHT: c = nx - 1; r = 0; nc = 3; nr = ny + 2; break;
	case BND_BOTTOM: c = r = 0; nc = nx + 2; nr = 3; break;
	case BND_TOP: c = 0; r = ny - 1; nc = nx + 2; nr = 3; break;
	}
	return new SWE_BlockGhost(*h.getBlockProxy(c, r, nc, nr), *b.getBlockProxy(c, r, nc, nr),
							  *hu.getBlockProxy(c, r, nc, nr), *hv.getBlockProxy(c, r, nc, nr));
}

SWE_BlockGhost* SWE_BlockAMR::getFineProxyCopyLayer_multiLayer(BoundaryEdge edge) {
	int _nx = block_neighbour[edge]->getNx();
	int _ny = block_neighbour[edge]->getNy();
	int _nghosts = block_neighbour[edge]->getNghosts();
	int _rx = block_neighbour[edge]->getRefinementLevel() / getRefinementLevel();
	int _ry = block_neighbour[edge]->getRefinementLevel() / getRefinementLevel();
	int c, r, nc, nr;

	switch (edge) {
	case BND_LEFT:
		c = nghosts - 1;
		r = 0;
		nc = _nghosts*_rx + 2;
		nr = ny + 2*nghosts;
		break;
	case BND_RIGHT:
		c = nx+nghosts - _nghosts*_rx - 1;
		r = 0;
		nc = _nghosts*_rx + 2;
		nr = ny + 2*nghosts;
		break;
	case BND_BOTTOM:
		c = 0;
		r = nghosts - 1;
		nc = nx + 2*nghosts;
		nr = _nghosts*_ry + 2;
		break;
	case BND_TOP:
		c = 0;
		r = ny+nghosts - _nghosts*_ry - 1;
		nc = nx + 2*nghosts;
		nr = _nghosts*_ry + 2;
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
SWE_BlockGhost* SWE_BlockAMR::getCoarseProxyCopyLayer_singleLayer(BoundaryEdge edge) {
	int c, r, nc, nr;
	int _rx = getRefinementLevel() / block_neighbour[edge]->getRefinementLevel();
	int _ry = getRefinementLevel() / block_neighbour[edge]->getRefinementLevel();

	switch (edge) {
	case BND_LEFT: c = r = 1; nc = _rx; nr = ny; break;
	case BND_RIGHT: c = nx+1-_rx; r = 1; nc = _rx; nr = ny; break;
	case BND_BOTTOM: c = r = 1; nc = nx; nr = _ry; break;
	case BND_TOP: c = 1; r = ny+1-_ry; nc = nx; nr = _ry; break;
	}

	return new SWE_BlockGhost(*h.getBlockProxy(c, r, nc, nr), *b.getBlockProxy(c, r, nc, nr),
							  *hu.getBlockProxy(c, r, nc, nr), *hv.getBlockProxy(c, r, nc, nr));
}

SWE_BlockGhost* SWE_BlockAMR::getCoarseProxyCopyLayer_multiLayer(BoundaryEdge edge) {
	int _ny = block_neighbour[edge]->getNy();
	int _nghosts = block_neighbour[edge]->getNghosts();
	int c, r, nc, nr;
	int _rx = getRefinementLevel() / block_neighbour[edge]->getRefinementLevel();
	int _ry = getRefinementLevel() / block_neighbour[edge]->getRefinementLevel();

	switch (edge) {
	case BND_LEFT:
		c = nghosts;
		r = 0;
		nc = _nghosts*_rx;
		nr = ny+2*nghosts;
		break;
	case BND_RIGHT:
		c = nx+nghosts-_nghosts*_rx;
		r = 0;
		nc = _nghosts*_rx;
		nr = ny + 2*nghosts;
		break;
	case BND_BOTTOM:
		c = 0;
		r = nghosts;
		nc = nx + 2*nghosts;
		nr = _nghosts*_ry;
		break;
	case BND_TOP:
		c = 0;
		r = ny+nghosts-_nghosts*ry;
		nc = nx + 2*nghosts;
		nr = _nghosts*_ry;
		break;
	}

	return new SWE_BlockGhost(*h.getBlockProxy(c, r, nc, nr), *b.getBlockProxy(c, r, nc, nr),
							  *hu.getBlockProxy(c, r, nc, nr), *hv.getBlockProxy(c, r, nc, nr));
}


/**
 * @param timeStepping
 * @param edge
 * @param t
 * @param dt_c
 * @return
 */
void SWE_BlockAMR::synchCopyLayerBeforeRead(TimeSteppingType timeStepping, BoundaryEdge edge,
												float t, float dt_c) {
	if (boundary[edge] == CONNECT) {
		// do nothing for same resolution neighbour
		if (block_neighbour[edge]->getRefinementLevel() == getRefinementLevel())
			return;
		// neighbour is coarser
		if (block_neighbour[edge]->getRefinementLevel() < getRefinementLevel())
			setCopyLayerCoarse(edge);
		else {
			// neighbour is finer
			if (timeStepping == GTS)
				// for global time-stepping scheme
				// perform higher order reconstruction
				setCopyLayerFine(edge, copyLayer[edge]);
			else if (timeStepping == LTS) {
				// for local time-stepping scheme
				// do time interpolation ([APPROX_]TIME_SPACE) or
				// decrease the computational domain (SPACE) -> in block manager
				if (interpolationStrategy == SPACE)
					setCopyLayerFine(edge, copyLayer[edge]);
				else
					timeInterpolateCopyLayer(edge, t, dt_c);
			}
		}
	}
}

/**
 * @return
 */
void SWE_BlockAMR::synchBeforeRead() {
	for (int edge = 0; edge < 4; ++edge)
		if (boundary[edge] == CONNECT && block_neighbour[edge]->getRefinementLevel() > getRefinementLevel()) {
			setCopyLayerFine(edge, startCopyLayer[edge]);
			int _rx = block_neighbour[edge]->getRefinementLevel()/getRefinementLevel();
			int _ry = block_neighbour[edge]->getRefinementLevel()/getRefinementLevel();

			// APPROX_TIME_SPACE: set delta
			// the difference between coarse value at beginning and end of time-step
			if (interpolationStrategy == APPROX_TIME_SPACE) {
				switch (edge) {
				case BND_LEFT:
					for (int i=1; i<=ny; ++i)
						for (int l=(i-1)*_ry+1; l<=i*_ry; ++l) {
							delta[edge]->h[0][l] = -h[1][i];
							delta[edge]->hu[0][l] = -hu[1][i];
							delta[edge]->hv[0][l] = -hv[1][i];
						}
					break;
				case BND_RIGHT:
					for (int i=1; i<=ny; ++i)
						for (int l=(i-1)*_ry+1; l<=i*_ry; ++l) {
							delta[edge]->h[0][l] = -h[nx][i];
							delta[edge]->hu[0][l] = -hu[nx][i];
							delta[edge]->hv[0][l] = -hv[nx][i];
						}
					break;
				case BND_BOTTOM:
					for (int i=1; i<=nx; ++i)
						for (int k=(i-1)*_rx+1; k<=i*_rx; ++k) {
							delta[edge]->h[k][0] = -h[i][1];
							delta[edge]->hu[k][0] = -hu[i][1];
							delta[edge]->hv[k][0] = -hv[i][1];
						}
					break;
				case BND_TOP:
					for (int i=1; i<=nx; ++i)
						for (int k=(i-1)*_rx+1; k<=i*_rx; ++k) {
							delta[edge]->h[k][0] = -h[i][ny];
							delta[edge]->hu[k][0] = -hu[i][ny];
							delta[edge]->hv[k][0] = -hv[i][ny];
						}
					break;
				}
			}
		}
}

/**
 * @return
 */
void SWE_BlockAMR::synchAfterWrite() {
	for (int edge = 0; edge < 4; ++edge)
		if (boundary[edge] == CONNECT && block_neighbour[edge]->getRefinementLevel() > getRefinementLevel()) {
			int _rx = block_neighbour[edge]->getRefinementLevel()/getRefinementLevel();
			int _ry = block_neighbour[edge]->getRefinementLevel()/getRefinementLevel();

			// TIME_SPACE: refine copy layer at end of time-step
			if (interpolationStrategy == TIME_SPACE)
				setCopyLayerFine(edge, endCopyLayer[edge]);

			// APPROX_TIME_SPACE: set delta
			// the difference between coarse value at beginning and end of time-step
			if (interpolationStrategy == APPROX_TIME_SPACE) {
				switch (edge) {
				case BND_LEFT:
					for (int i=1; i<=ny; ++i)
						for (int l=(i-1)*_ry+1; l<=i*_ry; ++l) {
							delta[edge]->h[0][l] += h[1][i];
							delta[edge]->hu[0][l] += hu[1][i];
							delta[edge]->hv[0][l] += hv[1][i];
						}
					break;
				case BND_RIGHT:
					for (int i=1; i<=ny; ++i)
						for (int l=(i-1)*_ry+1; l<=i*_ry; ++l) {
							delta[edge]->h[0][l] += h[nx][i];
							delta[edge]->hu[0][l] += hu[nx][i];
							delta[edge]->hv[0][l] += hv[nx][i];
						}
					break;
				case BND_BOTTOM:
					for (int i=1; i<=nx; ++i)
						for (int k=(i-1)*_rx+1; k<=i*_rx; ++k) {
							delta[edge]->h[k][0] += h[i][1];
							delta[edge]->hu[k][0] += hu[i][1];
							delta[edge]->hv[k][0] += hv[i][1];
						}
					break;
				case BND_TOP:
					for (int i=1; i<=nx; ++i)
						for (int k=(i-1)*_rx+1; k<=i*_rx; ++k) {
							delta[edge]->h[k][0] += h[i][ny];
							delta[edge]->hu[k][0] += hu[i][ny];
							delta[edge]->hv[k][0] += hv[i][ny];
						}
					break;
				}
			}
		}
}

/**
 * @param edge
 * @param t
 * @param dt_c
 */
void SWE_BlockAMR::timeInterpolateCopyLayer(BoundaryEdge edge, float t, float dt_c) {
	if (interpolationStrategy == APPROX_TIME_SPACE)
		for (int i = 0; i < copyLayer[edge]->nx; ++i)
			for (int j = 0; j < copyLayer[edge]->ny; ++j) {
				copyLayer[edge]->h[i][j] = startCopyLayer[edge]->h[i][j] + delta[edge]->h[i][j] * t/dt_c;
				copyLayer[edge]->hu[i][j] = startCopyLayer[edge]->hu[i][j] + delta[edge]->h[i][j] * t/dt_c;
				copyLayer[edge]->hv[i][j] = startCopyLayer[edge]->hv[i][j] + delta[edge]->h[i][j] * t/dt_c;
			}
	else if (interpolationStrategy == TIME_SPACE)
		for (int i = 0; i < copyLayer[edge]->nx; ++i)
			for (int j = 0; j < copyLayer[edge]->ny; ++j) {
				copyLayer[edge]->h[i][j] = startCopyLayer[edge]->h[i][j]
						+ (endCopyLayer[edge]->h[i][j] - startCopyLayer[edge]->h[i][j]) * t / dt_c;
				copyLayer[edge]->hu[i][j] = startCopyLayer[edge]->hu[i][j]
						+ (endCopyLayer[edge]->hu[i][j] - startCopyLayer[edge]->hu[i][j]) * t / dt_c;
				copyLayer[edge]->hv[i][j] = startCopyLayer[edge]->hv[i][j]
						+ (endCopyLayer[edge]->hv[i][j] - startCopyLayer[edge]->hv[i][j]) * t / dt_c;
			}
}

/**
 * @param edge
 */
void SWE_BlockAMR::setCopyLayerCoarse(BoundaryEdge edge) {
	int _nx = copyLayer[edge]->nx;
	int _ny = copyLayer[edge]->ny;
	int _rx = getRefinementLevel() / block_neighbour[edge]->getRefinementLevel();
	int _ry = getRefinementLevel() / block_neighbour[edge]->getRefinementLevel();

	SWE_BlockGhost* tmp = proxyCopyLayer[edge]->coarsen(_rx, _ry);

	// copy the result in the copy layer
	if (interpolationStrategy == APPROX_TIME_SPACE) { // single layer
		switch (edge) {
		case BND_LEFT:
			for (int i=1; i<_ny-1; i++) {
				copyLayer[edge]->h[0][i] = tmp->h[0][i-1];
				copyLayer[edge]->hu[0][i] = tmp->hu[0][i-1];
				copyLayer[edge]->hv[0][i] = tmp->hv[0][i-1];
			}
			break;
		case BND_RIGHT:
			for (int i=1; i<_ny-1; i++) {
				copyLayer[edge]->h[0][i] = tmp->h[0][i-1];
				copyLayer[edge]->hu[0][i] = tmp->hu[0][i-1];
				copyLayer[edge]->hv[0][i] = tmp->hv[0][i-1];
			}
			break;
		case BND_BOTTOM:
			for (int i=1; i<_nx-1; i++) {
				copyLayer[edge]->h[i][0] = tmp->h[i-1][0];
				copyLayer[edge]->hu[i][0] = tmp->hu[i-1][0];
				copyLayer[edge]->hv[i][0] = tmp->hv[i-1][0];
			}
			break;
		case BND_TOP:
			for (int i=1; i<_nx-1; i++) {
				copyLayer[edge]->h[i][0] = tmp->h[i-1][0];
				copyLayer[edge]->hu[i][0] = tmp->hu[i-1][0];
				copyLayer[edge]->hv[i][0] = tmp->hv[i-1][0];
			}
			break;
		}
	} else { // multi layer
		// the result is always the same size as the copy layer
		for (int i=0; i<_nx; i++)
			for (int j=0; j<_ny; j++) {
				copyLayer[edge]->h[i][j] = tmp->h[i][j];
				copyLayer[edge]->hu[i][j] = tmp->hu[i][j];
				copyLayer[edge]->hv[i][j] = tmp->hv[i][j];
			}
	}

}

/**
 * @param edge
 * @param layer
 */
void SWE_BlockAMR::setCopyLayerFine(BoundaryEdge edge, SWE_BlockGhost* layer) {
	int _nx = layer->nx;
	int _ny = layer->ny;
	int _rx = block_neighbour[edge]->getRefinementLevel() / getRefinementLevel();
	int _ry = block_neighbour[edge]->getRefinementLevel() / getRefinementLevel();
	int _nghosts = block_neighbour[edge]->getNghosts();

//	SWE_BlockGhost* tmp = proxyCopyLayer[edge]->refine_constant(_rx, _ry);

	SWE_BlockGhost* tmp;
	switch (edge) {
	case BND_LEFT: tmp = proxyCopyLayer[edge]->refine(_rx, _ry, dx, dy, scene, offsetX, offsetY+dy); break;
	case BND_RIGHT: tmp = proxyCopyLayer[edge]->refine(_rx, _ry, dx, dy, scene, offsetX+nx*dx, offsetY+dy); break;
	case BND_BOTTOM: tmp = proxyCopyLayer[edge]->refine(_rx, _ry, dx, dy, scene, offsetX+dx, offsetY); break;
	case BND_TOP: tmp = proxyCopyLayer[edge]->refine(_rx, _ry, dx, dy, scene, offsetX+dx, offsetY+ny*dy); break;
	}

	// copy the result in the copy layer
	// the result is always 2*_nghosts elements smaller in one direction
	switch (edge) {
	case BND_LEFT:
		for (int i=0; i<_nx; i++)
			for (int j=_nghosts; j<_ny-_nghosts; j++) {
				layer->h[i][j] = tmp->h[i][j-_nghosts];
				layer->hu[i][j] = tmp->hu[i][j-_nghosts];
				layer->hv[i][j] = tmp->hv[i][j-_nghosts];
			}
		break;
	case BND_RIGHT:
		for (int i=0; i<_nx; i++)
			for (int j=_nghosts; j<_ny-_nghosts; j++) {
				layer->h[i][j] = tmp->h[i][j-_nghosts];
				layer->hu[i][j] = tmp->hu[i][j-_nghosts];
				layer->hv[i][j] = tmp->hv[i][j-_nghosts];
			}
		break;
	case BND_BOTTOM:
		for (int i=_nghosts; i<_nx-_nghosts; i++)
			for (int j=0; j<_ny; j++) {
				layer->h[i][j] = tmp->h[i-_nghosts][j];
				layer->hu[i][j] = tmp->hu[i-_nghosts][j];
				layer->hv[i][j] = tmp->hv[i-_nghosts][j];
			}
		break;
	case BND_TOP:
		for (int i=_nghosts; i<_nx-_nghosts; i++)
			for (int j=0; j<_ny; j++) {
				layer->h[i][j] = tmp->h[i-_nghosts][j];
				layer->hu[i][j] = tmp->hu[i-_nghosts][j];
				layer->hv[i][j] = tmp->hv[i-_nghosts][j];
			}
		break;
	}
}

/**
 * @param edge
 * @return
 */
static BoundaryEdge SWE_BlockAMR::getOppositeEdge(BoundaryEdge edge) {
	switch (edge) {
	case BND_LEFT:
		return BND_RIGHT;
	case BND_RIGHT:
		return BND_LEFT;
	case BND_BOTTOM:
		return BND_TOP;
	case BND_TOP:
		return BND_BOTTOM;
	}
}

/**
 *
 */
void SWE_BlockAMR::increaseComputationalDomain() {
	if (boundary[BND_LEFT] == CONNECT) nxint_s--;
	if (boundary[BND_BOTTOM] == CONNECT) nyint_s--;
	if (boundary[BND_RIGHT] == CONNECT) nxint_e++;
	if (boundary[BND_TOP] == CONNECT) nyint_e++;
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
	nxint_s = nyint_s = 1;
	nxint_e = nx + 2 * nghosts - 2;
	nyint_e = ny + 2 * nghosts - 2;
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
 * @param edge	the boundary edge where the ghost layers must be set
 */
void SWE_BlockAMR::setGhostLayerEdge(BoundaryEdge edge) {

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

	switch (edge) {
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

