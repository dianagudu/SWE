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

#ifndef _SWE_BlockAMR_HH_
#define _SWE_BlockAMR_HH_

#include <iostream>
#include <stdio.h>
#include <fstream>
#include "tools/help.hh"
#include "SWE_Block.hh"
#include "SWE_BlockGhost.hh"

using namespace std;

/**
 * SWE_BlockAMR extends the base class SWE_Block, and provides
 * an implementation of a simple shallow water model augmented
 * with a block-adaptive mesh refinement
 */
class SWE_BlockAMR : public SWE_Block {
  public:
	// constructor
	SWE_BlockAMR(float _offsetX, float _offsetY, int _nx, int _ny, float _dx, float _dy, int _rx=1, int _ry=1,
					 int _nghosts=1, InterpolationType _interpolationScheme = APPROX_TIME_SPACE );

	// object methods
	virtual SWE_BlockGhost* registerCopyLayer(BoundaryEdge edge);
    virtual void synchCopyLayerBeforeRead(TimeSteppingType timeStepping, BoundaryEdge edge, float t, float dt_c);
    virtual void synchAfterWrite();
    virtual void synchBeforeRead();
    void initScenario(SWE_Scenario &i_scenario, const bool i_multipleBlocks);

    /**
     * Gets the refinement level of this block
     * @return the refinement factor
     */
    int getRefinementLevel() {return rx;};

    void decreaseComputationalDomain();
    void resetComputationalDomainMax();
    void resetComputationalDomainMin();

    void setGhostLayerEdge(BoundaryEdge edge);
    void setGhostLayerEdge(BoundaryEdge i_edge, SWE_BlockGhost* i_neighbour);

    void setNeighbourRefinementLevel(int i_rxy, BoundaryEdge edge);

#ifdef NOMPI
    void setBlockNeighbour(SWE_BlockAMR* neigh, BoundaryEdge edge);

    /**
     * Method that gets the block neighbour at a given boundary edge
     *
     * @param edge	one of the four boundaries of a block
     * @return		a pointer to a neighbouring block
     */
    SWE_BlockAMR* getNeighbour(BoundaryEdge edge) { return block_neighbour[edge]; };
#endif

  private:
    SWE_BlockGhost* getProxyCopyLayer(BoundaryEdge edge);
    // methods for single ghost layer
    SWE_BlockGhost* getFineProxyCopyLayer_singleLayer(BoundaryEdge edge);
    SWE_BlockGhost* getCoarseProxyCopyLayer_singleLayer(BoundaryEdge edge);
    // methods for multi ghost layer
    SWE_BlockGhost* getFineProxyCopyLayer_multiLayer(BoundaryEdge edge);
    SWE_BlockGhost* getCoarseProxyCopyLayer_multiLayer(BoundaryEdge edge);

    void setCopyLayerCoarse(BoundaryEdge edge);
    void setCopyLayerFine(BoundaryEdge edge, SWE_BlockGhost* layer);
    void setCopyLayerFineConstant(BoundaryEdge edge);
    void timeInterpolateCopyLayer(BoundaryEdge edge, float t, float dt_c);

  protected:
    /// refinement factor in the x-direction
    int rx;
    /// refinement factor in the y-direction
    int ry;
    /// copy layers for all edges
    SWE_BlockGhost* copyLayer[4];
    /// proxy copy layers for all edges
    SWE_BlockGhost* proxyCopyLayer[4];

    // only used for (approximate) time space interpolation
    /// copy layers at beginning of time-step
    SWE_BlockGhost* startCopyLayer[4];
    /// copy layers at the end of time-step
    SWE_BlockGhost* endCopyLayer[4];
    /// changes in copy layers in one time-step
    SWE_BlockGhost* delta[4];

    /// refinement levels of neighbours
    int neighbourRefinementLevel[4];
#ifdef NOMPI
    /// pointers to neighbouring blocks in each direction (only used in serial version)
    SWE_BlockAMR* block_neighbour[4];
#endif

    /// scenario
    SWE_Scenario* scene;
    /// type of interpolation used for the ghost layers; determines the number of ghost layers
    InterpolationType interpolationStrategy;
};

#endif /* _SWE_BlockAMR_HH_ */
