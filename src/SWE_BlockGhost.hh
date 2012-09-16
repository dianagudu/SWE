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

#ifndef _SWE_BlockGhost_HH_
#define _SWE_BlockGhost_HH_

#include <iostream>
#include "tools/help.hh"
#include "scenarios/SWE_Scenario.h"

/**
 * enum type: available types of time-stepping strategies used for adaptive refinement:
 * global/local time-stepping
 */
typedef enum TimeSteppingType {
	GTS, LTS
} TimeSteppingType;

/**
 * enum type: available interpolation schemes for ghost cell exchange in LTS
 * this type determines the number of ghost layers in a block
 */
typedef enum InterpolationType {
	APPROX_TIME_SPACE, TIME_SPACE, SPACE
} InterpolationType;

/**
 * SWE_Block1D is a simple struct that can represent a single line or row of
 * SWE_Block unknowns (using the Float1D proxy class).
 * It is intended to unify the implementation of inflow and periodic boundary
 * conditions, as well as the ghost/copy-layer connection between several SWE_Block
 * grids.
 */
struct SWE_Block1D {
    SWE_Block1D(const Float1D& _h, const Float1D& _hu, const Float1D& _hv)
    : h(_h), hu(_hu), hv(_hv) {};
    SWE_Block1D(float* _h, float* _hu, float* _hv, int _size, int _stride=1)
    : h(_h,_size,_stride), hu(_hu,_size,_stride), hv(_hv,_size,_stride) {};
    SWE_Block1D(const int len)
    : h(len), hu(len), hv(len) {};

    Float1D h;
    Float1D hu;
    Float1D hv;
};

/**
 * SWE_BlockGhost is a class that models 2D ghost/copy layers containing all the unknowns
 * It contains methods for coarsening and refining.
 */
class SWE_BlockGhost {
public:
	// constructing an object given the arrays
    SWE_BlockGhost(const Float2D& _h, const Float2D& _b, const Float2D& _hu, const Float2D& _hv)
    : nx(_h.getCols()), ny(_h.getRows()), h(_h), b(_b), hu(_hu), hv(_hv) {};

    // constructing an empty object given the dimensions
    SWE_BlockGhost(int _nx, int _ny)
    : nx(_nx), ny(_ny), h(_nx, _ny), b(_nx, _ny), hu(_nx, _ny), hv(_nx, _ny) {};

    // Coarsens whole the block by a factor of rx in x-direction and ry in y-direction
    SWE_BlockGhost* coarsen(int rx, int ry);

    // Refines the inner block by a factor of rx in x-direction and ry in y-direction
    // The outer block consists of one cell at each boundary and is needed for the interpolation scheme
    SWE_BlockGhost* refine(int rx, int ry, float dx, float dy, SWE_Scenario *scene, float offsetX, float offsetY);

    // Refines the inner block by a factor of rx in x-direction and ry in y-direction
    // The outer block consists of one cell at each boundary and is only present for consistency with the
    // other refining strategy
    SWE_BlockGhost* refine_constant(int rx, int ry);

    // methods to access only a certain part of the current block
    // (either a column, row or 2D block)
    // using a proxy class
    SWE_Block1D* getColProxy(int c);
    SWE_Block1D* getRowProxy(int r);
    SWE_BlockGhost* getBlockProxy(int c, int r, int nc, int nr);

private:
    // helper methods for the bilinear interpolation
    inline float sigma_x(const Float2D& v, int i, int j, float dx);
	inline float sigma_y(const Float2D& v, int i, int j, float dy);
	inline float sigma_x_2(const Float2D& v1, const Float2D& v2, int i, int j, float dx);
	inline float sigma_y_2(const Float2D& v1, const Float2D& v2, int i, int j, float dy);
	inline float interpolate(float v, float sigma_x, float sigma_y, float x, float xc, float y, float yc);

public:
	// size of the current block in x- and y-direction
    int nx;
    int ny;
    // 2D arrays to store the unknowns and the bathymetry
    Float2D h;
    Float2D b;
    Float2D hu;
    Float2D hv;
};

#endif //_SWE_BlockGhost_HH_
