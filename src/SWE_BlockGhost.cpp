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

#include <stdlib.h>
#include "SWE_BlockGhost.hh"

using namespace std;

/**
 * Performs a coarsening of the block
 * by integrating on each coarse cell
 * @param rx	coarsening factor in the x-direction
 * @param ry	coarsening factor in the y-direction
 * @return		a new SWE_BlockGhost object containing the coarsened values
 */
SWE_BlockGhost* SWE_BlockGhost::coarsen(int rx, int ry) {
	int _nx = nx/rx;
	int _ny = ny/ry;
	Float2D* _h = new Float2D(_nx, _ny);
	Float2D* _b = new Float2D(_nx, _ny);
	Float2D* _hu = new Float2D(_nx, _ny);
	Float2D* _hv = new Float2D(_nx, _ny);

	// compute the coarse values as an average of the corresponding fine values
	for (int i=0; i<_nx; ++i)
		for (int j=0; j<_ny; ++j) {
			for (int k=i*rx; k<(i+1)*rx; ++k)
				for (int l=j*ry; l<(j+1)*ry; ++l) {
					(*_h)[i][j] += h[k][l];
					(*_hu)[i][j] += hu[k][l];
					(*_hv)[i][j] += hv[k][l];
				}
			(*_h)[i][j]  /= rx*ry;
			(*_hu)[i][j] /= rx*ry;
			(*_hv)[i][j] /= rx*ry;
		}

	return new SWE_BlockGhost(*_h, *_b, *_hu, *_hv);
}

inline float SWE_BlockGhost::sigma_x(const Float2D& Q, int i, int j, float dx) {
	return minmod(Q[i][j] - Q[i-1][j], Q[i+1][j] - Q[i][j])/dx;
}

inline float SWE_BlockGhost::sigma_y(const Float2D& Q, int i, int j, float dy) {
	return minmod(Q[i][j] - Q[i][j-1], Q[i][j+1] - Q[i][j])/dy;
}

inline float SWE_BlockGhost::sigma_x_2(const Float2D& Q1, const Float2D& Q2, int i, int j, float dx) {
	return minmod( 	(Q1[i][j]+Q2[i][j]) - (Q1[i-1][j]+Q2[i-1][j]),
				(Q1[i+1][j]+Q2[i+1][j]) - (Q1[i][j]+Q2[i][j]) )/dx;
}

inline float SWE_BlockGhost::sigma_y_2(const Float2D& Q1, const Float2D& Q2, int i, int j, float dy) {
	return minmod(	(Q1[i][j]+Q2[i][j]) - (Q1[i][j-1]+Q2[i][j-1]),
				(Q1[i][j+1]+Q2[i][j+1]) - (Q1[i][j]+Q2[i][j]) )/dy;
}

inline float SWE_BlockGhost::interpolate(float Q, float sigma_x, float sigma_y, float x, float xc, float y, float yc) {
	return Q + sigma_x*(x-xc) + sigma_y*(y-yc);
}

/**
 * Performs a higher-order reconstruction of the block
 * using bilinear interpolation
 * @param rx		refinement factor in the x-direction
 * @param ry		refinement factor in the y-direction
 * @param dx		mesh size of the current block in x-direction
 * @param dy		mesh size of the current block in y-direction
 * @param scene		pointer to scenario used for initial and boundary conditions
 * @param offsetX	offset of the current block in x-direction
 * @param offsetY	offset of the current block in y-direction
 * @return			a new SWE_BlockGhost object containing the refined values
 */
SWE_BlockGhost* SWE_BlockGhost::refine(int rx, int ry, float dx, float dy,
							SWE_Scenario *scene, float offsetX, float offsetY) {
	int _nx = (nx-2)*rx;
	int _ny = (ny-2)*ry;
	float _dx = dx/rx;
	float _dy = dy/ry;
	float 	sigma_h_x, sigma_h_y,
			sigma_b_x, sigma_b_y,
			sigma_hu_x, sigma_hu_y,
			sigma_hv_x, sigma_hv_y;
	float x, xc, y, yc;		// cell centres of fine and coarse cells
	float u_max, v_max;		// maximum fine velocities over one coarse cell
	float u_max_c, v_max_c;	// maximum coarse velocities over 3x3 neighbouring domain around a cell
	Float2D* _h = new Float2D(_nx, _ny);
	Float2D* _b = new Float2D(_nx, _ny);
	Float2D* _hu = new Float2D(_nx, _ny);
	Float2D* _hv = new Float2D(_nx, _ny);

	// bilinear interpolation
	for (int i=1; i<=nx-2; ++i)
		for (int j=1; j<=ny-2; ++j) {
			sigma_h_x = sigma_x_2(h, b, i, j, dx);
			sigma_h_y = sigma_y_2(h, b, i, j, dy);
			sigma_b_x = sigma_x(b, i, j, dx);
			sigma_b_y = sigma_y(b, i, j, dy);
			sigma_hu_x = sigma_x(hu, i, j, dx);
			sigma_hu_y = sigma_y(hu, i, j, dy);
			sigma_hv_x = sigma_x(hv, i, j, dx);
			sigma_hv_y = sigma_y(hv, i, j, dy);
			xc = offsetX + (i-0.5f)*dx;
			yc = offsetY + (j-0.5f)*dy;
			// maximum velocities in the fine cells inside coarse cell (i,j) in x and y direction
			u_max = 0.0;
			v_max = 0.0;
			for (int k=(i-1)*rx; k<i*rx; ++k)
				for (int l=(j-1)*ry; l<j*ry; ++l) {
					x = offsetX + (k+0.5f)*_dx;
					y = offsetY + (l+0.5f)*_dy;
					// TODO: interpolate bathymetry as well instead of getting it from the scenario
					// TODO: remove scene and other unnecessary parameters
					// (*_b)[k][l] = scene->getBathymetry(x, y);
					(*_b)[k][l] = interpolate(b[i][j], sigma_b_x, sigma_b_y, x, xc, y, yc);
					(*_h)[k][l] = interpolate(h[i][j]+b[i][j], sigma_h_x, sigma_h_y, x, xc, y, yc) - (*_b)[k][l];
					(*_hu)[k][l] = interpolate(hu[i][j], sigma_hu_x, sigma_hu_y, x, xc, y, yc);
					(*_hv)[k][l] = interpolate(hv[i][j], sigma_hv_x, sigma_hv_y, x, xc, y, yc);
					if ((*_hu)[k][l] / (*_h)[k][l] > u_max)
						u_max = (*_hu)[k][l] / (*_h)[k][l];
					if ((*_hv)[k][l] / (*_h)[k][l] > v_max)
						v_max = (*_hv)[k][l] / (*_h)[k][l];
				}

			// maximum velocities in x and y directions in the coarse cell (i,j) and neighbours
			// we consider all nine coarse cells including and surrounding cell (i,j)
			u_max_c = 0.0;
			v_max_c = 0.0;
			for (int k=i-1; k<=i+1; ++k) {
				for (int l=j-1; l<=j+1; ++l) {
					if (hu[k][l] / h[k][l] > u_max_c)
						u_max_c = hu[k][l] / h[k][l];
					if (hv[k][l] / h[k][l] > v_max_c)
						v_max_c = hv[k][l] / h[k][l];
				}
			}
			// if new velocity extrema are introduced, limit all the velocities in the fine cells
			if (u_max > u_max_c)
				for (int k=(i-1)*rx; k<i*rx; ++k)
					for (int l=(j-1)*ry; l<j*ry; ++l)
						(*_hu)[k][l] = (*_h)[k][l] * hu[i][j]/h[i][j];
			if (v_max > v_max_c)
				for (int k=(i-1)*rx; k<i*rx; ++k)
					for (int l=(j-1)*ry; l<j*ry; ++l)
						(*_hv)[k][l] = (*_h)[k][l] * hv[i][j]/h[i][j];
		}
	return new SWE_BlockGhost(*_h, *_b, *_hu, *_hv);
}

/**
 * Performs a higher-order reconstruction of the block
 * using a piecewise constant interpolation
 * @param rx	refinement factor in the x-direction
 * @param ry	refinement factor in the y-direction
 * @return		a new SWE_BlockGhost object containing the refined values
 */
SWE_BlockGhost* SWE_BlockGhost::refine_constant(int rx, int ry) {
	int _nx = (nx-2)*rx;
	int _ny = (ny-2)*ry;
	Float2D* _h = new Float2D(_nx, _ny);
	Float2D* _b = new Float2D(_nx, _ny);
	Float2D* _hu = new Float2D(_nx, _ny);
	Float2D* _hv = new Float2D(_nx, _ny);

	// piecewise constant
	for (int i=1; i<=nx-2; ++i)
		for (int j=1; j<=ny-2; ++j)
			for (int k=(i-1)*rx; k<i*rx; ++k)
				for (int l=(j-1)*ry; l<j*ry; ++l) {
					(*_h)[k][l] = h[i][j];
					(*_hu)[k][l] = hu[i][j];
					(*_hv)[k][l] = hv[i][j];
				}
	return new SWE_BlockGhost(*_h, *_b, *_hu, *_hv);
}

/**
 * Method to retrieve a proxy to a certain column
 * @param c column number
 * @return 	the requested column
 */
SWE_Block1D* SWE_BlockGhost::getColProxy(int c) {
	return new SWE_Block1D( h.getColProxy(c), hu.getColProxy(c), hv.getColProxy(c));
}

/**
 * Method to retrieve a proxy to a certain row
 * @param r row number
 * @return	the requested row
 */
SWE_Block1D* SWE_BlockGhost::getRowProxy(int r) {
	return new SWE_Block1D( h.getRowProxy(r), hu.getRowProxy(r), hv.getRowProxy(r));
}

/**
 * Method to retrieve a proxy to a certain 2D sub-block
 * @param c	 starting index on the x-axis
 * @param r	 starting index on the y-axis
 * @param nc number of columns in the sub-block
 * @param nr number of rows in the sub-block
 * @return	 the requested sub-block
 */
SWE_BlockGhost* SWE_BlockGhost::getBlockProxy(int c, int r, int nc, int nr) {
	return new SWE_BlockGhost(	*h.getBlockProxy(c, r, nc, nr),
								*b.getBlockProxy(c, r, nc, nr),
								*hu.getBlockProxy(c, r, nc, nr),
								*hv.getBlockProxy(c, r, nc, nr)
							  );
}
