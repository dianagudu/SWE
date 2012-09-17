/**
 * @file
 * This file is part of SWE.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
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
 * SWE_BlockAMR, which uses solvers in the wave propagation formulation.
 */

#include <cassert>
#include <string>
#include <limits>
#include "SWE_BlockAMR.hh"
#include "SWE_WavePropagationAMR.hh"
#ifdef LOOP_OPENMP
#include <omp.h>
#endif

/**
 * Constructor of a SWE_WavePropagationAMR.
 *
 * Allocates the variables for the simulation:
 *   unknowns h,hu,hv,b are defined on grid indices [0,..,nx+2*nghosts-1]*[0,..,ny+2*nghosts-1]
 *   (-> Abstract class SWE_BlockAMR)
 *     -> computational domain is [nxint_s,..,nxint_e]*[nyint_s,..,nyint_e]
 *     -> some ghost cell layers can also be part of the computational domain
 *
 *   net-updates are defined for edges with indices [nxint_s-1,..,nxint_e]*[nyint_s-1,..,nyint_e-1]
 *   or [nxint_s-1,..,nxint_e-1]*[nyint_s-1,..,nyint_e] (for horizontal/vertical edges)
 *
 *   A left/right net update with index (i-1,j-1) is located on the edge between
 *   cells with index (i-1,j) and (i,j):
 *
 *   *********************
 *   *         *         *
 *   * (i-1,j) *  (i,j)  *
 *   *         *         *
 *   *********************
 *
 *             *
 *            ***
 *           *****
 *             *
 *             *
 *   NetUpdatesLeft(i-1,j-1)
 *             or
 *   NetUpdatesRight(i-1,j-1)
 *
 *
 *   A below/above net update with index (i-1, j-1) is located on the edge between
 *   cells with index (i, j-1) and (i,j):
 *
 *   ***********
 *   *         *
 *   * (i, j)  *   *
 *   *         *  **  NetUpdatesBelow(i-1,j-1)
 *   *********** *****         or
 *   *         *  **  NetUpdatesAbove(i-1,j-1)
 *   * (i,j-1) *   *
 *   *         *
 *   ***********
 *
 */
SWE_WavePropagationAMR::SWE_WavePropagationAMR( float i_offsetX,
                                                    	float i_offsetY,
                                                    	int i_nx,
                                                    	int i_ny,
                                                    	float i_dx,
                                                    	float i_dy,
                                                    	int i_rx,
                                                    	int i_ry,
                                                    	int i_nghosts,
                                                    	InterpolationType i_interpolationStrategy):
  SWE_BlockAMR(i_offsetX,i_offsetY,i_nx,i_ny,i_dx,i_dy,i_rx,i_ry,i_nghosts,i_interpolationStrategy),
  wavePropagationSolver(SIMULATION_TSUNAMI_ZERO_THRESHOLD, 1.0),
  hNetUpdatesLeft  (i_nx+2*i_nghosts-1, i_ny+2*i_nghosts-2),
  hNetUpdatesRight (i_nx+2*i_nghosts-1, i_ny+2*i_nghosts-2),
  huNetUpdatesLeft (i_nx+2*i_nghosts-1, i_ny+2*i_nghosts-2),
  huNetUpdatesRight(i_nx+2*i_nghosts-1, i_ny+2*i_nghosts-2),
  #if WAVE_PROPAGATION_SOLVER==3
  hvNetUpdatesLeft (i_nx+2*i_nghosts-1, i_ny+2*i_nghosts-2),
  hvNetUpdatesRight(i_nx+2*i_nghosts-1, i_ny+2*i_nghosts-2),
  #endif

  hNetUpdatesBelow (i_nx+2*i_nghosts-2, i_ny+2*i_nghosts-1),
  hNetUpdatesAbove (i_nx+2*i_nghosts-2, i_ny+2*i_nghosts-1),
  #if WAVE_PROPAGATION_SOLVER==3
  huNetUpdatesBelow(i_nx+2*i_nghosts-2, i_ny+2*i_nghosts-1),
  huNetUpdatesAbove(i_nx+2*i_nghosts-2, i_ny+2*i_nghosts-1),
  #endif
  hvNetUpdatesBelow(i_nx+2*i_nghosts-2, i_ny+2*i_nghosts-1),
  hvNetUpdatesAbove(i_nx+2*i_nghosts-2, i_ny+2*i_nghosts-1)
{}

/**
 * Compute net updates for the block.
 * The member variable #maxTimestep will be updated with the 
 * maximum allowed time step size
 */
void SWE_WavePropagationAMR::computeNumericalFluxes() {

#ifdef DBG
	cout << "Compute numerical fluxes " << endl << flush;
#endif //DBG

  //maximum (linearized) wave speed within one iteration
  float maxWaveSpeed = (float) 0.;

  //compute the net-updates for the vertical edges
  #ifdef LOOP_OPENMP
  #pragma omp parallel
  {
  float l_maxWaveSpeed = (float) 0.;
  solver::Hybrid<float> wavePropagationSolver;
  #pragma omp for
  #endif
  for(int i = nxint_s; i <= nxint_e + 1; i++) {
    for(int j = nyint_s; j <= nyint_e; j++) {
      float maxEdgeSpeed;
      #if WAVE_PROPAGATION_SOLVER!=3
      wavePropagationSolver.computeNetUpdates( h[i-1][j], h[i][j],
                                               hu[i-1][j], hu[i][j],
                                               b[i-1][j], b[i][j],
                                               hNetUpdatesLeft[i-1][j-1], hNetUpdatesRight[i-1][j-1],
                                               huNetUpdatesLeft[i-1][j-1], huNetUpdatesRight[i-1][j-1],
                                               maxEdgeSpeed );
     #else
     //TODO: implement again.
     assert(false);
     #endif

      #ifdef LOOP_OPENMP
      //update the maximum wave speed
      l_maxWaveSpeed = std::max(l_maxWaveSpeed, maxEdgeSpeed);
      #else
      //update the maximum wave speed
      maxWaveSpeed = std::max(maxWaveSpeed, maxEdgeSpeed);
      #endif
    }
  }

  //compute the net-updates for the horizontal edges
  #ifdef LOOP_OPENMP
  #pragma omp for
  #endif
  for(int i = nxint_s; i <= nxint_e; i++) {
    for(int j = nyint_s; j <= nyint_e + 1; j++) {
      float maxEdgeSpeed;
      #if WAVE_PROPAGATION_SOLVER!=3
      wavePropagationSolver.computeNetUpdates( h[i][j-1], h[i][j],
                                               hv[i][j-1], hv[i][j],
                                               b[i][j-1], b[i][j],
                                               hNetUpdatesBelow[i-1][j-1], hNetUpdatesAbove[i-1][j-1],
                                               hvNetUpdatesBelow[i-1][j-1], hvNetUpdatesAbove[i-1][j-1],
                                               maxEdgeSpeed );
      #else
      //TODO: implement again.
      assert(false);
      #endif

      #ifdef LOOP_OPENMP
      //update the maximum wave speed
      l_maxWaveSpeed = std::max(l_maxWaveSpeed, maxEdgeSpeed);
      #else
      //update the maximum wave speed
      maxWaveSpeed = std::max(maxWaveSpeed, maxEdgeSpeed);
      #endif
    }
  }
  #ifdef LOOP_OPENMP
  #pragma omp critical
  {
  maxWaveSpeed = std::max(l_maxWaveSpeed, maxWaveSpeed);
  }
  } // end of parallel for block
  #endif

  if(maxWaveSpeed > SIMULATION_TSUNAMI_ZERO_THRESHOLD) { //TODO zeroTol
    //compute the time step width
    //CFL-Codition
    //(max. wave speed) * dt / dx < .5
    // => dt = .5 * dx/(max wave speed)
    maxTimestep = std::min( dx/maxWaveSpeed, dy/maxWaveSpeed );

//    #if WAVE_PROPAGATION_SOLVER!=3
    maxTimestep *= (float) .4; //CFL-number = .5
//    #else
//    dt *= (float) .8; //CFL-number = 1. (wave limiters)
//    #endif
  }
  else
    maxTimestep = std::numeric_limits<float>::max(); //might happen in dry cells

}

/**
 * Updates the unknowns with the already computed net-updates.
 *
 * @param dt time step width used in the update.
 */
void SWE_WavePropagationAMR::updateUnknowns(float dt) {
#ifdef DBG
	cout << "Update unknowns " << endl << flush;
#endif //DBG
  //update cell averages with the net-updates
  #ifdef LOOP_OPENMP
  #pragma omp parallel for
  #endif
  for(int i = nxint_s; i <= nxint_e; i++) {
      for(int j = nyint_s; j <= nyint_e; j++) {
        h[i][j] -=   dt/dx * (hNetUpdatesRight[i-1][j-1] + hNetUpdatesLeft[i][j-1])
                   + dt/dy * (hNetUpdatesAbove[i-1][j-1] + hNetUpdatesBelow[i-1][j]);
        hu[i][j] -= dt/dx * (huNetUpdatesRight[i-1][j-1] + huNetUpdatesLeft[i][j-1]);
        hv[i][j] -= dt/dy * (hvNetUpdatesAbove[i-1][j-1] + hvNetUpdatesBelow[i-1][j]);
        #if WAVE_PROPAGATION_SOLVER==3
        hv[i][j] -= dt/dx * (hvNetUpdatesRight[i-1][j-1] + hvNetUpdatesLeft[i][j-1]);
        hu[i][j] -= dt/dy * (huNetUpdatesAbove[i-1][j-1] + huNetUpdatesBelow[i-1][j]);
        #endif

        //TODO: dryTol
        if(h[i][j] < 0) {
          if(h[i][j] < -SIMULATION_TSUNAMI_ZERO_THRESHOLD) {
            std::cerr << "Warning, negative height: (i,j)=(" << i << "," << j << ")=" << h[i][j] << std::endl;
            std::cerr << "         b: " << b[i][j] << std::endl;
          }
          //zero (small) negative depths
          h[i][j] = hu[i][j] = hv[i][j] = 0.;
        }
        else if(h[i][j] < SIMULATION_TSUNAMI_ZERO_THRESHOLD)
          hu[i][j] = hv[i][j] = 0.; //no water, no speed!
      }
  }
}

/**
 * Update the bathymetry values with the displacement corresponding to the current time step.
 *
 * @param i_asagiScenario the corresponding ASAGI-scenario
 */
#ifdef DYNAMIC_DISPLACEMENTS
bool SWE_WavePropagationAMR::updateBathymetryWithDynamicDisplacement(scenarios::Asagi &i_asagiScenario, const float i_time) {
  if (!i_asagiScenario.dynamicDisplacementAvailable(i_time))
    return false;

  // update the bathymetry
  for(int i=0; i<nx+2*nghosts; i++) {
    for(int j=0; j<ny+2*nghosts; j++) {
      b[i][j] = i_asagiScenario.getBathymetryAndDynamicDisplacement( offsetX + (i-0.5f)*dx,
                                                                     offsetY + (j-0.5f)*dy,
                                                                     i_time
                                                                   );
    }
  }
  return true;
}
#endif

/**
 * Executes a single timestep.
 *  * compute net updates for every edge
 *  * update cell values with the net updates
 *
 * @param dt	time step width of the update
 */
void SWE_WavePropagationAMR::simulateTimestep(float dt) {
  computeNumericalFluxes();
  updateUnknowns(dt);
}

/**
 * Runs the simulation until i_tEnd is reached.
 *
 * @param i_tStart time when the simulation starts
 * @param i_tEnd  time when the simulation should end
 * @return time we reached after the last update step, in general a bit later than i_tEnd
 */
float SWE_WavePropagationAMR::simulate(float i_tStart,float i_tEnd) {
  float t = i_tStart;
  do {
     //set values in ghost cells
     setGhostLayer();

     // compute net updates for every edge
     computeNumericalFluxes();
     //execute a wave propagation time step
     updateUnknowns(maxTimestep);
     t += maxTimestep;

     std::cout << "Simulation at time " << t << std::endl << std::flush;
  } while(t < i_tEnd);

  return t;
}

#if WAVE_PROPAGATION_SOLVER==0
#ifndef NDEBUG
#ifndef LOOP_OPENMP
/**
 * Gets the current solver statistics for a hybrid solver.
 *
 * @param o_counterFWave will be set to: times the f-Wave solver was used.
 * @param o_counterAugRie will be set to: times the Augmented Riemann solver was used.
 */
void SWE_WavePropagationAMR::getSolverStats( long &o_counterFWave, long &o_counterAugRie ) {
  wavePropagationSolver.getStats( o_counterFWave, o_counterAugRie );
}
#endif
#endif
#endif
