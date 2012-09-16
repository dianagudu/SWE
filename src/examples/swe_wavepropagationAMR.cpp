/**
 * @file
 * This file is part of SWE.
 *
 * @author Alexander Breuer (breuera AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *         Michael Bader (bader AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Univ.-Prof._Dr._Michael_Bader)
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
 * Very basic setting of SWE, which uses a wave propagation solver with
 * an adaptively refined grid and an artificial scenario on a single block.
 * There is also the possibility of using the ASAGI library for setting the scenario.
 */

#include "../tools/help.hh"
#include <cstdlib>
#include <string>

#include "../SWE_BlockAMR.hh"
#include "../SWE_BlockGhost.hh"

#ifndef CUDA
#include "../SWE_WavePropagationAMR.hh"
#else
#include "../SWE_WavePropagationBlockCuda.hh"
#endif

#include "../tools/Logger.hpp"

#ifdef WRITENETCDF
#include "../tools/NetCdfWriter.hh"
#endif

#ifdef ASAGI
#include "../scenarios/SWE_AsagiScenario.hpp"
#else
#include "../scenarios/SWE_simple_scenarios.h"
#include "../scenarios/SWE_SingleWaveOnSimpleBeach.h"
#endif

#ifndef STATICLOGGER
#define STATICLOGGER
#include "../tools/Logger.hpp"
static tools::Logger s_sweLogger;
#endif

#ifdef BENCHMARKING
#include "../benchmarking/BenchmarkingDataReceiver.hpp"
#endif

#include "../SWE_BlockManager.hh"

/**
 * Main program for the simulation on a single SWE_WavePropagationBlock.
 */
int main( int argc, char** argv ) {
  /**
   * Initialization.
   */
  // check if the necessary command line input parameters are given
	if(argc != 9) {
		std::cout << "Aborting ... please provide proper input parameters." << std::endl
				  << "Example: ./SWE_parallel 10 25 200 300 /work/openmp_out /work/config <time stepping> <interpolation scheme>" << std::endl
				  << "\tfor 2000 * 7500 cells with 10 * 25 blocks of size 200 * 300" << std::endl
				  << "\ttime stepping can be GTS / LTS for global/local time stepping" << std::endl
				  << "\tinterpolation scheme can be 0 / 1 / 2 for approximate time-space / time-space / space interpolation" << std::endl;
		return 1;
	}

  //! number of grid cells in x- and y-direction.
  int l_nX, l_nY;

  //! number of patches in x- and y-direction
  int l_blockX, l_blockY;

  //! l_baseName of the plots.
  std::string l_baseName;

  //! time stepping strategy: global or local
  TimeSteppingType l_timeSteppingStrategy;

  //! interpolation scheme for the ghost layers
  InterpolationType l_interpolationScheme;

  //! number of ghost layers
  int l_nghosts;

  // read command line parameters
  l_blockX = atoi(argv[1]);
  l_blockY = atoi(argv[2]);

  l_nY = l_nX = atoi(argv[3]);
  l_nY = atoi(argv[4]);
  l_baseName = std::string(argv[5]);

  // read patch configuration with refinement levels
  std::string l_configName = std::string(argv[6]);
  std::ifstream configFile(argv[6]);

  if(!configFile) {
    cout<<endl<<"Failed to open file"<<l_configName;
    return 1;
  }

  int** l_rxy = new int*[l_blockX];
  for (int i=0; i<l_blockX; ++i)
	  l_rxy[i] = new int[l_blockY];
  for (int j=l_blockY-1; j>=0; --j)
	  for (int i=0; i<l_blockX; ++i)
		  configFile>>l_rxy[i][j];

  configFile.close();


  // read time-stepping strategy
  std::string l_ts = std::string(argv[7]);
  if (l_ts.compare("GTS") == 0)
	  l_timeSteppingStrategy = GTS;
  else if (l_ts.compare("LTS") == 0)
	  l_timeSteppingStrategy = LTS;
  else {
	  cout<<endl<<"Wrong time-stepping strategy"<<l_ts<<endl;
	  return 1;
  }

  // read interpolation scheme
  int l_is = atoi(argv[8]);
  switch (l_is) {
  case 0: l_interpolationScheme = APPROX_TIME_SPACE; break;
  case 1: l_interpolationScheme = TIME_SPACE; break;
  case 2: l_interpolationScheme = SPACE; break;
  default: cout<<endl<<"Wrong interpolation scheme"<<l_is<<endl; return 1;
  }


#ifdef ASAGI
  /* Information about the example bathymetry grid (tohoku_gebco_ucsb3_500m_hawaii_bath.nc):
     *
     * Pixel node registration used [Cartesian grid]
     * Grid file format: nf = GMT netCDF format (float)  (COARDS-compliant)
     * x_min: -500000 x_max: 6500000 x_inc: 500 name: x nx: 14000
     * y_min: -2500000 y_max: 1500000 y_inc: 500 name: y ny: 8000
     * z_min: -6.48760175705 z_max: 16.1780223846 name: z
     * scale_factor: 1 add_offset: 0
     * mean: 0.00217145586762 stdev: 0.245563641735 rms: 0.245573241263
     */

  //simulation area
  float simulationArea[4];
  simulationArea[0] = -450000;
  simulationArea[1] = 3450000;
  simulationArea[2] = -2450000;
  simulationArea[3] = 1450000;

  SWE_AsagiScenario l_scenario( "/home/diana/workspace_c++/geo_information/tohoku_gebco_ucsb3_500m_hawaii_bath.nc",
                            	"/home/diana/workspace_c++/geo_information/tohoku_gebco_ucsb3_500m_hawaii_displ.nc",
                            	(float) 28800., simulationArea);
#else
  //create a simple artificial scenario
//  SWE_BathymetryDamBreakScenario l_scenario;
//  SWE_SplashingPoolScenario l_scenario;
//  SWE_SplashingConeScenario l_scenario;
//  SWE_SeaAtRestScenario l_scenario;
  SWE_SingleWaveOnSimpleBeach l_scenario;
#endif
  //! number of checkpoints for visualization (at each checkpoint in time, an output file is written).
  int l_numberOfCheckPoints = 40;

  //! size of a single cell in x- and y-direction
  float l_dX, l_dY;

  // compute the size of a single cell
  l_dX = (l_scenario.getBoundaryPos(BND_RIGHT) - l_scenario.getBoundaryPos(BND_LEFT) )/(l_blockX*l_nX);
  l_dY = (l_scenario.getBoundaryPos(BND_TOP) - l_scenario.getBoundaryPos(BND_BOTTOM) )/(l_blockY*l_nY);

  //! origin of the simulation domain in x- and y-direction
  float l_originX = l_scenario.getBoundaryPos(BND_LEFT);
  float l_originY = l_scenario.getBoundaryPos(BND_BOTTOM);

#ifndef CUDA
  // define array of Cartesian grid blocks+1
  // and initialise it to given scenario
  SWE_WavePropagationAMR*** l_wavePropagationBlock = new SWE_WavePropagationAMR**[l_blockX];
  for(int i=0; i<l_blockX; ++i)
    l_wavePropagationBlock[i] = new SWE_WavePropagationAMR*[l_blockY];

  for(int i=0; i<l_blockX; ++i) {
     for(int j=0; j<l_blockY; ++j) {

    	 if (l_interpolationScheme == APPROX_TIME_SPACE)
    		l_nghosts = 1;
    	else
    		l_nghosts = 2*l_rxy[i][j];

        l_wavePropagationBlock[i][j] = new SWE_WavePropagationAMR(
        		l_originX + i*l_nX*l_dX,
        		l_originY + j*l_nY*l_dY,
        		l_nX*l_rxy[i][j], l_nY*l_rxy[i][j],
        		l_dX/l_rxy[i][j], l_dY/l_rxy[i][j],
        		l_rxy[i][j], l_rxy[i][j],
        		l_nghosts, // size of ghost layer
        		l_interpolationScheme); // interpolation scheme for ghost layers
        l_wavePropagationBlock[i][j]->initScenario(l_scenario, true);
     };
  };

  // set outflow boundaries
  for (i=0; i<l_blockX; i++) {
	  l_wavePropagationBlock[i][0]->setBoundaryType(BND_BOTTOM, OUTFLOW, NULL);
	  l_wavePropagationBlock[i][l_blockY-1]->setBoundaryType(BND_TOP, OUTFLOW, NULL);
  }
  for (i=0; i<l_blockY; i++) {
  	  l_wavePropagationBlock[0][i]->setBoundaryType(BND_LEFT, OUTFLOW, NULL);
  	  l_wavePropagationBlock[l_blockX-1][i]->setBoundaryType(BND_RIGHT, OUTFLOW, NULL);
  }

  // connect SWE blocks at boundaries
  // first left and right boundaries
  SWE_BlockGhost*** leftOutflow  = new SWE_BlockGhost**[l_blockX-1];
  SWE_BlockGhost*** rightOutflow = new SWE_BlockGhost**[l_blockX-1];
  for(int i=0; i<l_blockX-1; ++i) {
     leftOutflow[i]  = new SWE_BlockGhost*[l_blockY];
     rightOutflow[i] = new SWE_BlockGhost*[l_blockY];
     for(int j=0; j<l_blockY; ++j) {
		l_wavePropagationBlock[i+1][j]->setBlockNeighbour(l_wavePropagationBlock[i][j], BND_LEFT);
		l_wavePropagationBlock[i][j]->setBlockNeighbour(l_wavePropagationBlock[i+1][j], BND_RIGHT);
        leftOutflow[i][j] = l_wavePropagationBlock[i+1][j]->registerCopyLayer(BND_LEFT);
        l_wavePropagationBlock[i][j]->setBoundaryType(BND_RIGHT, CONNECT, leftOutflow[i][j]);
        rightOutflow[i][j] = l_wavePropagationBlock[i][j]->registerCopyLayer(BND_RIGHT);
        l_wavePropagationBlock[i+1][j]->setBoundaryType(BND_LEFT, CONNECT, rightOutflow[i][j]);
     };
  };

  // then bottom and top boundaries
  SWE_BlockGhost*** botOutflow = new SWE_BlockGhost**[l_blockX];
  SWE_BlockGhost*** topOutflow = new SWE_BlockGhost**[l_blockX];
  for(int i=0; i<l_blockX; ++i) {
     botOutflow[i] = new SWE_BlockGhost*[l_blockY-1];
     topOutflow[i] = new SWE_BlockGhost*[l_blockY-1];
     for(int j=0; j<l_blockY-1; ++j) {
		l_wavePropagationBlock[i][j+1]->setBlockNeighbour(l_wavePropagationBlock[i][j], BND_BOTTOM);
		l_wavePropagationBlock[i][j]->setBlockNeighbour(l_wavePropagationBlock[i][j+1], BND_TOP);
        botOutflow[i][j] = l_wavePropagationBlock[i][j+1]->registerCopyLayer(BND_BOTTOM);
        l_wavePropagationBlock[i][j]->setBoundaryType(BND_TOP, CONNECT, botOutflow[i][j]);
        topOutflow[i][j] = l_wavePropagationBlock[i][j]->registerCopyLayer(BND_TOP);
        l_wavePropagationBlock[i][j+1]->setBoundaryType(BND_BOTTOM, CONNECT, topOutflow[i][j]);
     };
  };

  #else
  // get the origin from the scenario
  l_originX = l_scenario.getBoundaryPos(BND_LEFT);
  l_originY = l_scenario.getBoundaryPos(BND_BOTTOM);
  // create a single wave propagation block
  SWE_WavePropagationBlockCuda l_wavePropgationBlock(l_originX, l_originY);
  // initialize the wave propgation block
  l_wavePropgationBlock.initScenario(l_scenario);
  #endif

  //! time when the simulation ends.
  float l_endSimulation = l_scenario.endSimulation();

  //! checkpoints when output files are written.
  float* l_checkPoints = new float[l_numberOfCheckPoints+1];

  // compute the checkpoints in time
  for(int cp = 0; cp <= l_numberOfCheckPoints; ++cp) {
     l_checkPoints[cp] = cp*(l_endSimulation/l_numberOfCheckPoints);
  }

#ifdef WRITENETCDF
  for(int i=0; i<l_blockX; ++i) {
     for(int j=0; j<l_blockY; ++j) {
       //construct a NetCdfWriter
       std::string fileName = generateFileName(l_baseName,0,i,j, ".nc");
       s_sweLogger.printOutputFileCreation(fileName,i,j);

       io::NetCdfWriter netCdfWriter(fileName ,
    		   	   	   	   	   	   	 l_wavePropagationBlock[i][j]->getNx(),
    		   	   	   	   	   	   	 l_wavePropagationBlock[i][j]->getNy());
       //create the netCDF-file
       netCdfWriter.createNetCdfFile(l_wavePropagationBlock[i][j]->getDx(),
									 l_wavePropagationBlock[i][j]->getDy(),
						        	 l_originX + i*l_nX*l_dX,
						        	 l_originY + j*l_nY*l_dY
									 );
       int l_boundarySize[4];
       if (l_interpolationScheme == APPROX_TIME_SPACE)
    	   l_boundarySize[0] = l_boundarySize[1] = l_boundarySize[2] = l_boundarySize[3] = 1;
       else
    	   l_boundarySize[0] = l_boundarySize[1] = l_boundarySize[2] = l_boundarySize[3] = 2*l_rxy[i][j];
       netCdfWriter.writeBathymetry(l_wavePropagationBlock[i][j]->getBathymetry(), l_boundarySize);
       netCdfWriter.writeUnknowns(l_wavePropagationBlock[i][j]->getWaterHeight(),
							   	  l_wavePropagationBlock[i][j]->getDischarge_hu(),
							   	  l_wavePropagationBlock[i][j]->getDischarge_hv(),
							   	  l_boundarySize, (float) 0.);
     }
  }
#else
  // write the output at time zero
  l_wavePropagationBlock.writeVTKFileXML(generateFileName(l_baseName,0,0,0), l_nX, l_nY);
#endif

  // print some variables
  s_sweLogger.printCellSize(l_dX, l_dY);
  s_sweLogger.printNumberOfCells(l_nX, l_nY);
  s_sweLogger.printNumberOfBlocks(l_blockX, l_blockY);
  s_sweLogger.printTimeStepping(l_timeSteppingStrategy);
  s_sweLogger.printInterpolation(l_interpolationScheme);

  /**
   * Simulation.
   */
  // print the start message and reset the wall clock time
  s_sweLogger.printStartMessage();
  s_sweLogger.initWallClockTime(time(NULL));

  SWE_BlockManager mgr(l_wavePropagationBlock, l_blockX, l_blockY, l_interpolationScheme);

#ifdef BENCHMARKING
  // init benchmarking data receiver
  mgr.initBenchmarkingDataReceiver(std::string("single_wave_on_simple_beach"));

  /* add times and positions to collect benchmarking data
   *
   * water level profiles at:
   * t = 25(d/g)^1/2, t = 35(d/g)^1/2,
   * t = 45(d/g)^1/2, t = 55(d/g)^1/2,
   * t = 65(d/g)^1/2
   *
   * water level dynamics at (for some fixed y):
   * x/d = 0.25 and x/d = 9.95
   *
   * d = 1.
   * g = 1.
   *
   */
  mgr.addSpatialData(25.);
  mgr.addSpatialData(35.);
  mgr.addSpatialData(45.);
  mgr.addSpatialData(55.);
  mgr.addSpatialData(65.);

  mgr.addTimeSeries(0.25, 0.0);
  mgr.addTimeSeries(9.95, 0.0);
#endif

  //! simulation time.
  float l_t = 0.0;

  // loop over checkpoints
  for(int c=1; c<=l_numberOfCheckPoints; ++c) {
    // reset the cpu clock
    s_sweLogger.resetCpuClockToCurrentTime();

    // do time steps until next checkpoint is reached
    while( l_t < l_checkPoints[c] ) {
#ifndef CUDA
    if (l_timeSteppingStrategy == LTS)
      l_t += mgr.simulate(l_checkPoints[c] - l_checkPoints[c-1]);
    else
      l_t += mgr.simulate_gts(l_checkPoints[c] - l_checkPoints[c-1]);

#else
      // set values in ghost cells:
      l_wavePropgationBlock.setGhostLayer();

      // compute numerical flux on each edge
      l_wavePropgationBlock.computeNumericalFluxes();

      //! maximum allowed time step width.
      float l_maxTimeStepWidth = l_wavePropgationBlock.getMaxTimestep();

      // update the cell values
      l_wavePropgationBlock.updateUnknowns(l_maxTimeStepWidth);

      // update simulation time with time step width.
      l_t += l_maxTimeStepWidth;

#endif
      // print the current simulation time
      s_sweLogger.printSimulationTime(l_t);
    }

    // update the cpu time in the logger
    s_sweLogger.updateCpuTime();

    // print current simulation time
    s_sweLogger.printOutputTime(l_t);

#ifdef WRITENETCDF
    for(int i=0; i<l_blockX; ++i)
         for(int j=0; j<l_blockY; ++j) {
        	 //construct a NetCdfWriter
        	 io::NetCdfWriter netCdfWriter( generateFileName(l_baseName,0,i,j, ".nc"),
        		   l_wavePropagationBlock[i][j]->getNx(),
        		   l_wavePropagationBlock[i][j]->getNy());
        	 int l_boundarySize[4];
             if (l_interpolationScheme == APPROX_TIME_SPACE)
          	   l_boundarySize[0] = l_boundarySize[1] = l_boundarySize[2] = l_boundarySize[3] = 1;
             else
          	   l_boundarySize[0] = l_boundarySize[1] = l_boundarySize[2] = l_boundarySize[3] = 2*l_rxy[i][j];
        	 netCdfWriter.writeUnknowns(l_wavePropagationBlock[i][j]->getWaterHeight(),
        		   	   	   	   	   	    l_wavePropagationBlock[i][j]->getDischarge_hu(),
        		   	   	   	   	   	    l_wavePropagationBlock[i][j]->getDischarge_hv(),
        		   	   	   	   	   	    l_boundarySize, l_t );
         }
#else
    // write vtk output
    l_wavePropagationBlock.writeVTKFileXML(generateFileName(l_baseName,c,0,0), l_nX, l_nY);
#endif
  }

  /**
   * Finalize.
   */
  // write the statistics message
  s_sweLogger.printStatisticsMessage();

  // print the cpu time
  s_sweLogger.printCpuTime("CPU/GPU time");

  // print the wall clock time (includes plotting)
  s_sweLogger.printWallClockTime(time(NULL));

  return 0;
}
