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
 * Very basic setting of SWE, which uses a wave propagation solver with
 * an adaptively refined grid and an artificial scenario.
 * There is also the possibility of using the ASAGI library for setting the scenario.
 */

#include <cstdlib>
#include <string>
#include <mpi.h>

#include "../SWE_BlockAMR.hh"
#include "../SWE_BlockGhost.hh"
#include "../SWE_WavePropagationAMR.hh"

#include "../tools/Logger.hpp"
#include "../tools/help.hh"

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

#include "swe_mpiHelper.hh"

/**
 * Main program for the simulation on a single SWE_WavePropagationBlockAMR per MPI process.
 */
int main(int argc, char** argv) {
	/**
	 * Initialization.
	 */
	//! MPI Rank of a process.
	int l_mpiRank;
	//! number of MPI processes.
	int l_numberOfProcesses;

	// initialize MPI
	if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
		std::cerr << "MPI_Init failed." << std::endl;
	}

	// determine local MPI rank
	MPI_Comm_rank(MPI_COMM_WORLD, &l_mpiRank);
	// determine total number of processes
	MPI_Comm_size(MPI_COMM_WORLD, &l_numberOfProcesses);

	// initialize a logger for every MPI process
	tools::Logger l_sweLogger(l_mpiRank);

	// print the welcome message
	l_sweLogger.printWelcomeMessage();

	// set current wall clock time within the solver
	l_sweLogger.initWallClockTime(MPI_Wtime());
	//print the number of processes
	l_sweLogger.printNumberOfProcesses(l_numberOfProcesses);

	// check if the necessary command line input parameters are given
	if (argc != 9) {
		std::cout << "Aborting ... please provide proper input parameters."
				<< std::endl
				<< "Example: ./SWE_parallel 10 25 200 300 /work/openmp_out /work/config <time stepping> <interpolation scheme>"
				<< std::endl
				<< "\tfor 2000 * 7500 cells with 10 * 25 blocks of size 200 * 300"
				<< std::endl
				<< "\ttime stepping can be GTS / LTS for global/local time stepping"
				<< std::endl
				<< "\tinterpolation scheme can be 0 / 1 / 2 for approximate time-space / time-space / space interpolation"
				<< std::endl;
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

	if (!configFile) {
		cout << endl << "Failed to open file" << l_configName;
		return 1;
	}

	int** l_rxy = new int*[l_blockX];
	for (int i = 0; i < l_blockX; ++i)
		l_rxy[i] = new int[l_blockY];
	for (int j = l_blockY - 1; j >= 0; --j)
		for (int i = 0; i < l_blockX; ++i)
			configFile >> l_rxy[i][j];

	configFile.close();

	/**
	 * TODO:
	 * do not read refinement levels from file
	 * instead, get them from the other processes:
	 * 1. allgather
	 * MPI_Allgather(&own_refinement_level, 1, MPI_INT, l_rxy, 1, MPI_INT, MPI_COMM_WORLD);
	 * 2. communicate only with neighbours
	 * send/recv
	 */

	// read time-stepping strategy
	std::string l_ts = std::string(argv[7]);
	if (l_ts.compare("GTS") == 0) l_timeSteppingStrategy = GTS;
	else if (l_ts.compare("LTS") == 0) l_timeSteppingStrategy = LTS;
	else { cout << endl << "Wrong time-stepping strategy" << l_ts << endl; return 1; }

	// read interpolation scheme
	int l_is = atoi(argv[8]);
	switch (l_is) {
	case 0: l_interpolationScheme = APPROX_TIME_SPACE; break;
	case 1: l_interpolationScheme = TIME_SPACE; break;
	case 2: l_interpolationScheme = SPACE; break;
	default: cout << endl << "Wrong interpolation scheme" << l_is << endl; return 1;
	}

	//! local position of each MPI process in x- and y-direction.
	int l_blockPositionX, l_blockPositionY;

	// determine local block coordinates of each SWE_Block
	l_blockPositionX = l_mpiRank / l_blockY;
	l_blockPositionY = l_mpiRank % l_blockY;

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
  SWE_BathymetryDamBreakScenario l_scenario;
//	SWE_SingleWaveOnSimpleBeach l_scenario;
#endif
	//! number of checkpoints for visualization (at each checkpoint in time, an output file is written).
	int l_numberOfCheckPoints = 20;

	//! refinement level in x- and y-direction
	float l_rX, l_rY;

	l_rX = l_rxy[l_blockPositionX][l_blockPositionY];
	l_rY = l_rxy[l_blockPositionX][l_blockPositionY];

	//! number of grid cells in x- and y-direction per process.
	int l_nXLocal, l_nYLocal;

	l_nXLocal = l_nX*l_rX;
	l_nYLocal = l_nY*l_rY;

	//! size of a single cell in x- and y-direction
	float l_dX, l_dY;

	// compute the size of a single cell
	l_dX = (l_scenario.getBoundaryPos(BND_RIGHT) - l_scenario.getBoundaryPos(BND_LEFT)) / (l_blockX * l_nX);
	l_dY = (l_scenario.getBoundaryPos(BND_TOP) - l_scenario.getBoundaryPos(BND_BOTTOM)) / (l_blockY * l_nY);

	//! origin of the simulation domain in x- and y-direction
	float l_originX, l_originY;

	l_originX = l_scenario.getBoundaryPos(BND_LEFT) + l_blockPositionX * l_nX * l_dX;
	l_originY = l_scenario.getBoundaryPos(BND_BOTTOM) + l_blockPositionY * l_nY * l_dY;

	//! compute the number of ghost cells depending on the interpolation scheme
	if (l_interpolationScheme == APPROX_TIME_SPACE)
		l_nghosts = 1;
	else
		l_nghosts = 2 * l_rX;

	//! initialise it to given scenario
	SWE_WavePropagationAMR l_wavePropagationBlock(l_originX, l_originY,
												  l_nXLocal, l_nYLocal,
												  l_dX/l_rX, l_dY/l_rY,
												  l_rX, l_rY,
												  l_nghosts, // size of ghost layer
												  l_interpolationScheme); // interpolation scheme for ghost layers
	l_wavePropagationBlock.initScenario(l_scenario, true);

	int l_leftNeighbourRefinementLevel=0;
	int l_rightNeighbourRefinementLevel=0;
	int l_bottomNeighbourRefinementLevel=0;
	int l_topNeighbourRefinementLevel=0;

	/*
	 * Connect SWE blocks at boundaries
	 */
	// left and right boundaries
	l_sweLogger.printString("Connecting SWE blocks at left boundaries.");
	if (l_blockPositionX > 0) {
		l_leftNeighbourRefinementLevel = l_rxy[l_blockPositionX-1][l_blockPositionY];
		l_wavePropagationBlock.setNeighbourRefinementLevel(l_leftNeighbourRefinementLevel, BND_LEFT);
	}
	SWE_BlockGhost* l_leftInflow = l_wavePropagationBlock.grabGhostLayer(BND_LEFT);
	SWE_BlockGhost* l_leftOutflow = l_wavePropagationBlock.registerCopyLayer(BND_LEFT);
	if (l_blockPositionX == 0)
		l_wavePropagationBlock.setBoundaryType(BND_LEFT, l_scenario.getBoundaryType(BND_LEFT));

	l_sweLogger.printString("Connecting SWE blocks at right boundaries.");
	if (l_blockPositionX < l_blockX - 1) {
		l_rightNeighbourRefinementLevel = l_rxy[l_blockPositionX+1][l_blockPositionY];
		l_wavePropagationBlock.setNeighbourRefinementLevel(l_rightNeighbourRefinementLevel, BND_RIGHT);
	}
	SWE_BlockGhost* l_rightInflow = l_wavePropagationBlock.grabGhostLayer(BND_RIGHT);
	SWE_BlockGhost* l_rightOutflow = l_wavePropagationBlock.registerCopyLayer(BND_RIGHT);
	if (l_blockPositionX == l_blockX - 1)
		l_wavePropagationBlock.setBoundaryType(BND_RIGHT, l_scenario.getBoundaryType(BND_RIGHT));

	// bottom and top boundaries
	l_sweLogger.printString("Connecting SWE blocks at bottom boundaries.");
	if (l_blockPositionY > 0) {
		l_bottomNeighbourRefinementLevel = l_rxy[l_blockPositionX][l_blockPositionY-1];
		l_wavePropagationBlock.setNeighbourRefinementLevel(l_bottomNeighbourRefinementLevel, BND_BOTTOM);
	}
	SWE_BlockGhost* l_bottomInflow = l_wavePropagationBlock.grabGhostLayer(BND_BOTTOM);
	SWE_BlockGhost* l_bottomOutflow = l_wavePropagationBlock.registerCopyLayer(BND_BOTTOM);
	if (l_blockPositionY == 0)
		l_wavePropagationBlock.setBoundaryType(BND_BOTTOM, l_scenario.getBoundaryType(BND_BOTTOM));

	l_sweLogger.printString("Connecting SWE blocks at top boundaries.");
	if (l_blockPositionY < l_blockY - 1) {
		l_topNeighbourRefinementLevel = l_rxy[l_blockPositionX][l_blockPositionY+1];
		l_wavePropagationBlock.setNeighbourRefinementLevel(l_topNeighbourRefinementLevel, BND_TOP);
	}
	SWE_BlockGhost* l_topInflow = l_wavePropagationBlock.grabGhostLayer(BND_TOP);
	SWE_BlockGhost* l_topOutflow = l_wavePropagationBlock.registerCopyLayer(BND_TOP);
	if (l_blockPositionY == l_blockY - 1)
		l_wavePropagationBlock.setBoundaryType(BND_TOP, l_scenario.getBoundaryType(BND_TOP));

	/*
	 * Create the receive data types for each boundary -> ghost layer type
	 * The grid is stored column wise in memory:
	 *
	 *        ************************** . . . **********
	 *        *       *  ny+2 *2(ny+2)*         * (ny+1)*
	 *        *  ny+1 * +ny+1 * +ny+1 *         * (ny+2)*
	 *        *       *       *       *         * +ny+1 *
	 *        ************************** . . . **********
	 *        *       *       *       *         *       *
	 *        .       .       .       .         .       .
	 *        .       .       .       .         .       .
	 *        .       .       .       .         .       .
	 *        *       *       *       *         *       *
	 *        ************************** . . . **********
	 *        *       *  ny+2 *2(ny+2)*         * (ny+1)*
	 *        *   1   *   +1  *   +1  *         * (ny+2)*
	 *        *       *       *       *         *   +1  *
	 *        ************************** . . . **********
	 *        *       *  ny+2 *2(ny+2)*         * (ny+1)*
	 *        *   0   *   +0  *   +0  *         * (ny+2)*
	 *        *       *       *       *         *   +0  *
	 *        ************************** . . . ***********
	 *
	 *
	 *  -> The stride for a row is ny+2*nghosts, because we have to jump over a whole column
	 *     for every row-element. This holds only in the CPU-version, in CUDA a buffer is implemented.
	 *     See SWE_BlockCUDA.hh/.cu for details.
	 *  -> The stride for a column is 1, because we can access the elements linear in memory.
	 */
	//! MPI block of rows: l_nXLocal+2*l_nghosts blocks, l_nghosts elements per block, stride of l_nYLocal+2*l_nghosts
	MPI_Datatype l_mpiRows;
	MPI_Type_vector(l_nXLocal+2*l_nghosts,	l_nghosts,				l_nYLocal+2*l_nghosts,	MPI_FLOAT, &l_mpiRows);
	MPI_Type_commit(&l_mpiRows);

	//! MPI block of columns: l_nghosts blocks, l_nYLocal+2*l_nghosts elements per block, stride of l_nYLocal+2*l_nghosts
	MPI_Datatype l_mpiCols;
	MPI_Type_contiguous(l_nghosts*(l_nYLocal+2*l_nghosts), MPI_FLOAT, &l_mpiCols);
	MPI_Type_commit(&l_mpiCols);

	/**
	 * If the blocks have the same resolution, then a proxy block is used as a copy layer.
	 * If the blocks have different resolutions, then a newly allocated block is used as a copy layer.
	 * Solution: either create different datatypes types for each boundary, depending on the
	 * implementation of the copy layer, or allocate the copy layer even in the first case (same resolution blocks)
	 * The first option is implemented below.
	 * (TODO) Note: these types, as well as the ones above, have to be destroyed and created again
	 * if the resolution of a block changes during the program execution
	 */
	MPI_Datatype l_mpiSendType[4];
	// left boundary type
	if (l_leftNeighbourRefinementLevel == l_rX)
		l_mpiSendType[BND_LEFT] = l_mpiCols;
	else {
		MPI_Type_contiguous(l_leftOutflow->nx * l_leftOutflow->ny, MPI_FLOAT, &l_mpiSendType[BND_LEFT]);
		MPI_Type_commit(&l_mpiSendType[BND_LEFT]);
	}

	// right boundary type
	if (l_rightNeighbourRefinementLevel == l_rX)
		l_mpiSendType[BND_RIGHT] = l_mpiCols;
	else {
		MPI_Type_contiguous(l_rightOutflow->nx * l_rightOutflow->ny, MPI_FLOAT, &l_mpiSendType[BND_RIGHT]);
		MPI_Type_commit(&l_mpiSendType[BND_RIGHT]);
	}

	// bottom boundary type
	if (l_bottomNeighbourRefinementLevel == l_rX)
		l_mpiSendType[BND_BOTTOM] = l_mpiRows;
	else {
		MPI_Type_contiguous(l_bottomOutflow->nx * l_bottomOutflow->ny, MPI_FLOAT, &l_mpiSendType[BND_BOTTOM]);
		MPI_Type_commit(&l_mpiSendType[BND_BOTTOM]);
	}

	// top boundary type
	if (l_topNeighbourRefinementLevel == l_rX)
		l_mpiSendType[BND_TOP] = l_mpiRows;
	else {
		MPI_Type_contiguous(l_topOutflow->nx * l_topOutflow->ny, MPI_FLOAT, &l_mpiSendType[BND_TOP]);
		MPI_Type_commit(&l_mpiSendType[BND_TOP]);
	}

	//! MPI ranks of the neighbors
	int l_leftNeighborRank, l_rightNeighborRank, l_bottomNeighborRank, l_topNeighborRank;

	// compute MPI ranks of the neighbour processes
	l_leftNeighborRank   = (l_blockPositionX > 0) ? l_mpiRank-l_blockY : MPI_PROC_NULL;
	l_rightNeighborRank  = (l_blockPositionX < l_blockX-1) ? l_mpiRank+l_blockY : MPI_PROC_NULL;
	l_bottomNeighborRank = (l_blockPositionY > 0) ? l_mpiRank-1 : MPI_PROC_NULL;
	l_topNeighborRank    = (l_blockPositionY < l_blockY-1) ? l_mpiRank+1 : MPI_PROC_NULL;

	// print the MPI grid
	l_sweLogger.cout() << "neighbors: "
					 << l_leftNeighborRank << " (left), "
					 << l_rightNeighborRank << " (right), "
					 << l_bottomNeighborRank << " (bottom), "
					 << l_topNeighborRank << " (top)" << std::endl;

	//! time when the simulation ends.
	float l_endSimulation = l_scenario.endSimulation();

	//! checkpoints when output files are written.
	float* l_checkPoints = new float[l_numberOfCheckPoints + 1];

	// compute the checkpoints in time
	for (int cp = 0; cp <= l_numberOfCheckPoints; ++cp) {
		l_checkPoints[cp] = cp * (l_endSimulation / l_numberOfCheckPoints);
	}

	// intially exchange ghost and copy layers
//	exchangeLeftRightGhostLayers( &l_wavePropagationBlock, 0.0, 0.0, l_timeSteppingStrategy,
//								  l_leftNeighborRank,  l_leftInflow,  l_leftOutflow, l_leftNeighbourRefinementLevel,
//								  l_rightNeighborRank, l_rightInflow, l_rightOutflow,l_rightNeighbourRefinementLevel,
//								  l_rX, l_mpiCols, l_mpiSendType );
//	exchangeBottomTopGhostLayers( &l_wavePropagationBlock, 0.0, 0.0, l_timeSteppingStrategy,
//			  	  	  	  	  	  l_bottomNeighborRank, l_bottomInflow, l_bottomOutflow,l_bottomNeighbourRefinementLevel,
//								  l_topNeighborRank,    l_topInflow,    l_topOutflow,	l_topNeighbourRefinementLevel,
//								  l_rX, l_mpiRows, l_mpiSendType );
	// set values in ghost cells
	l_wavePropagationBlock.setGhostLayer();

	// write the output at time zero
	l_sweLogger.printOutputTime(0);

#ifdef WRITENETCDF
	std::string l_fileName = generateFileName(l_baseName,l_blockPositionX,l_blockPositionY);
	io::NetCdfWriter netCdfWriter(l_fileName ,
			l_wavePropagationBlock.getNx(),
			l_wavePropagationBlock.getNy());

	//create the netCDF-file
	netCdfWriter.createNetCdfFile(l_wavePropagationBlock.getDx(),
								  l_wavePropagationBlock.getDy(),
								  l_originX, l_originY);
	int l_boundarySize[4];
	if (l_interpolationScheme == APPROX_TIME_SPACE)
		l_boundarySize[0] = l_boundarySize[1] = l_boundarySize[2] = l_boundarySize[3] = 1;
	else
		l_boundarySize[0] = l_boundarySize[1] = l_boundarySize[2] = l_boundarySize[3] = 2*l_rX;

	netCdfWriter.writeBathymetry(l_wavePropagationBlock.getBathymetry(), l_boundarySize);
	netCdfWriter.writeUnknowns(l_wavePropagationBlock.getWaterHeight(),
							   l_wavePropagationBlock.getDischarge_hu(),
							   l_wavePropagationBlock.getDischarge_hv(),
							   l_boundarySize, (float) 0.);
#else
	  l_wavePropagationBlock.writeVTKFileXML( generateFileName(l_baseName,0,l_blockPositionX,l_blockPositionY),
	                                         l_nXLocal*l_blockPositionX, l_nYLocal*l_blockPositionY);
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

	//! if local time-stepping is used, create an MPI group for each refinement level
	MPI_Comm l_refinementLevelComm;
	if (l_timeSteppingStrategy == LTS) {
		MPI_Comm_split( MPI_COMM_WORLD, l_rX, l_mpiRank, &l_refinementLevelComm );
		/**
		 * TODO: how to destroy a communicator:
		 * MPI_Comm_free( &l_refinementLevelComm );
		 * if (l_refinementLevelComm != MPI_COMM_NULL) {
		 * 	 printf( "Freed comm was not set to COMM_NULL\n" );
		 * }
		 */
	}

	//! simulation time.
	float l_t = 0.0;

	// loop over checkpoints
	for (int c = 1; c <= l_numberOfCheckPoints; ++c) {
		// reset the cpu clock
		s_sweLogger.resetCpuClockToCurrentTime();

		if (l_timeSteppingStrategy == GTS) { // global time-stepping

			// do time steps until next checkpoint is reached
			while( l_t < l_checkPoints[c] ) {
				// exchange ghost and copy layers

				// left - right direction
				exchangeLeftRightGhostLayers( &l_wavePropagationBlock, l_t, l_t, l_timeSteppingStrategy,
											  l_leftNeighborRank,  l_leftInflow,  l_leftOutflow,	l_leftNeighbourRefinementLevel,
											  l_rightNeighborRank, l_rightInflow, l_rightOutflow,	l_rightNeighbourRefinementLevel,
											  l_rX, l_mpiCols, l_mpiSendType );
				// bottom - top direction
				exchangeBottomTopGhostLayers( &l_wavePropagationBlock, l_t, l_t, l_timeSteppingStrategy,
											  l_bottomNeighborRank, l_bottomInflow, l_bottomOutflow,l_bottomNeighbourRefinementLevel,
											  l_topNeighborRank,    l_topInflow,    l_topOutflow,	l_topNeighbourRefinementLevel,
											  l_rX, l_mpiRows, l_mpiSendType );

				// reset the cpu clock
				l_sweLogger.resetCpuClockToCurrentTime();

				// set values in ghost cells
				l_wavePropagationBlock.setGhostLayer();

				// compute numerical flux on each edge
				l_wavePropagationBlock.computeNumericalFluxes();

				//! maximum allowed time step width within a block.
				float l_maxTimeStepWidth = l_wavePropagationBlock.getMaxTimestep();

				//! maximum allowed time steps of all blocks
				float l_maxTimeStepWidthGlobal;

				// determine smallest time step of all blocks
				MPI_Allreduce(&l_maxTimeStepWidth, &l_maxTimeStepWidthGlobal, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);

				// update the cell values
				l_wavePropagationBlock.updateUnknowns(l_maxTimeStepWidthGlobal);

				// update simulation time with time step width.
				l_t += l_maxTimeStepWidthGlobal;
			}
			// print the current simulation time
			//l_sweLogger.printSimulationTime(l_t);
		} else { // local time-stepping
			if (l_interpolationScheme == SPACE) {
				// do time steps until next checkpoint is reached
				while( l_t < l_checkPoints[c] ) {

					// exchange ghost and copy layers
					// left - right direction
					exchangeLeftRightGhostLayers( &l_wavePropagationBlock, l_t, l_t, l_timeSteppingStrategy,
												  l_leftNeighborRank,  l_leftInflow,  l_leftOutflow,	l_leftNeighbourRefinementLevel,
												  l_rightNeighborRank, l_rightInflow, l_rightOutflow,	l_rightNeighbourRefinementLevel,
												  l_rX, l_mpiCols, l_mpiSendType );
					// bottom - top direction
					exchangeBottomTopGhostLayers( &l_wavePropagationBlock, l_t, l_t, l_timeSteppingStrategy,
												  l_bottomNeighborRank, l_bottomInflow, l_bottomOutflow,l_bottomNeighbourRefinementLevel,
												  l_topNeighborRank,    l_topInflow,    l_topOutflow,	l_topNeighbourRefinementLevel,
												  l_rX, l_mpiRows, l_mpiSendType );

					// reset the cpu clock
					l_sweLogger.resetCpuClockToCurrentTime();

					// set values in ghost cells
					l_wavePropagationBlock.setGhostLayer();

					// reset computational domain to include all but one ghost layers
					l_wavePropagationBlock.resetComputationalDomainMax();

					// number of time-steps performed
					int l_numTs = 0;

					//! maximum allowed time step width within a block multiplied by the refinement level
					float l_maxTimeStepWidth = l_wavePropagationBlock.getMaxTimestep() * l_rX;

					//! maximum allowed time steps of all blocks on a refinement level
					float l_maxTimeStepWidthGlobal;

					// determine smallest time step of all blocks
					MPI_Allreduce(&l_maxTimeStepWidth, &l_maxTimeStepWidthGlobal, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);

					// do time-stepping until all ghost layers become invalid
					while (l_numTs < l_rX) {
						// set values in ghost cells
						l_wavePropagationBlock.setGhostLayer();

						// compute numerical flux on each edge
						l_wavePropagationBlock.computeNumericalFluxes();

						// update the cell values
						l_wavePropagationBlock.updateUnknowns(l_maxTimeStepWidthGlobal/l_rX);

						// decrease the computational domain
						l_wavePropagationBlock.decreaseComputationalDomain();

						// increase number of time-steps executed since the last ghost cell exchange
						l_numTs++;
					}

					// update simulation time
					l_t += l_maxTimeStepWidthGlobal;
				}
			} else { // l_interpolationScheme is APPROX_TIME_SPACE or TIME_SPACE
				// do time steps until next checkpoint is reached
				// use fixed time-steps for now
				while( l_t < l_checkPoints[c] ) {
					// reset the cpu clock
					l_sweLogger.resetCpuClockToCurrentTime();
					// simulate one coarse time-step
					l_t += simulateLTSTimeSpace(&l_wavePropagationBlock, l_t, l_checkPoints[c],
							  l_leftNeighborRank, l_leftInflow, l_leftOutflow, l_leftNeighbourRefinementLevel,
							  l_rightNeighborRank, l_rightInflow, l_rightOutflow, l_rightNeighbourRefinementLevel,
							  l_bottomNeighborRank, l_bottomInflow, l_bottomOutflow, l_bottomNeighbourRefinementLevel,
							  l_topNeighborRank, l_topInflow, l_topOutflow, l_topNeighbourRefinementLevel,
							  l_rX, l_mpiCols, l_mpiRows, l_mpiSendType, l_refinementLevelComm);
				}
				// print the current simulation time
				l_sweLogger.printSimulationTime(l_t);
			}
		}

		// update the cpu time in the logger
		s_sweLogger.updateCpuTime();

		// print current simulation time
		s_sweLogger.printOutputTime(l_t);

#ifdef WRITENETCDF
	std::string l_fileName = generateFileName(l_baseName,l_blockPositionX,l_blockPositionY);
	io::NetCdfWriter netCdfWriter(l_fileName ,
			l_wavePropagationBlock.getNx(),
			l_wavePropagationBlock.getNy());
	int l_boundarySize[4];
	if (l_interpolationScheme == APPROX_TIME_SPACE)
		l_boundarySize[0] = l_boundarySize[1] = l_boundarySize[2] = l_boundarySize[3] = 1;
	else
		l_boundarySize[0] = l_boundarySize[1] = l_boundarySize[2] = l_boundarySize[3] = 2*l_rX;
	netCdfWriter.writeUnknowns(l_wavePropagationBlock.getWaterHeight(),
							   l_wavePropagationBlock.getDischarge_hu(),
							   l_wavePropagationBlock.getDischarge_hv(),
							   l_boundarySize, l_t );
#else
	  l_wavePropagationBlock.writeVTKFileXML( generateFileName(l_baseName, c,l_blockPositionX,l_blockPositionY),
	                                         l_nXLocal*l_blockPositionX,
	                                         l_nYLocal*l_blockPositionY);
#endif
	}

	/**
	 * Finalize.
	 */
	// write the statistics message
	s_sweLogger.printStatisticsMessage();

	// print the cpu time
	s_sweLogger.printCpuTime("CPU time");

	// print the wall clock time (includes plotting)
	s_sweLogger.printWallClockTime(time(NULL));

	// finalize MPI execution
	MPI_Finalize();

	return 0;
}
