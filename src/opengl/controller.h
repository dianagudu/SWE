#ifndef CONTROLLER_H
#define CONTROLLER_H
// =====================================================================
// This file is part of SWE_CUDA (see file SWE_Block.cu for details).
// 
// Copyright (C) 2010,2011 Tobias Schnabel
// 
// SWE_CUDA is free software: you can redristribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// SWE_CUDA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with SWE_CUDA.  If not, see <http://www.gnu.org/licenses/>.
// =====================================================================
#include <SDL.h>
#include "simulation.h"
#include "visualization.h"

class Controller {
public:
	Controller(Simulation* sim, Visualization* vis);
	
	// Process new events
	bool handleEvents();

	// Return if window is enabled
	bool hasFocus();
	// Return if program is paused
	bool isPaused();
private:
	// Internal functions
	bool isActive;
	bool done;
	bool paused;
	bool allowStep;

	// References to other classes
	Simulation* simulation;
	Visualization* visualization;
	
	// Handle keyboard events
	bool handleKeyPress( SDL_keysym *keysym);
};

#endif