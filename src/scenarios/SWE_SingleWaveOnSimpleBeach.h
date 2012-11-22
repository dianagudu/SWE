#include "../benchmarking/sierpi_datasets/subsimulation_tsunami/tsunami_benchmarks/CSingleWaveOnSimpleBeach.hpp"

/**
 * Scenario "Single wave on simple beach":
 * uniform bathymetry for a part of the domain,
 * sloping beach in the other part of the domain,
 * single wave
 */
class SWE_SingleWaveOnSimpleBeach : public SWE_Scenario {

  CSingleWaveOnSimpleBeach<float> singleWaveOnSimpleBeach;

  public:

	/** Get the bathymetry in a given point in the domain
	 *
	 * @param x point coordinate in x-direction
	 * @param y point coordinate in y-direction
	 * @return bathymetry value at given point
	 */
  	float getBathymetry(float x, float y) {
      return singleWaveOnSimpleBeach.p_getBathymetryData(x, y, 0);
    };

    /** Get the water height
     * in a given point in the domain
	 *
	 * @param x point coordinate in x-direction
	 * @param y point coordinate in y-direction
	 * @return water height value in given point
	 */
	float getWaterHeight( float x, float y ) {
		return singleWaveOnSimpleBeach.p_getWaterSufaceData(x, y, 0);
    }

    /** Get the velocity in x-direction
     * in a given point in the domain
	 *
	 * @param x point coordinate in x-direction
	 * @param y point coordinate in y-direction
	 * @return velocity value in x-direction
	 */
    float getVeloc_u(float x, float y) {
    	if (singleWaveOnSimpleBeach.p_getWaterSufaceData(x, y, 0) < 0.000001) return 0;
		return singleWaveOnSimpleBeach.p_getMomentum(x, y, 0) / singleWaveOnSimpleBeach.p_getWaterSufaceData(x, y, 0);
	}

    /** Get the velocity in y-direction
     * in a given point in the domain
	 *
	 * @param x point coordinate in x-direction
	 * @param y point coordinate in y-direction
	 * @return velocity value in y-direction
	 */
	float getVeloc_v(float x, float y) {
    	return 0;
	}

    /**
     * Get the ending time of the simulation
     * @return time when simulation ends
     */
	float endSimulation() { return 120.; };

    /** Get the boundary type
	 *
	 * @param i_edge which edge
	 * @return type of the given edge
	 */
	BoundaryType getBoundaryType(BoundaryEdge edge) {
    	return OUTFLOW;
    }

    /** Get the boundary positions
	 *
	 * @param i_edge which edge
	 * @return value in the corresponding dimension
	 */
    float getBoundaryPos(BoundaryEdge edge) {
    	switch (edge) {
		case BND_LEFT: return -5.0f;
		case BND_RIGHT: return 75.0f;
		case BND_BOTTOM: return -10.0f;
		case BND_TOP: return 10.0f;
		default: return 0.0;
		}
	};
};
