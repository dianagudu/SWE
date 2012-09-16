#include "../benchmarking/sierpi_datasets/subsimulation_tsunami/tsunami_benchmarks/CSingleWaveOnSimpleBeach.hpp"

class SWE_SingleWaveOnSimpleBeach : public SWE_Scenario {

  CSingleWaveOnSimpleBeach<float> singleWaveOnSimpleBeach;

  public:

    float getBathymetry(float x, float y) {
      return singleWaveOnSimpleBeach.p_getBathymetryData(x, y, 0);
    };

	float endSimulation() { return 120.; };

    float getWaterHeight( float x, float y ) {
		return singleWaveOnSimpleBeach.p_getWaterSufaceData(x, y, 0);
    }

    float getVeloc_u(float x, float y) {
    	if (singleWaveOnSimpleBeach.p_getWaterSufaceData(x, y, 0) < 0.000001) return 0;
		return singleWaveOnSimpleBeach.p_getMomentum(x, y, 0) / singleWaveOnSimpleBeach.p_getWaterSufaceData(x, y, 0);
	}

    float getVeloc_v(float x, float y) {
    	return 0;
	}

    BoundaryType getBoundaryType(BoundaryEdge edge) { return OUTFLOW; }

    float getBoundaryPos(BoundaryEdge edge) {
    	singleWaveOnSimpleBeach.outputVerboseInformation();
    	switch (edge) {
		case BND_LEFT: return -10.0f;
		case BND_RIGHT: return 50.0f;
		case BND_BOTTOM: return -10.0f;
		case BND_TOP: return 10.0f;
		}
	};
};
