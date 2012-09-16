
#define CONFIG_DEFAULT_FLOATING_POINT_TYPE float

#include "CHyperbolicTypes.hpp"
#include "../../datasets_common/CParameters_Datasets.hpp"

typedef CParameters_Datasets	CParameters;


#include "../CDatasets.hpp"


int main(int argc, char *argv[])
{
	int verbosity_level;

	CParameters_Datasets cParameters_Datasets;
	CDatasets cDatasets(cParameters_Datasets, verbosity_level);

	return 1;
}
