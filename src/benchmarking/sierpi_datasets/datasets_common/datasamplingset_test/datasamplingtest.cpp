#include <iostream>

#include "../CDataSamplingSet_MultiResolution.hpp"


void run()
{
	typedef float T;

	int array_size_x = 3;
	int array_size_y = 6;

	float domain_start_x = -10;
	float domain_start_y = -20;

	float domain_length_x = 100;
	float domain_length_y = 200;

	std::cout << "array_size_x: " << array_size_x << std::endl;
	std::cout << "array_size_y: " << array_size_y << std::endl;
	std::cout << "domain_start_x: " << domain_start_x << std::endl;
	std::cout << "domain_start_y: " << domain_start_y << std::endl;
	std::cout << "domain_length_x: " << domain_length_x << std::endl;
	std::cout << "domain_length_y: " << domain_length_y << std::endl;


	CDataSamplingSet_MultiResolution<T> d;
	T *data = d.allocateMultiresolutionLevels(array_size_x, array_size_y);

	/**
	 * dummy values for array
	 */
	auto arrayValue = [](float x, float y) -> float
	{
		return x+y;
	};

	/*
	 * fill array with values
	 */
	for (int y = 0; y < array_size_y; y++)
		for (int x = 0; x < array_size_x; x++)
			data[x+y*array_size_x] = arrayValue(x, y);

	/*
	 * setup domain parameters
	 */
	d.setupDomainParameters(domain_start_x, domain_start_y, domain_start_x+domain_length_x, domain_start_y+domain_length_y);

	/*
	 * determine number of multi-resolution levels
	 */
	int max_levels = d.getLevels();
	d.setupMultiresolutionMap();


	float cell_size_x = domain_length_x / (float)array_size_x;
	float cell_size_y = domain_length_y / (float)array_size_y;


	float sampling_displ_x = 0.0;
	float sampling_displ_y = 0.0;

	/*
	 * output values for levels
	 */
	for (int l = 0; l < max_levels; l++)
	{
		std::cout << "Level " << l << std::endl;

		std::cout << "RAW Array:" << std::endl;
		d.debugPrintArray(l);

		std::cout << "Sampled Array:" << std::endl;
		for (int y = 0; y < array_size_y; y++)
		{
			for (int x = 0; x < array_size_x; x++)
			{
				float sx = ((float)x*cell_size_x) + domain_start_x + cell_size_x*0.5f + sampling_displ_x;
				float sy = ((float)y*cell_size_y) + domain_start_y + cell_size_y*0.5f + sampling_displ_y;

				float value = d.getFloat2D(	sx, sy, l);
				std::cout << value << "\t";
			}

			std::cout << std::endl;
		}
	}
}


int main()
{
	run();

	return 0;
}
