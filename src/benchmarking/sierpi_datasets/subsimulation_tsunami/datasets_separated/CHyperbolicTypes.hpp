#ifndef CHYPERBOLIC_TYPES_HPP
#define CHYPERBOLIC_TYPES_HPP


class CHyperbolicTypes
{
public:
	class CSimulationTypes
	{
	public:
		typedef CONFIG_DEFAULT_FLOATING_POINT_TYPE T;

		/**
		 * Degree of freedoms for one point in tsunami simulation
		 */
		class CNodeData
		{
		public:
			typedef CONFIG_DEFAULT_FLOATING_POINT_TYPE T;
			T h, hu, hv, b;
		};
	};
};

#endif
