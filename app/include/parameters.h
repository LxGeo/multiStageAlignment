#pragma once
#include "defs.h"
#include "cli/base_parameters.h"
#include "CLI/CLI.hpp"


namespace LxGeo
{
	namespace MSA
	{
		class Parameters: public GeometryFactoryShared::baseParameters
		{
		public:
			Parameters(int argc, char* argv[]);

			~Parameters();

			bool initialized();

		protected:
			void init();

		public:

			std::string imd1_path;

		};

		extern Parameters* params;
	}
}
