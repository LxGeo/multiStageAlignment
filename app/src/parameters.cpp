#include "parameters.h"

namespace LxGeo
{
	namespace MSA
	{
		Parameters::Parameters(int argc, char* argv[])
		{
			init();
			CLI::App app{ "MultiStage alignment" };
			app.add_option("--imd1", imd1_path, "Metadata file respective to template image (template shapefile)!")->check(CLI::ExistingFile);

			try {
				\
					(app).parse((argc), (argv));\
			}
			catch (const CLI::ParseError& e) {
				\
					(app).exit(e);\
			}
			//app.parse(argc, argv);
		}


		Parameters::~Parameters()
		{
		}


		bool Parameters::initialized()
		{
			return !imd1_path.empty();
		}


		void Parameters::init()
		{
			imd1_path.clear();
		}

		
		Parameters* params = nullptr;
	}
}