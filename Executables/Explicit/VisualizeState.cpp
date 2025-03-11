/**
 * @file VisualizeState.cpp
 * @brief Simple general executable to visualize a state file for a model
 */

// System includes
//

// Project includes
//
#include "Environment/QuICCEnv.hpp"
#include "QuICC/Timers/StageTimer.hpp"
#include "Profiler/Interface.hpp"
#include "QuICC/Model/RunApplication.hpp"
#include "QuICC/Model/VisualizationGeneratorFactory.hpp"
#include "Model/Boussinesq/Sphere/RTC/Explicit/PhysicalModel.hpp"
#include "Model/Boussinesq/Sphere/RTC/IRTCVisualization.hpp"

/**
 * @brief Setup and run the simulation
 */
int run()
{
   namespace model_ns = QuICC::Model::Boussinesq::Sphere::RTC;
   typedef model_ns::Explicit::PhysicalModel<model_ns::IRTCVisualization> Model;
   int status = QuICC::run_application<QuICC::VisualizationGeneratorFactory,Model>();

   return status;
}

/**
 * @brief General main, setting up MPI if required
 *
 * The actual program is in run to make sure MPI initialisations
 * are called before anything else and finalisation after destruction
 */
int main(int argc, char* argv[])
{
   // Environment fixture
   QuICC::QuICCEnv();

   // Initialise everything that can't be done inside a class
   QuICC::StageTimer::allowIo(QuICC::QuICCEnv().allowsIO());
   QuICC::Profiler::Initialize();

   // Compute simulation
   QuICC::Profiler::RegionStart<1>("Walltime");
   auto code = run();
   QuICC::Profiler::RegionStop<1>("Walltime");

   // Finalise everything that can't be done inside a class
   QuICC::Profiler::Finalize();

   return code;
}
