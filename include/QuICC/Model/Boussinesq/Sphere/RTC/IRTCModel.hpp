/**
 * @file IRTCModel.hpp
 * @brief Implementation of the Boussinesq rotating thermal convection in a sphere (Toroidal/Poloidal formulation)
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_IRTCMODEL_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_IRTCMODEL_HPP

// Configuration includes
//

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "QuICC/Simulation/Simulation.hpp"
#include "QuICC/Generator/StateGenerator.hpp"
#include "QuICC/Generator/VisualizationGenerator.hpp"
#include "QuICC/Model/IPhysicalPyModel.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace RTC {

   /**
    * @brief Implementation of the Boussinesq rotating thermal convection sphere model (Toroidal/Poloidal formulation)
    */
   class IRTCModel: public IPhysicalPyModel
   {
      public:
         /**
          * @brief Constructor
          */
         IRTCModel() = default;

         /**
          * @brief Destructor
          */
         virtual ~IRTCModel() = default;

         /// Formulation used for vector fields
         virtual VectorFormulation::Id SchemeFormulation() override;

         /**
          * @brief Add the required equations
          *
          * @param spSim   Shared simulation object
          */
         virtual void addEquations(SharedSimulation spSim) override;

         /**
          * @brief Add the initial state generation equations
          *
          * @param spGen   Shared generator object
          */
         virtual void addStates(SharedStateGenerator spGen) override;

         /**
          * @brief Add the visualization generation equations
          *
          * @param spGen   Shared visualization generator
          */
         virtual void addVisualizers(SharedVisualizationGenerator spVis) override;

         /**
          * @brief Add the required ASCII output files
          *
          * @param spSim   Shared simulation object
          */
         virtual void addAsciiOutputFiles(SharedSimulation spSim) override;

      protected:
         /**
          * @brieg Register named IDs
          */
         virtual void registerNames() override;

      private:
   };

}
}
}
}
}

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_IRTCMODEL_HPP
