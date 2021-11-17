/**
 * @file Momentum.hpp
 * @brief Implementation of the vector Navier-Stokes equation for the Boussinesq rotating thermal convection in a sphere
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_MOMENTUM_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_MOMENTUM_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "QuICC/Typedefs.hpp"
#include "QuICC/Equations/IVectorEquation.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Sphere {

namespace RTC {

   /**
    * @brief Implementation of the vector Navier-Stokes equation for the Boussinesq rotating thermal convection in a sphere
    */
   class Momentum: public IVectorEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         Momentum(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple empty destructor
          */
         virtual ~Momentum();

         /**
          * @brief Initialize constraint kernel
          */
         virtual void initConstraintKernel();

         /**
          * @brief Initialize nonlinear interaction kernel
          */
         virtual void initNLKernel(const bool force = false);

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements();

         /**
          * @brief Set the equation coupling information
          */
         virtual void setCoupling();

         /**
          * @brief Set the nonlinear integration components
          */
         virtual void setNLComponents();

      private:
   };

}
}
}
}
}

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_MOMENTUM_HPP
