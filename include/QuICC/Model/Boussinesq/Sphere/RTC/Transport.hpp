/**
 * @file Transport.hpp
 * @brief Implementation of the transport equation for the Boussinesq rotating thermal convection in a sphere
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_TRANSPORT_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_TRANSPORT_HPP

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
#include "QuICC/Equations/IScalarEquation.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Sphere {

namespace RTC {

   /**
    * @brief Implementation of the transport equation for the Boussinesq rotating thermal convection in a sphere
    */
   class Transport: public IScalarEquation
   {
      public:
         /**
          * @brief Simple constructor
          *
          * @param spEqParams  Shared equation parameters
          */
         Transport(SharedEquationParameters spEqParams, SpatialScheme::SharedCISpatialScheme spScheme, std::shared_ptr<Model::IModelBackend> spBackend);

         /**
          * @brief Simple empty destructor
          */
         virtual ~Transport();

         /**
          * @brief Initialize nonlinear interaction kernel
          */
         virtual void initNLKernel(const bool force = false) override;

         /**
          * @brief Generic boundary value implementation
          *
          * @param compId  ID of the spectral component
          * @param i       Fastest index
          * @param j       Second index
          * @param k       Slowest index
          */
         MHDVariant boundaryValue(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const final;

      protected:
         /**
          * @brief Set variable requirements
          */
         virtual void setNLComponents() override;

         /**
          * @brief Set variable requirements
          */
         virtual void setRequirements() override;

         /**
          * @brief Set the equation coupling information
          */
         virtual void setCoupling() override;

      private:
   };

}
}
}
}
}

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_TRANSPORT_HPP
