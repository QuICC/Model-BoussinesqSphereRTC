/**
 * @file Transport.hpp
 * @brief Implementation of the transport equation for the Boussinesq rotating
 * thermal convection in a sphere
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_TRANSPORT_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_TRANSPORT_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Equations/IScalarEquation.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Sphere {

namespace RTC {

/**
 * @brief Implementation of the transport equation for the Boussinesq rotating
 * thermal convection in a sphere
 */
class Transport : public IScalarEquation
{
public:
   /**
    * @brief Simple constructor
    *
    * @param spEqParams  Shared equation parameters
    */
   Transport(SharedEquationParameters spEqParams,
      SpatialScheme::SharedCISpatialScheme spScheme,
      std::shared_ptr<Model::IModelBackend> spBackend);

   /**
    * @brief Simple empty destructor
    */
   ~Transport() = default;

   /**
    * @brief Initialize nonlinear interaction kernel
    */
   void initNLKernel(const bool force = false) final;

   /**
    * @brief Generic boundary value implementation
    *
    * @param compId  ID of the spectral component
    * @param i       Fastest index
    * @param j       Second index
    * @param k       Slowest index
    */
   virtual MHDVariant boundaryValue(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const;

protected:
   /**
    * @brief Set variable requirements
    */
   void setNLComponents() final;

   /**
    * @brief Set variable requirements
    */
   void setRequirements() final;

   /**
    * @brief Set the equation coupling information
    */
   void setCoupling() final;

private:
};

} // namespace RTC
} // namespace Sphere
} // namespace Boussinesq
} // namespace Equations
} // namespace QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_TRANSPORT_HPP
