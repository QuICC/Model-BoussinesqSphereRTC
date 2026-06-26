/**
 * @file IExponentialModel.cpp
 * @brief Source of the Boussinesq rotating thermal convection in a sphere
 * (Toroidal/Poloidal formulation)
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Sphere/RTC/Exponential/IExponentialModel.hpp"
#include "Model/Boussinesq/Sphere/RTC/Momentum.hpp"
#include "Model/Boussinesq/Sphere/RTC/Transport.hpp"
#include "Model/Boussinesq/Sphere/RTC/Exponential/MomentumJacobian.hpp"
#include "Model/Boussinesq/Sphere/RTC/Exponential/TransportJacobian.hpp"
#include "Model/Boussinesq/Sphere/RTC/gitHash.hpp"
#include "QuICC/PhysicalNames/JacobianTemperature.hpp"
#include "QuICC/PhysicalNames/JacobianVelocity.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace RTC {

namespace Exponential {

std::vector<std::size_t> IExponentialModel::excludedFieldIds() const
{
   std::vector<std::size_t> fields = {
      PhysicalNames::JacobianVelocity::id(),
      PhysicalNames::JacobianTemperature::id(),
   };

   return fields;
}

void IExponentialModel::addEquations(SharedSimulation spSim)
{
   auto optZero = std::make_shared<Equations::EquationOptions>(0, false, false);

   // Add transport equation
   spSim->addEquation<Equations::Boussinesq::Sphere::RTC::Transport>(
      this->spBackend(), optZero);

   // Add Navier-Stokes equation
   spSim->addEquation<Equations::Boussinesq::Sphere::RTC::Momentum>(
      this->spBackend(), optZero);

   auto optOne = std::make_shared<Equations::EquationOptions>(1, false, false, false);

   // Add transport jacobian equation
   spSim->addEquation<Equations::Boussinesq::Sphere::RTC::Exponential::TransportJacobian>(
      this->spBackend(), optOne);

   // Add Navier-Stokes jacobian equation
   spSim->addEquation<Equations::Boussinesq::Sphere::RTC::Exponential::MomentumJacobian>(
      this->spBackend(), optOne);

   #ifdef QUICC_USE_MLIR_GRAPH
      #error "MLIR Graph setup not implemented for InertialessDynamo"
   #endif
}

} // namespace Exponential
} // namespace RTC
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
