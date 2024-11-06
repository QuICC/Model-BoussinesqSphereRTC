/**
 * @file Transport.cpp
 * @brief Source of the implementation of the transport equation in the
 * Boussinesq rotating thermal convection in a sphere
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Sphere/RTC/Utils.hpp"
#include "Model/Boussinesq/Sphere/RTC/Transport.hpp"
#include "Model/Boussinesq/Sphere/RTC/TransportKernel.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/SolveTiming/Prognostic.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Sphere {

namespace RTC {

Transport::Transport(SharedEquationParameters spEqParams,
   SpatialScheme::SharedCISpatialScheme spScheme,
   std::shared_ptr<Model::IModelBackend> spBackend) :
    IScalarEquation(spEqParams, spScheme, spBackend)
{
   // Set the variable requirements
   this->setRequirements();
}

void Transport::setCoupling()
{
   auto features = defaultCouplingFeature();
   features.at(CouplingFeature::Nonlinear) = true;

   this->defineCoupling(FieldComponents::Spectral::SCALAR,
      CouplingInformation::PROGNOSTIC, 0, features);
}

void Transport::setNLComponents()
{
   using Model::Boussinesq::Sphere::RTC::getPathId;
   std::size_t pathId = getPathId(this->mspBcIds, this->name(), true, FieldComponents::Spectral::SCALAR);

   this->addNLComponent(FieldComponents::Spectral::SCALAR,
      pathId);
}

std::vector<Transform::TransformPath> Transport::backwardPaths()
{
   using Model::Boussinesq::Sphere::RTC::getPathId;
   std::size_t pathId = getPathId(this->mspBcIds, this->name(), false);

   return this->defaultBackwardPaths(pathId);
}

void Transport::initNLKernel(const bool force)
{
   // Initialize if empty or forced
   if (force || !this->mspNLKernel)
   {
      // Initialize the physical kernel
      auto spNLKernel = std::make_shared<Physical::Kernel::TransportKernel>();
      spNLKernel->setScalar(this->name(), this->spUnknown());
      spNLKernel->setVector(PhysicalNames::Velocity::id(),
         this->spVector(PhysicalNames::Velocity::id()));
      spNLKernel->init(1.0);
      this->mspNLKernel = spNLKernel;
   }
}

void Transport::setRequirements()
{
   // Set temperatur as equation unknown
   this->setName(PhysicalNames::Temperature::id());

   // Set solver timing
   this->setSolveTiming(SolveTiming::Prognostic::id());

   // Forward transform generates nonlinear RHS
   this->setForwardPathsType(FWD_IS_NONLINEAR);

   // Get reference to spatial scheme
   const auto& ss = this->ss();

   // Add temperature to requirements: is scalar?, need spectral?, need
   // physical?, need diff?
   auto& tempReq =
      this->mRequirements.addField(PhysicalNames::Temperature::id(),
         FieldRequirement(true, ss.spectral(), ss.physical()));
   tempReq.enableSpectral();
   tempReq.enableGradient();

   // Add velocity to requirements: is scalar?, need spectral?, need physical?,
   // need diff?(, need curl?)
   auto& velReq = this->mRequirements.addField(PhysicalNames::Velocity::id(),
      FieldRequirement(false, ss.spectral(), ss.physical()));
   velReq.enableSpectral();
   velReq.enablePhysical();
}

} // namespace RTC
} // namespace Sphere
} // namespace Boussinesq
} // namespace Equations
} // namespace QuICC
