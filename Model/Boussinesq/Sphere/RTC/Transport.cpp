/**
 * @file Transport.cpp
 * @brief Source of the implementation of the transport equation in the
 * Boussinesq rotating thermal convection in a sphere
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Sphere/RTC/Transport.hpp"
#include "Model/Boussinesq/Sphere/RTC/TransportKernel.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/SolveTiming/Prognostic.hpp"
#include "QuICC/Transform/Path/I2ScalarNl.hpp"
#include "QuICC/Transform/Path/Scalar.hpp"

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
   features.at(CouplingFeature::BoundaryValue) = true;

   this->defineCoupling(FieldComponents::Spectral::SCALAR,
      CouplingInformation::PROGNOSTIC, 0, features);
}

void Transport::setNLComponents()
{
   this->addNLComponent(FieldComponents::Spectral::SCALAR,
      Transform::Path::I2ScalarNl::id());
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

MHDVariant Transport::boundaryValue(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
{
   assert(compId == FieldComponents::Spectral::SCALAR);

   MHDComplex val = 0;

   // Boundary condition is applied to first row
   if(i == 0)
   {
      int l, m;
      const auto& tRes = *this->spRes()->cpu()->dim(Dimensions::Transform::SPECTRAL);
      if(this->spRes()->sim().ss().has(SpatialScheme::Feature::SpectralOrdering132))
      {
         l = tRes.idx<Dimensions::Data::DAT3D>(k);
         m = tRes.idx<Dimensions::Data::DAT2D>(j,k);
      } else
      {
         l = tRes.idx<Dimensions::Data::DAT2D>(j,k);
         m = tRes.idx<Dimensions::Data::DAT3D>(k);
      }

      // Set value at boundary for l,m mode
      if(l == 3 && m == 3)
      {
         val = MHDComplex(42,42); // 42 + 42j
      }
   }

   return val;
}

} // namespace RTC
} // namespace Sphere
} // namespace Boussinesq
} // namespace Equations
} // namespace QuICC
