/**
 * @file Momentum.cpp
 * @brief Source of the implementation of the vector Navier-Stokes equation in
 * the Boussinesq rotating thermal convection in a sphere model
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Sphere/RTC/Momentum.hpp"
#include "Model/Boussinesq/Sphere/RTC/MomentumKernel.hpp"
#include "QuICC/Bc/Name/NoSlip.hpp"
#include "QuICC/Bc/Name/StressFree.hpp"
#include "QuICC/NonDimensional/Ekman.hpp"
#include "QuICC/NonDimensional/Rayleigh.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/SolveTiming/Prognostic.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"
#include "QuICC/SpectralKernels/Sphere/ConserveAngularMomentum.hpp"
#include "QuICC/Transform/Path/NoSlipTorPol.hpp"
#include "QuICC/Transform/Path/StressFreeTorPol.hpp"
#include "QuICC/Transform/Path/InsulatingTorPol.hpp"
#include "QuICC/Transform/Path/ValueCurlNl.hpp"
#include "QuICC/Transform/Path/StressFreeCurlNl.hpp"
#include "QuICC/Transform/Path/ValueBc1NegCurlCurlNl.hpp"
#include "QuICC/Transform/Path/InsulatingBc2NegCurlCurlNl.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Sphere {

namespace RTC {

Momentum::Momentum(SharedEquationParameters spEqParams,
   SpatialScheme::SharedCISpatialScheme spScheme,
   std::shared_ptr<Model::IModelBackend> spBackend) :
    IVectorEquation(spEqParams, spScheme, spBackend)
{
   // Set the variable requirements
   this->setRequirements();
}

void Momentum::setCoupling()
{
   int start;
   if (this->ss().has(SpatialScheme::Feature::SpectralOrdering132))
   {
      start = 1;
   }
   else if (this->ss().has(SpatialScheme::Feature::SpectralOrdering123))
   {
      start = 0;
   }
   else
   {
      throw std::logic_error(
         "Unknown spatial scheme was used to setup equations!");
   }

   auto features = defaultCouplingFeature();
   features.at(CouplingFeature::Nonlinear) = true;

   this->defineCoupling(FieldComponents::Spectral::TOR,
      CouplingInformation::PROGNOSTIC, start, features);

   this->defineCoupling(FieldComponents::Spectral::POL,
      CouplingInformation::PROGNOSTIC, start, features);
}

void Momentum::setNLComponents()
{
   std::size_t torPathId;
   std::size_t polPathId;
   auto bcId = this->bcIds().bcId(this->name());
   if(bcId == Bc::Name::NoSlip::id())
   {
#if defined QUICC_BESSEL_VELOCITY_BC_VALUE_TOR_VALUE_POL
      torPathId = Transform::Path::ValueCurlNl::id();
      polPathId = Transform::Path::ValueBc1NegCurlCurlNl::id();
#elif defined QUICC_BESSEL_VELOCITY_BC_VALUE_TOR_INSULATING_POL
      torPathId = Transform::Path::ValueCurlNl::id();
      polPathId = Transform::Path::InsulatingBc2NegCurlCurlNl::id();
#else
#error "Unknown basis setup for Velocity field"
#endif
   }
   else if(bcId == Bc::Name::StressFree::id())
   {
      torPathId = Transform::Path::StressFreeCurlNl::id();
      polPathId = Transform::Path::ValueBc1NegCurlCurlNl::id();
   }
   else
   {
      throw std::logic_error("Unknown Boundary condition");
   }

   this->addNLComponent(FieldComponents::Spectral::TOR, torPathId);

   if (this->couplingInfo(FieldComponents::Spectral::POL).isSplitEquation())
   {
      this->addNLComponent(FieldComponents::Spectral::POL, polPathId);
   }
   else
   {
      this->addNLComponent(FieldComponents::Spectral::POL, polPathId);
   }
}

std::vector<Transform::TransformPath> Momentum::backwardPaths()
{
   std::size_t pathId;
   auto bcId = this->bcIds().bcId(this->name());
   if(bcId == Bc::Name::NoSlip::id())
   {
#if defined QUICC_BESSEL_VELOCITY_BC_VALUE_TOR_VALUE_POL
      pathId = Transform::Path::ValueTorPol::id();
#elif defined QUICC_BESSEL_VELOCITY_BC_VALUE_TOR_INSULATING_POL
      pathId = Transform::Path::InsulatingTorPol::id();
#else
#error "Unknown basis setup for Velocity field"
#endif
   }
   else if(bcId == Bc::Name::StressFree::id())
   {
      pathId = Transform::Path::StressFreeTorPol::id();
   }
   else
   {
      throw std::logic_error("Unknown Boundary condition");
   }

   return this->defaultBackwardPaths(pathId);
}

void Momentum::initNLKernel(const bool force)
{
   // Initialize if empty or forced
   if (force || !this->mspNLKernel)
   {
      // Initialize the physical kernel
      auto spNLKernel = std::make_shared<Physical::Kernel::MomentumKernel>();
      spNLKernel->setVelocity(this->name(), this->spUnknown());
      spNLKernel->setTemperature(PhysicalNames::Temperature::id(),
         this->spScalar(PhysicalNames::Temperature::id()));
      auto T = 1.0 / this->eqParams().nd(NonDimensional::Ekman::id());
      auto Ra = this->eqParams().nd(NonDimensional::Rayleigh::id());
      spNLKernel->init(1.0, T, Ra * T);
      this->mspNLKernel = spNLKernel;
   }
}

void Momentum::initConstraintKernel(const std::shared_ptr<std::vector<Array>>)
{
   if (this->bcIds().bcId(this->name()) == Bc::Name::StressFree::id())
   {
      // Initialize the physical kernel
      auto spConstraint =
         std::make_shared<Spectral::Kernel::Sphere::ConserveAngularMomentum>(
            this->ss().has(SpatialScheme::Feature::ComplexSpectrum));
      spConstraint->setField(this->name(), this->spUnknown());
      spConstraint->setResolution(this->spRes());
      spConstraint->init(
         this->ss().has(SpatialScheme::Feature::SpectralOrdering123));
      this->setConstraintKernel(FieldComponents::Spectral::TOR, spConstraint);
   }
}

void Momentum::setRequirements()
{
   // Set velocity as equation unknown
   this->setName(PhysicalNames::Velocity::id());

   // Set solver timing
   this->setSolveTiming(SolveTiming::Prognostic::id());

   // Forward transform generates nonlinear RHS
   this->setForwardPathsType(FWD_IS_NONLINEAR);

   // Get reference to spatial scheme
   const auto& ss = this->ss();

   // Add velocity to requirements: is scalar?
   auto& velReq = this->mRequirements.addField(PhysicalNames::Velocity::id(),
      FieldRequirement(false, ss.spectral(), ss.physical()));
   velReq.enableSpectral();
   velReq.enablePhysical();
   velReq.enableCurl();

   // Add temperature to requirements: is scalar?
   auto& tempReq =
      this->mRequirements.addField(PhysicalNames::Temperature::id(),
         FieldRequirement(true, ss.spectral(), ss.physical()));
   tempReq.enableSpectral();
   tempReq.enablePhysical();
}

} // namespace RTC
} // namespace Sphere
} // namespace Boussinesq
} // namespace Equations
} // namespace QuICC
