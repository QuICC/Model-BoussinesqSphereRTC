/**
 * @file Utils.cpp
 * @brief Source of the implementation utilities
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Sphere/RTC/Utils.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/Transform/Path/JlK0ScalarNl.hpp"
#include "QuICC/Transform/Path/JlK0Scalar.hpp"
#include "QuICC/Transform/Path/JlK2ScalarNl.hpp"
#include "QuICC/Transform/Path/JlK2Scalar.hpp"
#include "QuICC/Transform/Path/JlK0JlK2TorPol.hpp"
#include "QuICC/Transform/Path/JlK0TorPol.hpp"
#include "QuICC/Transform/Path/JlK2TorPol.hpp"
#include "QuICC/Transform/Path/JlK0CurlNl.hpp"
#include "QuICC/Transform/Path/JlK2CurlNl.hpp"
#include "QuICC/Transform/Path/JlK0NegCurlCurlNl.hpp"
#include "QuICC/Transform/Path/JlK2NegCurlCurlNl.hpp"
#include "QuICC/Bc/Name/FixedTemperature.hpp"
#include "QuICC/Bc/Name/FixedFlux.hpp"
#include "QuICC/Bc/Name/NoSlip.hpp"
#include "QuICC/Bc/Name/StressFree.hpp"

//#define TEMPERATURE_USE_JLK2_BASIS
//#define VELOCITY_TOR_USE_JLK2_BASIS
//#define VELOCITY_POL_USE_JLK2_BASIS

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace RTC {

std::size_t getPathId(std::shared_ptr<SimulationBoundary> spBcs, const std::size_t fieldId, const bool isNl, const FieldComponents::Spectral::Id comp)
{
   std::size_t pathId;

   // Temperature
   if(fieldId == PhysicalNames::Temperature::id())
   {
      if(spBcs->bcId(fieldId) == Bc::Name::FixedTemperature::id())
      {
         if(isNl)
         {
#ifdef TEMPERATURE_USE_JLK2_BASIS
            pathId = Transform::Path::JlK2ScalarNl::id();
#else
            pathId = Transform::Path::JlK0ScalarNl::id();
#endif
         }
         else
         {
#ifdef TEMPERATURE_USE_JLK2_BASIS
            pathId = Transform::Path::JlK2Scalar::id();
#else
            pathId = Transform::Path::JlK0Scalar::id();
#endif
         }
      }
      else
      {
         throw std::logic_error("Boundary condition for Temperature not implemented");
      }
   }
   // Velocity
   else if(fieldId == PhysicalNames::Velocity::id())
   {
      if(spBcs->bcId(fieldId) == Bc::Name::NoSlip::id())
      {
         if(isNl)
         {
            if(comp == FieldComponents::Spectral::TOR)
            {
#ifdef VELOCITY_TOR_USE_JLK2_BASIS
               pathId = Transform::Path::JlK2CurlNl::id();
#else
               pathId = Transform::Path::JlK0CurlNl::id();
#endif
            }
            else
            {
#ifdef VELOCITY_POL_USE_JLK2_BASIS
               pathId = Transform::Path::JlK2NegCurlCurlNl::id();
#else
               pathId = Transform::Path::JlK0NegCurlCurlNl::id();
#endif
            }
         }
         else
         {
#if defined(VELOCITY_TOR_USE_JLK2_BASIS) && defined(VELOCITY_POL_USE_JLK2_BASIS)
            pathId = Transform::Path::JlK2TorPol::id();
#elif defined(VELOCITY_POL_USE_JLK2_BASIS)
            pathId = Transform::Path::JlK0JlK2TorPol::id();
#else
            pathId = Transform::Path::JlK0TorPol::id();
#endif
         }
      }
      else if(spBcs->bcId(fieldId) == Bc::Name::StressFree::id())
      {
         throw std::logic_error("Stress-free boundary condition path are not implemented");
      }
   }

   return pathId;
}

SparseSM::Bessel::BesselKind bKind(const SpectralFieldId& fId)
{
   if (fId == std::make_pair(PhysicalNames::Velocity::id(),
                 FieldComponents::Spectral::TOR))
   {
#ifdef VELOCITY_TOR_USE_JLK2_BASIS
      return SparseSM::Bessel::BesselKind::JlK2;
#else
      return SparseSM::Bessel::BesselKind::JlK0;
#endif
   }
   else if (fId == std::make_pair(PhysicalNames::Velocity::id(),
                 FieldComponents::Spectral::POL))
   {
#ifdef VELOCITY_POL_USE_JLK2_BASIS
      return SparseSM::Bessel::BesselKind::JlK2;
#else
      return SparseSM::Bessel::BesselKind::JlK0;
#endif
   }
   else if (fId == std::make_pair(PhysicalNames::Temperature::id(),
                 FieldComponents::Spectral::SCALAR))
   {
#ifdef TEMPERATURE_USE_JLK2_BASIS
      return SparseSM::Bessel::BesselKind::JlK2;
#else
      return SparseSM::Bessel::BesselKind::JlK0;
#endif
   }
   else
   {
      throw std::logic_error("Unknown spectral field");
   }
}

Polynomial::SphericalBessel::sphjnl_basis_t bBasis(const SparseSM::Bessel::BesselKind kind)
{
   Polynomial::SphericalBessel::sphjnl_basis_t basis;
   if(kind == SparseSM::Bessel::BesselKind::JlK0)
   {
      basis = Polynomial::SphericalBessel::jl_t();
   }
   else if(kind == SparseSM::Bessel::BesselKind::Jlm1K0)
   {
      basis = Polynomial::SphericalBessel::jlm1_t();
   }
   else if(kind == SparseSM::Bessel::BesselKind::JlK2)
   {
      basis = Polynomial::SphericalBessel::jl_k2_t();
   }
   else
   {
      throw std::logic_error("Unknown Bessel kind");
   }

   return basis;
}

} // namespace RTC
} // namespace Sphere
} // namespace Boussinesq
} // namespace Equations
} // namespace QuICC
