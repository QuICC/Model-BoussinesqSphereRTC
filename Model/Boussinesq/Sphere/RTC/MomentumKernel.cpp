/**
 * @file MomentumKernel.cpp
 * @brief Source of physical space kernel for the Momentum equation
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Sphere/RTC/MomentumKernel.hpp"
#include "QuICC/PhysicalOperators/Cross.hpp"
#include "QuICC/PhysicalOperators/SphericalBuoyancy.hpp"
#include "QuICC/PhysicalOperators/SphericalCoriolis.hpp"
#include "QuICC/SpatialScheme/ISpatialScheme.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

std::size_t MomentumKernel::name() const
{
   return this->mName;
}

void MomentumKernel::setVelocity(std::size_t name,
   Framework::Selector::VariantSharedVectorVariable spField)
{
   // Safety assertion
   assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

   this->mName = name;

   this->setField(name, spField);
}

void MomentumKernel::setTemperature(std::size_t name,
   Framework::Selector::VariantSharedScalarVariable spField)
{
   // Safety assertion
   assert(this->mScalars.count(name) + this->mVectors.count(name) == 0);

   this->mTempName = name;

   this->setField(name, spField);
}

void MomentumKernel::init(const MHDFloat inertia, const MHDFloat coriolis,
   const MHDFloat buoyancy)
{
   // Set scaling constants
   this->mInertia = inertia;
   this->mCoriolis = coriolis;
   this->mBuoyancy = buoyancy;
}

void MomentumKernel::setMesh(std::shared_ptr<std::vector<Array>> spMesh)
{
   IPhysicalKernel::setMesh(spMesh);

   this->mRadius = this->mspMesh->at(0);
   if (std::visit(
          [&](auto&& v) -> bool
          {
             return (v->dom(0).res().sim().ss().has(
                SpatialScheme::Feature::SpectralOrdering132));
          },
          this->vector(this->name())))
   {
      this->mCosTheta = this->mspMesh->at(1).array().cos();
      this->mSinTheta = this->mspMesh->at(1).array().sin();
   }
}

void MomentumKernel::compute(Framework::Selector::PhysicalScalarField& rNLComp,
   FieldComponents::Physical::Id id) const
{
   ///
   /// Compute \f$\left(\nabla\wedge\vec u\right)\wedge\vec u\f$
   ///
   std::visit(
      [&](auto&& v)
      {
         switch (id)
         {
         case (FieldComponents::Physical::R):
            Physical::Cross<FieldComponents::Physical::THETA,
               FieldComponents::Physical::PHI>::set(rNLComp, v->dom(0).curl(),
               v->dom(0).phys(), this->mInertia);
            break;
         case (FieldComponents::Physical::THETA):
            Physical::Cross<FieldComponents::Physical::PHI,
               FieldComponents::Physical::R>::set(rNLComp, v->dom(0).curl(),
               v->dom(0).phys(), this->mInertia);
            break;
         case (FieldComponents::Physical::PHI):
            Physical::Cross<FieldComponents::Physical::R,
               FieldComponents::Physical::THETA>::set(rNLComp, v->dom(0).curl(),
               v->dom(0).phys(), this->mInertia);
            break;
         default:
            assert(false);
            break;
         }
      },
      this->vector(this->name()));

   std::visit(
      [&](auto&& s)
      {
         Physical::SphericalBuoyancy::sub(rNLComp, id, s->dom(0).res(),
            this->mRadius, s->dom(0).phys(), this->mBuoyancy);
      },
      this->scalar(this->mTempName));

   std::visit(
      [&](auto&& v)
      {
         if (v->dom(0).res().sim().ss().has(
                SpatialScheme::Feature::SpectralOrdering132))
         {
            ///
            /// Compute Coriolis term
            ///
            Physical::SphericalCoriolis::add(rNLComp, id, v->dom(0).res(),
               this->mCosTheta, this->mSinTheta, v->dom(0).phys(),
               this->mCoriolis);
         }
      },
      this->vector(this->name()));
}

} // namespace Kernel
} // namespace Physical
} // namespace QuICC
