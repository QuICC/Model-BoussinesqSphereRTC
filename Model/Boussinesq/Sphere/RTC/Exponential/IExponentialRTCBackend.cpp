/**
 * @file IExponentialRTCBackend.cpp
 * @brief Source of the interface for model backend
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "Model/Boussinesq/Sphere/RTC/Exponential/IExponentialRTCBackend.hpp"
#include "QuICC/Bc/Name/FixedFlux.hpp"
#include "QuICC/Bc/Name/FixedTemperature.hpp"
#include "QuICC/Bc/Name/NoSlip.hpp"
#include "QuICC/Bc/Name/StressFree.hpp"
#include "QuICC/Bc/Name/QuasiInverseOnly.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/PhysicalNames/JacobianVelocity.hpp"
#include "QuICC/PhysicalNames/JacobianTemperature.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/Polynomial/Worland/WorlandTypes.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/SparseSM/Worland/Boundary/D1.hpp"
#include "QuICC/SparseSM/Worland/Boundary/D2.hpp"
#include "QuICC/SparseSM/Worland/Boundary/Operator.hpp"
#include "QuICC/SparseSM/Worland/Boundary/R1D1DivR1.hpp"
#include "QuICC/SparseSM/Worland/Boundary/Value.hpp"
#include "QuICC/SparseSM/Worland/Id.hpp"
#include "QuICC/SparseSM/Worland/Stencil/D1.hpp"
#include "QuICC/SparseSM/Worland/Stencil/R1D1DivR1.hpp"
#include "QuICC/SparseSM/Worland/Stencil/Value.hpp"
#include "QuICC/SparseSM/Worland/Stencil/ValueD1.hpp"
#include "QuICC/SparseSM/Worland/Stencil/ValueD2.hpp"
#include "DenseSM/Worland/Stencil/OrthogonalValue.hpp"
#include "DenseSM/Worland/Stencil/OrthogonalValueD1.hpp"
#include "DenseSM/Worland/Stencil/OrthogonalValueD2.hpp"

//#define QUICC_USE_ORTHOGONAL_STENCIL

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace RTC {

namespace Exponential {

std::vector<std::string> IExponentialRTCBackend::fieldNames() const
{
   std::vector<std::string> names = {
      PhysicalNames::Velocity().tag(),
      PhysicalNames::Temperature().tag(),
      PhysicalNames::JacobianVelocity().tag(),
      PhysicalNames::JacobianTemperature().tag()};

   return names;
}

int IExponentialRTCBackend::nBc(const SpectralFieldId& fId) const
{
   int nBc = 0;

   auto vel_tor = std::make_pair(PhysicalNames::Velocity::id(),
                   FieldComponents::Spectral::TOR);
   auto jvel_tor = std::make_pair(PhysicalNames::JacobianVelocity::id(),
                   FieldComponents::Spectral::TOR);
   auto vel_pol = std::make_pair(PhysicalNames::Velocity::id(),
                   FieldComponents::Spectral::POL);
   auto jvel_pol = std::make_pair(PhysicalNames::JacobianVelocity::id(),
                   FieldComponents::Spectral::POL);
   auto temp = std::make_pair(PhysicalNames::Temperature::id(),
                   FieldComponents::Spectral::SCALAR);
   auto jtemp = std::make_pair(PhysicalNames::JacobianTemperature::id(),
                   FieldComponents::Spectral::SCALAR);

   if (fId == vel_tor ||
       fId == jvel_tor ||
       fId == temp ||
       fId == jtemp)
   {
      nBc = 1;
   }
   else if (fId == vel_pol ||
            fId == jvel_pol)
   {
      nBc = 2;
   }
   else
   {
      nBc = 0;
   }

   return nBc;
}

void IExponentialRTCBackend::applyTau(SparseMatrix& mat, const SpectralFieldId& rowId,
   const SpectralFieldId& colId, const int l,
   std::shared_ptr<details::BlockOptions> opts, const int nN,
   const BcMap& bcs, const NonDimensional::NdMap& nds,
   const bool isSplitOperator) const
{
   auto a = Polynomial::Worland::worland_default_t::ALPHA;
   auto b = Polynomial::Worland::worland_default_t::DBETA;

   auto bcId = bcs.find(rowId.first)->second;

   SparseSM::Worland::Boundary::Operator bcOp(nN, nN, a, b, l);

   auto vel_tor = std::make_pair(PhysicalNames::Velocity::id(),
                   FieldComponents::Spectral::TOR);
   auto jvel_tor = std::make_pair(PhysicalNames::JacobianVelocity::id(),
                   FieldComponents::Spectral::TOR);
   auto vel_pol = std::make_pair(PhysicalNames::Velocity::id(),
                   FieldComponents::Spectral::POL);
   auto jvel_pol = std::make_pair(PhysicalNames::JacobianVelocity::id(),
                   FieldComponents::Spectral::POL);
   auto temp = std::make_pair(PhysicalNames::Temperature::id(),
                   FieldComponents::Spectral::SCALAR);
   auto jtemp = std::make_pair(PhysicalNames::JacobianTemperature::id(),
                   FieldComponents::Spectral::SCALAR);

   if ((rowId == vel_tor || rowId == jvel_tor) &&
       rowId == colId)
   {
      if (l > 0)
      {
         if (bcId == Bc::Name::NoSlip::id())
         {
            bcOp.addRow<SparseSM::Worland::Boundary::Value>();
         }
         else if (bcId == Bc::Name::StressFree::id())
         {
            bcOp.addRow<SparseSM::Worland::Boundary::R1D1DivR1>();
         }
         else
         {
            throw std::logic_error("Boundary conditions for Velocity "
                                   "Toroidal component not implemented");
         }
      }
   }
   else if ((rowId == vel_pol || rowId == jvel_pol) &&
            rowId == colId)
   {
      if (l > 0)
      {
         if (this->useSplitEquation())
         {
            if (isSplitOperator)
            {
               bcOp.addRow<SparseSM::Worland::Boundary::Value>();
            }
            else if (bcId == Bc::Name::NoSlip::id())
            {
               bcOp.addRow<SparseSM::Worland::Boundary::D1>();
            }
            else if (bcId == Bc::Name::StressFree::id())
            {
               bcOp.addRow<SparseSM::Worland::Boundary::D2>();
            }
            else
            {
               throw std::logic_error(
                  "Boundary conditions for Velocity Poloidal component "
                  "not implemented");
            }
         }
         else
         {
            if (bcId == Bc::Name::NoSlip::id())
            {
               bcOp.addRow<SparseSM::Worland::Boundary::Value>();
               bcOp.addRow<SparseSM::Worland::Boundary::D1>();
            }
            else if (bcId == Bc::Name::StressFree::id())
            {
               bcOp.addRow<SparseSM::Worland::Boundary::Value>();
               bcOp.addRow<SparseSM::Worland::Boundary::D2>();
            }
            else
            {
               throw std::logic_error(
                  "Boundary conditions for Velocity Poloidal component "
                  "not implemented");
            }
         }
      }
   }
   else if ((rowId == temp || rowId == jtemp) &&
            rowId == colId)
   {
      if (bcId == Bc::Name::FixedTemperature::id())
      {
         bcOp.addRow<SparseSM::Worland::Boundary::Value>();
      }
      else if (bcId == Bc::Name::FixedFlux::id())
      {
         bcOp.addRow<SparseSM::Worland::Boundary::D1>();
      }
      else
      {
         throw std::logic_error(
            "Boundary conditions for Temperature not implemented (" +
            std::to_string(bcId) + ")");
      }
   }

   mat.real() += bcOp.mat();
}

void IExponentialRTCBackend::stencil(SparseMatrix& mat, const SpectralFieldId& fieldId,
   const int l, const int nN, const bool makeSquare, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   auto a = Polynomial::Worland::worland_default_t::ALPHA;
   auto b = Polynomial::Worland::worland_default_t::DBETA;

   auto bcId = bcs.find(fieldId.first)->second;

   auto vel_tor = std::make_pair(PhysicalNames::Velocity::id(),
                   FieldComponents::Spectral::TOR);
   auto jvel_tor = std::make_pair(PhysicalNames::JacobianVelocity::id(),
                   FieldComponents::Spectral::TOR);
   auto vel_pol = std::make_pair(PhysicalNames::Velocity::id(),
                   FieldComponents::Spectral::POL);
   auto jvel_pol = std::make_pair(PhysicalNames::JacobianVelocity::id(),
                   FieldComponents::Spectral::POL);
   auto temp = std::make_pair(PhysicalNames::Temperature::id(),
                   FieldComponents::Spectral::SCALAR);
   auto jtemp = std::make_pair(PhysicalNames::JacobianTemperature::id(),
                   FieldComponents::Spectral::SCALAR);

   int s = this->nBc(fieldId);
   if(bcId == Bc::Name::QuasiInverseOnly::id())
   {
      SparseSM::Worland::Id qid(nN, nN - s, a, b, l);
      mat = qid.mat();
   }
   else
   {
      if (fieldId == vel_tor || fieldId == jvel_tor)
      {
         if (bcId == Bc::Name::NoSlip::id())
         {
#ifdef QUICC_USE_ORTHOGONAL_STENCIL
            DenseSM::Worland::Stencil::OrthogonalValue bc(nN, nN - s, a, b, l);
            mat = bc.spmat();
#else
            SparseSM::Worland::Stencil::Value bc(nN, nN - s, a, b, l);
            mat = bc.mat();
#endif
         }
         else if (bcId == Bc::Name::StressFree::id())
         {
#ifdef QUICC_USE_ORTHOGONAL_STENCIL
            throw std::logic_error("Orthogonal stencil for toroidal stress-free boundary not yet implemented");
#else
            SparseSM::Worland::Stencil::R1D1DivR1 bc(nN, nN - s, a, b, l);
            mat = bc.mat();
#endif
         }
         else
         {
            throw std::logic_error("Galerkin boundary conditions for Velocity "
                  "Toroidal component not implemented");
         }
      }
      else if (fieldId == vel_pol || fieldId == jvel_pol)
      {
         if (bcId == Bc::Name::NoSlip::id())
         {
#ifdef QUICC_USE_ORTHOGONAL_STENCIL
            DenseSM::Worland::Stencil::OrthogonalValueD1 bc(nN, nN - s, a, b, l);
            mat = bc.spmat();
#else
            SparseSM::Worland::Stencil::ValueD1 bc(nN, nN - s, a, b, l);
            mat = bc.mat();
#endif
         }
         else if (bcId == Bc::Name::StressFree::id())
         {
#ifdef QUICC_USE_ORTHOGONAL_STENCIL
            DenseSM::Worland::Stencil::OrthogonalValueD2 bc(nN, nN - s, a, b, l);
            mat = bc.spmat();
#else
            SparseSM::Worland::Stencil::ValueD2 bc(nN, nN - s, a, b, l);
            mat = bc.mat();
#endif
         }
         else
         {
            throw std::logic_error("Galerin boundary conditions for Velocity "
                  "Poloidal component not implemented");
         }
      }
      else if (fieldId == temp || fieldId == jtemp)
      {
         if (bcId == Bc::Name::FixedTemperature::id())
         {
#ifdef QUICC_USE_ORTHOGONAL_STENCIL
            DenseSM::Worland::Stencil::OrthogonalValue bc(nN, nN - s, a, b, l);
            mat = bc.spmat();
#else
            SparseSM::Worland::Stencil::Value bc(nN, nN - s, a, b, l);
            mat = bc.mat();
#endif
         }
         else if (bcId == Bc::Name::FixedFlux::id())
         {
#ifdef QUICC_USE_ORTHOGONAL_STENCIL
            throw std::logic_error("Orthogonal stencil for fixed-flux boundary not yet implemented");
#else
            SparseSM::Worland::Stencil::D1 bc(nN, nN - s, a, b, l);
            mat = bc.mat();
#endif
         }
         else
         {
            throw std::logic_error(
                  "Galerkin boundary conditions for Temperature not implemented");
         }
      }
   }

   if (makeSquare)
   {
      SparseSM::Worland::Id qId(nN - s, nN, a, b, l);
      mat = qId.mat() * mat;
   }
}

} // namespace Exponential
} // namespace RTC
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
