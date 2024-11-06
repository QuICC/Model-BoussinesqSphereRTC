/**
 * @file IRTCBackend.cpp
 * @brief Source of the interface for model backend
 */

// System includes
//
#include <exception>
#include <functional>
#include <stdexcept>

// Project includes
//
#include "Model/Boussinesq/Sphere/RTC/IRTCBackend.hpp"
#include "Model/Boussinesq/Sphere/RTC/Utils.hpp"
#include "QuICC/Bc/Name/FixedFlux.hpp"
#include "QuICC/Bc/Name/FixedTemperature.hpp"
#include "QuICC/Bc/Name/NoSlip.hpp"
#include "QuICC/Bc/Name/StressFree.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/ModelOperator/Boundary.hpp"
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
#include "QuICC/ModelOperator/ImplicitLinear.hpp"
#include "QuICC/ModelOperator/SplitBoundary.hpp"
#include "QuICC/ModelOperator/SplitImplicitLinear.hpp"
#include "QuICC/ModelOperator/Stencil.hpp"
#include "QuICC/ModelOperator/Time.hpp"
#include "QuICC/ModelOperatorBoundary/FieldToRhs.hpp"
#include "QuICC/ModelOperatorBoundary/SolverHasBc.hpp"
#include "QuICC/ModelOperatorBoundary/SolverNoTau.hpp"
#include "QuICC/ModelOperatorBoundary/Stencil.hpp"
#include "QuICC/NonDimensional/CflInertial.hpp"
#include "QuICC/NonDimensional/Ekman.hpp"
#include "QuICC/NonDimensional/Prandtl.hpp"
#include "QuICC/NonDimensional/Rayleigh.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/SparseSM/Bessel/BesselKind.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#include "QuICC/SparseSM/Id.hpp"
#include "QuICC/SparseSM/Bessel/Boundary/Operator.hpp"
#include "QuICC/SparseSM/Bessel/Boundary/Value.hpp"
#include "QuICC/SparseSM/Bessel/Boundary/D1.hpp"
#include "QuICC/SparseSM/Bessel/Boundary/D2.hpp"
#include "QuICC/SparseSM/Bessel/Boundary/R1D1DivR1.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace RTC {

std::vector<std::string> IRTCBackend::fieldNames() const
{
   std::vector<std::string> names = {PhysicalNames::Velocity().tag(),
      PhysicalNames::Temperature().tag()};

   return names;
}

std::vector<std::string> IRTCBackend::paramNames() const
{
   std::vector<std::string> names = {NonDimensional::Prandtl().tag(),
      NonDimensional::Rayleigh().tag(), NonDimensional::Ekman().tag()};

   return names;
}

std::vector<bool> IRTCBackend::isPeriodicBox() const
{
   std::vector<bool> periodic = {false, false, false};

   return periodic;
}

std::map<std::string, MHDFloat> IRTCBackend::automaticParameters(
   const std::map<std::string, MHDFloat>& cfg) const
{
   auto E = cfg.find(NonDimensional::Ekman().tag())->second;

   std::map<std::string, MHDFloat> params = {
      {NonDimensional::CflInertial().tag(), 0.1 * E}};

   return params;
}

int IRTCBackend::nBc(const SpectralFieldId& fId) const
{
   int nBc = 0;

   if (fId == std::make_pair(PhysicalNames::Velocity::id(),
                 FieldComponents::Spectral::TOR) ||
       fId == std::make_pair(PhysicalNames::Temperature::id(),
                 FieldComponents::Spectral::SCALAR))
   {
      if(bKind(fId) == SparseSM::Bessel::BesselKind::JlK0)
      {
         nBc = 0;
      }
      else
      {
         nBc = 1;
      }
   }
   else if (fId == std::make_pair(PhysicalNames::Velocity::id(),
                      FieldComponents::Spectral::POL))
   {
      if(bKind(fId) == SparseSM::Bessel::BesselKind::JlK2)
      {
         nBc = 2;
      }
      else
      {
         nBc = 1;
      }
   }
   else
   {
      nBc = 0;
   }

   return nBc;
}

void IRTCBackend::applyTau(SparseMatrix& mat, const SpectralFieldId& rowId,
   const SpectralFieldId& colId, const int l, const Resolution& res,
   const BcMap& bcs, const NonDimensional::NdMap& nds,
   const bool isSplitOperator) const
{
   auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);

   auto bcId = bcs.find(rowId.first)->second;

   if (rowId == std::make_pair(PhysicalNames::Velocity::id(),
                   FieldComponents::Spectral::TOR) &&
       rowId == colId)
   {
      if (l > 0)
      {
         if(this->nBc(colId) == 1)
         {
            bool atTop = false;
            SparseSM::Bessel::Boundary::Operator bcOp(nN, nN, bKind(colId), l, atTop);
            if (bcId == Bc::Name::NoSlip::id())
            {
               bcOp.addRow<SparseSM::Bessel::Boundary::Value>();
            }
            else if (bcId == Bc::Name::StressFree::id())
            {
               throw std::logic_error("Stress-free tau line is not implemented yet");
            }
            else
            {
               throw std::logic_error(
                  "Boundary conditions for Velocity Toroidal component "
                  "not implemented");
            }
            mat.real() += bcOp.mat();
         }
      }
   }
   else if (rowId == std::make_pair(PhysicalNames::Velocity::id(),
                        FieldComponents::Spectral::POL) &&
            rowId == colId)
   {
      if (l > 0)
      {
         bool atTop = false;
         SparseSM::Bessel::Boundary::Operator bcOp(nN, nN, bKind(colId), l, atTop);
         if (this->useSplitEquation())
         {
            if (isSplitOperator)
            {
               if(this->nBc(colId) == 2)
               {
                  bcOp.addRow<SparseSM::Bessel::Boundary::Value>();
               }

            }
            else if (bcId == Bc::Name::NoSlip::id())
            {
               bcOp.addRow<SparseSM::Bessel::Boundary::Value>();
            }
            else if (bcId == Bc::Name::StressFree::id())
            {
               bcOp.addRow<SparseSM::Bessel::Boundary::D2>();
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
            if(this->nBc(colId) == 2)
            {
               bcOp.addRow<SparseSM::Bessel::Boundary::Value>();
            }

            if (bcId == Bc::Name::NoSlip::id())
            {
               bcOp.addRow<SparseSM::Bessel::Boundary::D1>();
            }
            else if (bcId == Bc::Name::StressFree::id())
            {
               bcOp.addRow<SparseSM::Bessel::Boundary::D2>();
            }
            else
            {
               throw std::logic_error(
                  "Boundary conditions for Velocity Poloidal component "
                  "not implemented");
            }
         }
         mat.real() += bcOp.mat();
      }
   }
   else if (rowId == std::make_pair(PhysicalNames::Temperature::id(),
                        FieldComponents::Spectral::SCALAR) &&
            rowId == colId)
   {
      if(this->nBc(colId) == 1)
      {
         bool atTop = false;
         SparseSM::Bessel::Boundary::Operator bcOp(nN, nN, bKind(colId), l, atTop);
         if (bcId == Bc::Name::FixedTemperature::id())
         {
            bcOp.addRow<SparseSM::Bessel::Boundary::Value>();
         }
         else if (bcId == Bc::Name::FixedFlux::id())
         {
            bcOp.addRow<SparseSM::Bessel::Boundary::D1>();
         }
         else
         {
            throw std::logic_error(
                  "Boundary conditions for temperature "
                  "not implemented");
         }
         mat.real() += bcOp.mat();
      }
   }
}

void IRTCBackend::stencil(SparseMatrix& mat, const SpectralFieldId& fieldId,
   const int l, const Resolution& res, const bool makeSquare, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);

   auto bcId = bcs.find(fieldId.first)->second;

   int s = this->nBc(fieldId);
   if (fieldId == std::make_pair(PhysicalNames::Velocity::id(),
                     FieldComponents::Spectral::TOR))
   {
      SparseSM::Id qid(nN, nN - s);
      mat = qid.mat();
   }
   else if (fieldId == std::make_pair(PhysicalNames::Velocity::id(),
                          FieldComponents::Spectral::POL))
   {
      if (bcId == Bc::Name::NoSlip::id())
      {
         //SparseSM::Worland::Stencil::ValueD1 bc(nN, nN - s, a, b, l);
         //mat = bc.mat();
      }
      else if (bcId == Bc::Name::StressFree::id())
      {
         //SparseSM::Worland::Stencil::ValueD2 bc(nN, nN - s, a, b, l);
         //mat = bc.mat();
      }
      else
      {
         throw std::logic_error("Galerin boundary conditions for Velocity "
                                "Poloidal component not implemented");
      }
   }
   else if (fieldId == std::make_pair(PhysicalNames::Temperature::id(),
                          FieldComponents::Spectral::SCALAR))
   {
      SparseSM::Id qid(nN, nN - s);
      mat = qid.mat();
   }
}

void IRTCBackend::applyGalerkinStencil(SparseMatrix& mat,
   const SpectralFieldId& rowId, const SpectralFieldId& colId, const int lr,
   const int lc, const Resolution& res, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   auto nNr = res.counter().dimensions(Dimensions::Space::SPECTRAL, lr)(0);

   auto S = mat;
   this->stencil(S, colId, lc, res, false, bcs, nds);

   auto s = this->nBc(rowId);
   SparseSM::Id qId(nNr - s, nNr, 0, s);
   mat = qId.mat() * (mat * S);
}

} // namespace RTC
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
