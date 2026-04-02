/**
 * @file IRTCBackend.cpp
 * @brief Source of the interface for model backend
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "Model/Boussinesq/Sphere/RTC/IRTCBackend.hpp"
#include "QuICC/Bc/Name/FixedFlux.hpp"
#include "QuICC/Bc/Name/FixedTemperature.hpp"
#include "QuICC/Bc/Name/NoSlip.hpp"
#include "QuICC/Bc/Name/StressFree.hpp"
#include "QuICC/Bc/Name/QuasiInverseOnly.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/NonDimensional/CflInertial.hpp"
#include "QuICC/NonDimensional/Ekman.hpp"
#include "QuICC/NonDimensional/Prandtl.hpp"
#include "QuICC/NonDimensional/Rayleigh.hpp"
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
#include "QuICC/Tools/IdToHuman.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace RTC {

std::vector<std::string> IRTCBackend::fieldNames() const
{
   std::vector<std::string> names = {
      PhysicalNames::Velocity().tag(),
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

   auto vel_tor = std::make_pair(PhysicalNames::Velocity::id(),
                   FieldComponents::Spectral::TOR);
   auto vel_pol = std::make_pair(PhysicalNames::Velocity::id(),
                   FieldComponents::Spectral::POL);
   auto temp = std::make_pair(PhysicalNames::Temperature::id(),
                   FieldComponents::Spectral::SCALAR);

   if (fId == vel_tor ||
       fId == temp)
   {
      nBc = 1;
   }
   else if (fId == vel_pol)
   {
      nBc = 2;
   }
   else
   {
      nBc = 0;
   }

   return nBc;
}

int IRTCBackend::baseNn(const int l, const Resolution& res) const
{
   int nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);

   return nN;
}

void IRTCBackend::applyTau(SparseMatrix& mat, const SpectralFieldId& rowId,
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
   auto vel_pol = std::make_pair(PhysicalNames::Velocity::id(),
                   FieldComponents::Spectral::POL);
   auto temp = std::make_pair(PhysicalNames::Temperature::id(),
                   FieldComponents::Spectral::SCALAR);

   if ((rowId == vel_tor) &&
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
   else if ((rowId == vel_pol) &&
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
   else if ((rowId == temp) &&
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

void IRTCBackend::stencil(SparseMatrix& mat, const SpectralFieldId& fieldId,
   const int l, const int nN, const bool makeSquare, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   auto a = Polynomial::Worland::worland_default_t::ALPHA;
   auto b = Polynomial::Worland::worland_default_t::DBETA;

   auto bcId = bcs.find(fieldId.first)->second;

   auto vel_tor = std::make_pair(PhysicalNames::Velocity::id(),
                   FieldComponents::Spectral::TOR);
   auto vel_pol = std::make_pair(PhysicalNames::Velocity::id(),
                   FieldComponents::Spectral::POL);
   auto temp = std::make_pair(PhysicalNames::Temperature::id(),
                   FieldComponents::Spectral::SCALAR);

   int s = this->nBc(fieldId);
   if(bcId == Bc::Name::QuasiInverseOnly::id())
   {
      SparseSM::Worland::Id qid(nN, nN - s, a, b, l);
      mat = qid.mat();
   }
   else
   {
      if (fieldId == vel_tor)
      {
         if (bcId == Bc::Name::NoSlip::id())
         {
            SparseSM::Worland::Stencil::Value bc(nN, nN - s, a, b, l);
            mat = bc.mat();
         }
         else if (bcId == Bc::Name::StressFree::id())
         {
            SparseSM::Worland::Stencil::R1D1DivR1 bc(nN, nN - s, a, b, l);
            mat = bc.mat();
         }
         else
         {
            throw std::logic_error("Galerkin boundary conditions for Velocity "
                  "Toroidal component not implemented");
         }
      }
      else if (fieldId == vel_pol)
      {
         if (bcId == Bc::Name::NoSlip::id())
         {
            SparseSM::Worland::Stencil::ValueD1 bc(nN, nN - s, a, b, l);
            mat = bc.mat();
         }
         else if (bcId == Bc::Name::StressFree::id())
         {
            SparseSM::Worland::Stencil::ValueD2 bc(nN, nN - s, a, b, l);
            mat = bc.mat();
         }
         else
         {
            throw std::logic_error("Galerin boundary conditions for Velocity "
                  "Poloidal component not implemented");
         }
      }
      else if (fieldId == temp)
      {
         if (bcId == Bc::Name::FixedTemperature::id())
         {
            SparseSM::Worland::Stencil::Value bc(nN, nN - s, a, b, l);
            mat = bc.mat();
         }
         else if (bcId == Bc::Name::FixedFlux::id())
         {
            SparseSM::Worland::Stencil::D1 bc(nN, nN - s, a, b, l);
            mat = bc.mat();
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

void IRTCBackend::applyGalerkinStencil(SparseMatrix& mat,
   const SpectralFieldId& rowId, const SpectralFieldId& colId, const int lr,
   const int lc, std::shared_ptr<details::BlockOptions> opts,
   const int nNr, const int nNc, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
   auto a = Polynomial::Worland::worland_default_t::ALPHA;
   auto b = Polynomial::Worland::worland_default_t::DBETA;

   auto S = mat;
   this->stencil(S, colId, lc, nNc, false, bcs, nds);

   auto s = this->nBc(rowId);
   SparseSM::Worland::Id qId(nNr - s, nNr, a, b, lr, 0, s);
   mat = qId.mat() * (mat * S);
}

void IRTCBackend::operatorInfo(OperatorInfo& info, const SpectralFieldId& fId,
   const Resolution& res, const Equations::Tools::ICoupling& coupling,
   const BcMap& bcs) const
{
   // Loop overall matrices/eigs
   for (int idx = 0; idx < info.tauN.size(); ++idx)
   {
      auto eigs = coupling.getIndexes(res, idx);

      int tN, gN, rhs;
      ArrayI shift(3);

      auto nTauLines = this->nBc(fId);
      auto nN = this->baseNn(eigs.at(0), res);
      this->blockInfo(tN, gN, shift, rhs, nTauLines, nN, this->useGalerkin());

      info.tauN(idx) = tN;
      info.galN(idx) = gN;
      info.galShift.row(idx) = shift;
      info.rhsCols(idx) = rhs;

      // Compute system size
      int sN = 0;
      for (auto f: this->implicitFields(fId))
      {
         nTauLines = this->nBc(f);
         this->blockInfo(tN, gN, shift, rhs, nTauLines, nN, this->useGalerkin());
         sN += gN;
      }

      if (sN == 0)
      {
         sN = info.galN(idx);
      }

      info.sysN(idx) = sN;
   }
}

} // namespace RTC
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
