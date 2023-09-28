/**
 * @file ModelBackend.cpp
 * @brief Source of the interface for model backend
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/Enums/Dimensions.hpp"
#include "QuICC/Model/Boussinesq/Sphere/RTC/Implicit/ModelBackend.hpp"
#include "QuICC/ModelOperator/Time.hpp"
#include "QuICC/ModelOperator/ImplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"
#include "QuICC/ModelOperator/Stencil.hpp"
#include "QuICC/ModelOperator/Boundary.hpp"
#include "QuICC/ModelOperator/SplitImplicitLinear.hpp"
#include "QuICC/ModelOperator/SplitBoundary.hpp"
#include "QuICC/ModelOperator/SplitBoundaryValue.hpp"
#include "QuICC/ModelOperatorBoundary/FieldToRhs.hpp"
#include "QuICC/ModelOperatorBoundary/SolverHasBc.hpp"
#include "QuICC/ModelOperatorBoundary/SolverNoTau.hpp"
#include "QuICC/ModelOperatorBoundary/Stencil.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Bc/Name/FixedTemperature.hpp"
#include "QuICC/Bc/Name/FixedFlux.hpp"
#include "QuICC/Bc/Name/StressFree.hpp"
#include "QuICC/Bc/Name/NoSlip.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/NonDimensional/Prandtl.hpp"
#include "QuICC/NonDimensional/Rayleigh.hpp"
#include "QuICC/NonDimensional/Ekman.hpp"
#include "QuICC/NonDimensional/CflInertial.hpp"
#include "QuICC/Tools/IdToHuman.hpp"
#include "QuICC/Resolutions/Tools/IndexCounter.hpp"
#include "QuICC/SparseSM/Worland/Id.hpp"
#include "QuICC/SparseSM/Worland/I2.hpp"
#include "QuICC/SparseSM/Worland/I4.hpp"
#include "QuICC/SparseSM/Worland/I2Lapl.hpp"
#include "QuICC/SparseSM/Worland/I2Qp.hpp"
#include "QuICC/SparseSM/Worland/I2Qm.hpp"
#include "QuICC/SparseSM/Worland/I4Lapl.hpp"
#include "QuICC/SparseSM/Worland/I4Lapl2.hpp"
#include "QuICC/SparseSM/Worland/I4Qm.hpp"
#include "QuICC/SparseSM/Worland/I4Qp.hpp"
#include "QuICC/SparseSM/Worland/Boundary/Value.hpp"
#include "QuICC/SparseSM/Worland/Boundary/D1.hpp"
#include "QuICC/SparseSM/Worland/Boundary/D2.hpp"
#include "QuICC/SparseSM/Worland/Boundary/R1D1DivR1.hpp"
#include "QuICC/SparseSM/Worland/Boundary/Operator.hpp"
#include "QuICC/SparseSM/Worland/Stencil/Value.hpp"
#include "QuICC/SparseSM/Worland/Stencil/D1.hpp"
#include "QuICC/SparseSM/Worland/Stencil/R1D1DivR1.hpp"
#include "QuICC/SparseSM/Worland/Stencil/ValueD1.hpp"
#include "QuICC/SparseSM/Worland/Stencil/ValueD2.hpp"

#include "QuICC/Polynomial/Worland/WorlandBase.hpp"
#include "QuICC/Equations/CouplingIndexType.hpp"

//#define QUICC_INVISCID_4TH
#define QUICC_INVISCID_2ND
#include <iostream>

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace RTC {

namespace Implicit {

   ModelBackend::ModelBackend()
      : IRTCBackend(),
#ifdef QUICC_TRANSFORM_WORLAND_TRUNCATE_QI
      mcTruncateQI(true)
#else
      mcTruncateQI(false)
#endif // QUICC_TRANSFORM_WORLAND_TRUNCATE_QI
   {
   }

   void ModelBackend::enableSplitEquation(const bool flag)
   {
      if(flag)
      {
         throw std::logic_error("Split equation for implicit model is not implemented");
      }
      else
      {
         IRTCBackend::enableSplitEquation(flag);
      }
   }

   ModelBackend::SpectralFieldIds ModelBackend::implicitFields(const SpectralFieldId& fId) const
   {
      SpectralFieldId velTor = std::make_pair(PhysicalNames::Velocity::id(), FieldComponents::Spectral::TOR);
      SpectralFieldId velPol = std::make_pair(PhysicalNames::Velocity::id(), FieldComponents::Spectral::POL);
      SpectralFieldId temp = std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR);
      SpectralFieldIds fields = {velTor, velPol, temp};

      // sort the fields
      std::sort(fields.begin(), fields.end());
      return fields;
   }

   void ModelBackend::equationInfo(EquationInfo& info, const SpectralFieldId& fId, const Resolution& res) const
   {
      // Operators are real
      info.isComplex = true;

      // Splitting 4th poloidal equation into two systems
      if(fId == std::make_pair(PhysicalNames::Velocity::id(), FieldComponents::Spectral::POL))
      {
         info.isSplitEquation = this->useSplitEquation();
      }
      else
      {
         info.isSplitEquation = false;
      }

      // Implicit coupled fields
      info.im = this->implicitFields(fId);

      // Explicit linear terms
      info.exL.clear();

      // Explicit nonlinear terms
      info.exNL.clear();

      // Explicit nextstep terms
      info.exNS.clear();

      // Index mode
      info.indexMode = static_cast<int>(Equations::CouplingIndexType::SLOWEST_SINGLE_RHS);
   }

   void ModelBackend::operatorInfo(OperatorInfo& info, const SpectralFieldId& fId, const Resolution& res, const Equations::Tools::ICoupling& coupling, const BcMap& bcs) const
   {
      // Loop overall matrices/eigs
      for(int idx = 0; idx < info.tauN.size(); ++idx)
      {
         auto eigs = coupling.getIndexes(res, idx);

         int tN, gN, rhs;
         ArrayI shift(3);

         this->blockInfo(tN, gN, shift, rhs, fId, res, eigs.at(0), bcs);

         info.tauN(idx) = tN;
         info.galN(idx) = gN;
         info.galShift.row(idx) = shift;
         info.rhsCols(idx) = rhs;

         // Compute system size
         int sN = 0;
         for(auto f: this->implicitFields(fId))
         {
            this->blockInfo(tN, gN, shift, rhs, f, res, eigs.at(0), bcs);
            sN += gN;
         }

         if(sN == 0)
         {
            sN = info.galN(idx);
         }

         info.sysN(idx) = sN;
      }
   }

   void ModelBackend::buildBlock(DecoupledZSparse& decMat, const std::vector<internal::BlockDescription>& descr, const SpectralFieldId& rowId, const SpectralFieldId& colId, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds, const bool isSplitOperator) const
   {
      assert(eigs.size() == 1);
      int m = eigs.at(0);

      // Compute system size
      const auto sysInfo = systemInfo(rowId, colId, m, res, bcs, this->useGalerkin(), false);
      const auto& sysN = sysInfo.systemSize;
      const auto& baseRowShift = sysInfo.startRow;
      const auto& baseColShift = sysInfo.startCol;
      auto nL = res.counter().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL, m);

      // Resize matrices the first time
      if(decMat.real().size() == 0)
      {
         decMat.real().resize(sysN,sysN);
         decMat.imag().resize(sysN,sysN);
      }
      assert(decMat.real().rows() == sysN);
      assert(decMat.real().cols() == sysN);
      assert(decMat.imag().rows() == sysN);
      assert(decMat.imag().cols() == sysN);

      int tN, gN, rhs;
      ArrayI shift(3);

      bool needTau = (bcType == ModelOperatorBoundary::SolverHasBc::id());

      for(auto&& d: descr)
      {
         assert(d.nRowShift == 0 ||  d.nColShift == 0);

         // Shift starting row
         int rowShift = baseRowShift;
         for(int s = 0; s < d.nRowShift; s++)
         {
            this->blockInfo(tN, gN, shift, rhs, rowId, res, m+s, bcs);
            rowShift += gN;
         }

         // Shift starting col
         int colShift = baseColShift;
         for(int s = 0; s < d.nColShift; s++)
         {
            this->blockInfo(tN, gN, shift, rhs, colId, res, m+s, bcs);
            colShift += gN;
         }

         int lShift = -d.nRowShift + d.nColShift;

         for(int l = m + d.nRowShift; l < nL - d.nColShift; l++)
         {
            auto nNr = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
            auto nNc = res.counter().dimensions(Dimensions::Space::SPECTRAL, l + lShift)(0);

            //
            // Build real part of block
            if(d.realOp)
            {
               auto bMat = d.realOp(nNr, nNc, l, d.opts, nds);

               if(this->useGalerkin())
               {
                  this->applyGalerkinStencil(bMat, rowId, colId, l, l + lShift, res, bcs, nds);
               }
               else if(needTau)
               {
                  this->applyTau(bMat, rowId, colId, l + lShift, res, bcs, nds, isSplitOperator);
               }
               this->addBlock(decMat.real(), bMat, rowShift, colShift);
            }

            //
            // Build imaginary part of block
            if(d.imagOp)
            {
               auto bMat = d.imagOp(nNr, nNc, l, d.opts, nds);

               if(this->useGalerkin())
               {
                  this->applyGalerkinStencil(bMat, rowId, colId, l, l + lShift, res, bcs, nds);
               }
               else if(needTau)
               {
                  this->applyTau(bMat, rowId, colId, l + lShift, res, bcs, nds, isSplitOperator);
               }
               this->addBlock(decMat.imag(), bMat, rowShift, colShift);
            }

            // Shift to next block
            this->blockInfo(tN, gN, shift, rhs, rowId, res, l, bcs);
            rowShift += gN;

            this->blockInfo(tN, gN, shift, rhs, colId, res, l + lShift, bcs);
            colShift += gN;
         }
      }
   }

   std::vector<internal::BlockDescription> ModelBackend::implicitBlockBuilder(const SpectralFieldId& rowId, const SpectralFieldId& colId, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds, const bool isSplitOperator) const
   {
      std::vector<internal::BlockDescription> descr;

      // Create description with common options
      auto getDescription = [&]()->internal::BlockDescription&
      {
            descr.push_back({});
            auto& d = descr.back();
            auto opts = std::make_shared<internal::BlockOptionsImpl>();
            opts->a = Polynomial::Worland::WorlandBase::ALPHA_CHEBYSHEV;
            opts->b = Polynomial::Worland::WorlandBase::DBETA_CHEBYSHEV;
            opts->m = eigs.at(0);
            opts->bcId = bcs.find(colId.first)->second;
            opts->truncateQI = this->mcTruncateQI;
            opts->isSplitOperator = isSplitOperator;
            d.opts = opts;

            return d;
      };

      if(rowId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::TOR))
      {
         if(rowId == colId)
         {
            // Real part of operator
            auto realOp = [](const int nNr, const int nNc, const int l, std::shared_ptr<internal::BlockOptions> opts, const NonDimensional::NdMap& nds)
            {
               SparseMatrix bMat(nNr, nNc);

               if(l > 0)
               {
                  auto& o = *std::dynamic_pointer_cast<internal::BlockOptionsImpl>(opts);
                  SparseSM::Worland::I2Lapl i2lapl(nNr, nNc, o.a, o.b, l, 1*o.truncateQI);
                  bMat = i2lapl.mat();
               }

               return bMat;
            };

            // Imaginary part of operator
            auto imagOp = [](const int nNr, const int nNc, const int l, std::shared_ptr<internal::BlockOptions> opts, const NonDimensional::NdMap& nds)
            {
               SparseMatrix bMat(nNr, nNc);

               if(l > 0)
               {
                  auto& o = *std::dynamic_pointer_cast<internal::BlockOptionsImpl>(opts);

                  const auto T = 1.0/nds.find(NonDimensional::Ekman::id())->second->value();
                  const auto dl = static_cast<MHDFloat>(l);
                  const auto invlapl = 1.0/(dl*(dl + 1.0));

                  SparseSM::Worland::I2 i2(nNr, nNc, o.a, o.b, l, 1*o.truncateQI);
                  bMat = o.m*T*invlapl*i2.mat();
               }

               return bMat;
            };

            // Imaginary part of operator
            auto imagOp2 = [](const int nNr, const int nNc, const int l, std::shared_ptr<internal::BlockOptions> opts, const NonDimensional::NdMap& nds)
            {
               SparseMatrix bMat(nNr, nNc);

               if(l > 0)
               {
                  auto& o = *std::dynamic_pointer_cast<internal::BlockOptionsImpl>(opts);

                  const auto T = 1.0/nds.find(NonDimensional::Ekman::id())->second->value();
                  const auto dl = static_cast<MHDFloat>(l);
                  const auto invlapl = 1.0/(dl*(dl + 1.0));

                  SparseSM::Worland::I2 i2(nNr+1, nNc, o.a, o.b, l, 1*o.truncateQI);
                  SparseSM::Worland::Id qid(nNr, nNr+1, o.a, o.b, l, 0, 1);
                  bMat = o.m*T*invlapl*(qid.mat()*i2.mat());
               }

               return bMat;
            };

            // Create block diagonal operator
            auto& d = getDescription();
            d.nRowShift = 0;
            d.nColShift = 0;
#if defined QUICC_INVISCID_4TH
            d.realOp = nullptr;
            d.imagOp = imagOp;
#elif defined QUICC_INVISCID_2ND
            d.realOp = nullptr;
            d.imagOp = imagOp2;
#else
            d.realOp = realOp;
            d.imagOp = imagOp;
#endif
         }
         else if(colId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL))
         {
            // Real part of first lower diagonal
            auto realOpLower = [](const int nNr, const int nNc, const int l, std::shared_ptr<internal::BlockOptions> opts, const NonDimensional::NdMap& nds)
            {
               SparseMatrix bMat(nNr, nNc);

               if(l > 0)
               {
                  auto& o = *std::dynamic_pointer_cast<internal::BlockOptionsImpl>(opts);

                  auto coriolis = [](const int l, const int m){
                     return (l - MHD_MP(1.0))*(l + MHD_MP(1.0))*precision::sqrt(((l - m)*(l + m))/((MHD_MP(2.0)*l - MHD_MP(1.0))*(MHD_MP(2.0)*l + MHD_MP(1.0))));
                  };

                  const auto T = 1.0/nds.find(NonDimensional::Ekman::id())->second->value();
                  const auto dl = static_cast<QuICC::internal::MHDFloat>(l);
                  const auto invlapl = 1.0/(dl*(dl + 1.0));

                  SparseSM::Worland::I2Qm corQm(nNr, nNc, o.a, o.b, l, 1*o.truncateQI);
                  auto norm = coriolis(l, o.m);
                  bMat = -static_cast<MHDFloat>(norm*T*invlapl)*corQm.mat();
               }

               return bMat;
            };

            // Real part of first lower diagonal
            auto realOpLower2 = [](const int nNr, const int nNc, const int l, std::shared_ptr<internal::BlockOptions> opts, const NonDimensional::NdMap& nds)
            {
               SparseMatrix bMat(nNr, nNc);

               if(l > 0)
               {
                  auto& o = *std::dynamic_pointer_cast<internal::BlockOptionsImpl>(opts);

                  auto coriolis = [](const int l, const int m){
                     return (l - MHD_MP(1.0))*(l + MHD_MP(1.0))*precision::sqrt(((l - m)*(l + m))/((MHD_MP(2.0)*l - MHD_MP(1.0))*(MHD_MP(2.0)*l + MHD_MP(1.0))));
                  };

                  const auto T = 1.0/nds.find(NonDimensional::Ekman::id())->second->value();
                  const auto dl = static_cast<QuICC::internal::MHDFloat>(l);
                  const auto invlapl = 1.0/(dl*(dl + 1.0));

                  SparseSM::Worland::I2Qm corQm(nNr+1, nNc, o.a, o.b, l, 1*o.truncateQI);
                  int s = 0;
                  if(l%2 == 0)
                  {
                     s = 1;
                  }

                  SparseSM::Worland::Id qid(nNr, nNr+1, o.a, o.b, l, 0, s);
                  auto norm = coriolis(l, o.m);
                  bMat = -static_cast<MHDFloat>(norm*T*invlapl)*(qid.mat()*corQm.mat());
               }

               return bMat;
            };

            // Create first lower diagonal operator
            auto& dLow = getDescription();
            dLow.nRowShift = 1;
            dLow.nColShift = 0;
#ifdef QUICC_INVISCID_2ND
            dLow.realOp = realOpLower2;
            dLow.imagOp = nullptr;
#else
            dLow.realOp = realOpLower;
            dLow.imagOp = nullptr;
#endif

            // Real part of first upper diagonal
            auto realOpUpper = [](const int nNr, const int nNc, const int l, std::shared_ptr<internal::BlockOptions> opts, const NonDimensional::NdMap& nds)
            {
               SparseMatrix bMat(nNr, nNc);

               if(l > 0)
               {
                  auto& o = *std::dynamic_pointer_cast<internal::BlockOptionsImpl>(opts);

                  auto coriolis = [](const int l, const int m){
                     return (l - MHD_MP(1.0))*(l + MHD_MP(1.0))*precision::sqrt(((l - m)*(l + m))/((MHD_MP(2.0)*l - MHD_MP(1.0))*(MHD_MP(2.0)*l + MHD_MP(1.0))));
                  };

                  const auto T = 1.0/nds.find(NonDimensional::Ekman::id())->second->value();
                  const auto dl = static_cast<QuICC::internal::MHDFloat>(l);
                  const auto invlapl = MHD_MP(1.0)/(dl*(dl + MHD_MP(1.0)));
                  SparseSM::Worland::I2Qp corQp(nNr, nNc, o.a, o.b, l, 1*o.truncateQI);
                  auto norm = -coriolis(l+1, o.m);
                  bMat = -static_cast<MHDFloat>(norm*T*invlapl)*corQp.mat();
               }

               return bMat;
            };

            // Real part of first upper diagonal
            auto realOpUpper2 = [](const int nNr, const int nNc, const int l, std::shared_ptr<internal::BlockOptions> opts, const NonDimensional::NdMap& nds)
            {
               SparseMatrix bMat(nNr, nNc);

               if(l > 0)
               {
                  auto& o = *std::dynamic_pointer_cast<internal::BlockOptionsImpl>(opts);

                  auto coriolis = [](const int l, const int m){
                     return (l - MHD_MP(1.0))*(l + MHD_MP(1.0))*precision::sqrt(((l - m)*(l + m))/((MHD_MP(2.0)*l - MHD_MP(1.0))*(MHD_MP(2.0)*l + MHD_MP(1.0))));
                  };

                  const auto T = 1.0/nds.find(NonDimensional::Ekman::id())->second->value();
                  const auto dl = static_cast<QuICC::internal::MHDFloat>(l);
                  const auto invlapl = MHD_MP(1.0)/(dl*(dl + MHD_MP(1.0)));
                  SparseSM::Worland::I2Qp corQp(nNr+1, nNc, o.a, o.b, l, 1*o.truncateQI);
                  int s = 0;
                  if(l%2 == 0)
                  {
                     s = 1;
                  }
                  SparseSM::Worland::Id qid(nNr, nNr+1, o.a, o.b, l, 0, s);
                  auto norm = -coriolis(l+1, o.m);
                  bMat = -static_cast<MHDFloat>(norm*T*invlapl)*(qid.mat()*corQp.mat());
               }

               return bMat;
            };

            // Create first upper diagonal operator
            auto& dUp = getDescription();
            dUp.nRowShift = 0;
            dUp.nColShift = 1;
#ifdef QUICC_INVISCID_2ND
            dUp.realOp = realOpUpper2;
            dUp.imagOp = nullptr;
#else
            dUp.realOp = realOpUpper;
            dUp.imagOp = nullptr;
#endif
         }
      }
      else if(rowId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL))
      {
         if(rowId == colId)
         {
            // Real part of block
            auto realOp = [](const int nNr, const int nNc, const int l, std::shared_ptr<internal::BlockOptions> opts, const NonDimensional::NdMap& nds)
            {
               SparseMatrix bMat(nNr, nNc);

               if(l > 0)
               {
                  auto& o = *std::dynamic_pointer_cast<internal::BlockOptionsImpl>(opts);

                  SparseSM::Worland::I4Lapl2 i4lapl2(nNr, nNc, o.a, o.b, l, 2*o.truncateQI);
                  bMat = i4lapl2.mat();
               }

               return bMat;
            };

            // Imaginary part of block
            auto imagOp = [](const int nNr, const int nNc, const int l, std::shared_ptr<internal::BlockOptions> opts, const NonDimensional::NdMap& nds)
            {
               SparseMatrix bMat(nNr, nNc);

               if(l > 0)
               {
                  auto& o = *std::dynamic_pointer_cast<internal::BlockOptionsImpl>(opts);

                  const auto dl = static_cast<QuICC::internal::MHDFloat>(l);
                  const auto invlapl = MHD_MP(1.0)/(dl*(dl + MHD_MP(1.0)));
                  const auto T = 1.0/nds.find(NonDimensional::Ekman::id())->second->value();
                  SparseSM::Worland::I4Lapl coriolis(nNr, nNc, o.a, o.b, l, 2*o.truncateQI);

                  // Correct Laplacian for 4th order system according to:
                  // McFadden,Murray,Boisvert,
                  // Elimination of Spurious Eigenvalues in the
                  // Chebyshev Tau Spectral Method,
                  // JCP 91, 228-239 (1990)
                  // We simply drop the last column
                  if(o.bcId == Bc::Name::NoSlip::id())
                  {
                     SparseSM::Worland::Id qid(nNr, nNc, o.a, o.b, l, -1);
                     bMat = static_cast<MHDFloat>(o.m*T*invlapl)*coriolis.mat()*qid.mat();
                  }
                  else
                  {
                     bMat = static_cast<MHDFloat>(o.m*T*invlapl)*coriolis.mat();
                  }
               }

               return bMat;
            };

            // Imaginary part of block
            auto imagOp2 = [](const int nNr, const int nNc, const int l, std::shared_ptr<internal::BlockOptions> opts, const NonDimensional::NdMap& nds)
            {
               SparseMatrix bMat(nNr, nNc);

               if(l > 0)
               {
                  auto& o = *std::dynamic_pointer_cast<internal::BlockOptionsImpl>(opts);

                  const auto dl = static_cast<QuICC::internal::MHDFloat>(l);
                  const auto invlapl = MHD_MP(1.0)/(dl*(dl + MHD_MP(1.0)));
                  const auto T = 1.0/nds.find(NonDimensional::Ekman::id())->second->value();
                  SparseSM::Worland::I2Lapl coriolis(nNr, nNc, o.a, o.b, l, 2*o.truncateQI);

                  bMat = static_cast<MHDFloat>(o.m*T*invlapl)*coriolis.mat();
               }

               return bMat;
            };

            // Create diagonal block
            auto& d = getDescription();
            d.nRowShift = 0;
            d.nColShift = 0;
#if defined QUICC_INVISCID_4TH
            d.realOp = nullptr;
            d.imagOp = imagOp;
#elif defined QUICC_INVISCID_2ND
            d.realOp = nullptr;
            d.imagOp = imagOp2;
#else
            d.realOp = realOp;
            d.imagOp = imagOp;
#endif

         }
         else if(colId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::TOR))
         {
            // Create real part of block
            auto realOpLower = [](const int nNr, const int nNc, const int l, std::shared_ptr<internal::BlockOptions> opts, const NonDimensional::NdMap& nds)
            {
               SparseMatrix bMat(nNr, nNc);

               if(l > 0)
               {
                  auto& o = *std::dynamic_pointer_cast<internal::BlockOptionsImpl>(opts);

                  auto coriolis = [](const int l, const int m){
                     return (l - MHD_MP(1.0))*(l + MHD_MP(1.0))*precision::sqrt(((l - m)*(l + m))/((MHD_MP(2.0)*l - MHD_MP(1.0))*(MHD_MP(2.0)*l + MHD_MP(1.0))));
                  };

                  const auto T = 1.0/nds.find(NonDimensional::Ekman::id())->second->value();

                  const auto dl = static_cast<MHDFloat>(l);
                  const auto invlapl = 1.0/(dl*(dl + 1.0));
                  SparseSM::Worland::I4Qm corQm(nNr, nNc, o.a, o.b, l, 2*o.truncateQI);
                  auto norm = coriolis(l, o.m);
                  bMat = static_cast<MHDFloat>(norm*T*invlapl)*corQm.mat();
               }

               return bMat;
            };

            // Create real part of block
            auto realOpLower2 = [](const int nNr, const int nNc, const int l, std::shared_ptr<internal::BlockOptions> opts, const NonDimensional::NdMap& nds)
            {
               SparseMatrix bMat(nNr, nNc);

               if(l > 0)
               {
                  auto& o = *std::dynamic_pointer_cast<internal::BlockOptionsImpl>(opts);

                  auto coriolis = [](const int l, const int m){
                     return (l - MHD_MP(1.0))*(l + MHD_MP(1.0))*precision::sqrt(((l - m)*(l + m))/((MHD_MP(2.0)*l - MHD_MP(1.0))*(MHD_MP(2.0)*l + MHD_MP(1.0))));
                  };

                  const auto T = 1.0/nds.find(NonDimensional::Ekman::id())->second->value();

                  const auto dl = static_cast<MHDFloat>(l);
                  const auto invlapl = 1.0/(dl*(dl + 1.0));
                  SparseSM::Worland::I2Qm corQm(nNr, nNc, o.a, o.b, l, 1*o.truncateQI);
                  auto norm = coriolis(l, o.m);
                  bMat = static_cast<MHDFloat>(norm*T*invlapl)*corQm.mat();
               }

               return bMat;
            };

            // Create first lower diagonal operator
            auto& dLow = getDescription();
            dLow.nRowShift = 1;
            dLow.nColShift = 0;
#ifdef QUICC_INVISCID_2ND
            dLow.realOp = realOpLower2;
            dLow.imagOp = nullptr;
#else
            dLow.realOp = realOpLower;
            dLow.imagOp = nullptr;
#endif

            // Create real part of block
            auto realOpUpper = [](const int nNr, const int nNc, const int l, std::shared_ptr<internal::BlockOptions> opts, const NonDimensional::NdMap& nds)
            {
               SparseMatrix bMat(nNr, nNc);

               if(l > 0)
               {
                  auto& o = *std::dynamic_pointer_cast<internal::BlockOptionsImpl>(opts);

                  auto coriolis = [](const int l, const int m){
                     return (l - MHD_MP(1.0))*(l + MHD_MP(1.0))*precision::sqrt(((l - m)*(l + m))/((MHD_MP(2.0)*l - MHD_MP(1.0))*(MHD_MP(2.0)*l + MHD_MP(1.0))));
                  };

                  const auto T = 1.0/nds.find(NonDimensional::Ekman::id())->second->value();

                  const auto dl = static_cast<MHDFloat>(l);
                  const auto invlapl = 1.0/(dl*(dl + 1.0));
                  SparseSM::Worland::I4Qp corQp(nNr, nNc, o.a, o.b, l, 2*o.truncateQI);
                  auto norm = -coriolis(l+1, o.m);
                  bMat = static_cast<MHDFloat>(norm*T*invlapl)*corQp.mat();
               }

               return bMat;
            };

            // Create real part of block
            auto realOpUpper2 = [](const int nNr, const int nNc, const int l, std::shared_ptr<internal::BlockOptions> opts, const NonDimensional::NdMap& nds)
            {
               SparseMatrix bMat(nNr, nNc);

               if(l > 0)
               {
                  auto& o = *std::dynamic_pointer_cast<internal::BlockOptionsImpl>(opts);

                  auto coriolis = [](const int l, const int m){
                     return (l - MHD_MP(1.0))*(l + MHD_MP(1.0))*precision::sqrt(((l - m)*(l + m))/((MHD_MP(2.0)*l - MHD_MP(1.0))*(MHD_MP(2.0)*l + MHD_MP(1.0))));
                  };

                  const auto T = 1.0/nds.find(NonDimensional::Ekman::id())->second->value();

                  const auto dl = static_cast<MHDFloat>(l);
                  const auto invlapl = 1.0/(dl*(dl + 1.0));
                  SparseSM::Worland::I2Qp corQp(nNr, nNc, o.a, o.b, l, 1*o.truncateQI);
                  auto norm = -coriolis(l+1, o.m);
                  bMat = static_cast<MHDFloat>(norm*T*invlapl)*corQp.mat();
               }

               return bMat;
            };

            // Create first upper diagonal operator
            auto& dUp = getDescription();
            dUp.nRowShift = 0;
            dUp.nColShift = 1;
#ifdef QUICC_INVISCID_2ND
            dUp.realOp = realOpUpper2;
            dUp.imagOp = nullptr;
#else
            dUp.realOp = realOpUpper;
            dUp.imagOp = nullptr;
#endif
         }
         else if(this->useLinearized() &&
               colId == std::make_pair(PhysicalNames::Temperature::id(),FieldComponents::Spectral::SCALAR))
         {

            auto realOp = [](const int nNr, const int nNc, const int l, std::shared_ptr<internal::BlockOptions> opts, const NonDimensional::NdMap& nds)
            {
               SparseMatrix bMat(nNr, nNc);

               if(l > 0)
               {
                  auto& o = *std::dynamic_pointer_cast<internal::BlockOptionsImpl>(opts);

                  const auto Ra = nds.find(NonDimensional::Rayleigh::id())->second->value();
                  const auto T = 1.0/nds.find(NonDimensional::Ekman::id())->second->value();
                  const auto forcing = Ra*T;

                  SparseSM::Worland::I4 i4(nNr, nNc, o.a, o.b, l, 2*o.truncateQI);
                  bMat = -forcing*i4.mat();
               }

               return bMat;
            };

            // Create diagonal block
            auto& d = getDescription();
            d.nRowShift = 0;
            d.nColShift = 0;
            d.realOp = realOp;
            d.imagOp = nullptr;
         }
      }
      else if(rowId == std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR))
      {
         if(rowId == colId)
         {
            // Creat real part of block
            auto realOp = [](const int nNr, const int nNc, const int l, std::shared_ptr<internal::BlockOptions> opts, const NonDimensional::NdMap& nds)
            {
               auto& o = *std::dynamic_pointer_cast<internal::BlockOptionsImpl>(opts);

               const auto Pr = nds.find(NonDimensional::Prandtl::id())->second->value();

               SparseSM::Worland::I2Lapl i2lapl(nNr, nNc, o.a, o.b, l, 1*o.truncateQI);
               SparseMatrix bMat = (1.0/Pr)*i2lapl.mat();

               return bMat;
            };

            // Create diagonal block
            auto& d = getDescription();
            d.nRowShift = 0;
            d.nColShift = 0;
            d.realOp = realOp;
            d.imagOp = nullptr;
         }
         else if(this->useLinearized() &&
               colId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL))
         {
            // Create part of block
            auto realOp = [](const int nNr, const int nNc, const int l, std::shared_ptr<internal::BlockOptions> opts, const NonDimensional::NdMap& nds)
            {
               auto& o = *std::dynamic_pointer_cast<internal::BlockOptionsImpl>(opts);

               const auto dl = static_cast<MHDFloat>(l);
               const auto laplh = (dl*(dl + 1.0));
               SparseSM::Worland::I2 i2(nNr, nNc, o.a, o.b, l, 1*o.truncateQI);
               SparseMatrix bMat = laplh*i2.mat();

               return bMat;
            };

            // Create diagonal block
            auto& d = getDescription();
            d.nRowShift = 0;
            d.nColShift = 0;
            d.realOp = realOp;
            d.imagOp = nullptr;
         }
      }
      else
      {
         //throw std::logic_error("Equations are not setup properly");
      }

      return descr;
   }

   std::vector<internal::BlockDescription> ModelBackend::timeBlockBuilder(const SpectralFieldId& rowId, const SpectralFieldId& colId, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      assert(rowId == colId);
      auto fieldId = rowId;

      std::vector<internal::BlockDescription> descr;

      // Create description with common options
      auto getDescription = [&]()->internal::BlockDescription&
      {
            descr.push_back({});
            auto& d = descr.back();
            auto opts = std::make_shared<internal::BlockOptionsImpl>();
            opts->a = Polynomial::Worland::WorlandBase::ALPHA_CHEBYSHEV;
            opts->b = Polynomial::Worland::WorlandBase::DBETA_CHEBYSHEV;
            opts->m = eigs.at(0);
            opts->bcId = bcs.find(colId.first)->second;
            opts->truncateQI = this->mcTruncateQI;
            opts->isSplitOperator = false;
            d.opts = opts;

            return d;
      };

      if(fieldId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::TOR))
      {
         // Real part of operator
         auto realOp = [](const int nNr, const int nNc, const int l, std::shared_ptr<internal::BlockOptions> opts, const NonDimensional::NdMap& nds)
         {
            assert(nNr == nNc);

            SparseMatrix bMat;
            auto& o = *std::dynamic_pointer_cast<internal::BlockOptionsImpl>(opts);

            if(l > 0)
            {
               SparseSM::Worland::I2 spasm(nNr, nNc, o.a, o.b, l, 1*o.truncateQI);
               bMat = spasm.mat();
            }
            else
            {
               SparseSM::Worland::Id qid(nNr, nNc, o.a, o.b, l);
               bMat = qid.mat();
            }

            return bMat;
         };

         // Real part of operator
         auto realOp2 = [](const int nNr, const int nNc, const int l, std::shared_ptr<internal::BlockOptions> opts, const NonDimensional::NdMap& nds)
         {
            assert(nNr == nNc);

            SparseMatrix bMat;
            auto& o = *std::dynamic_pointer_cast<internal::BlockOptionsImpl>(opts);

            if(l > 0)
            {
               SparseSM::Worland::I2 spasm(nNr+1, nNc, o.a, o.b, l, 1*o.truncateQI);
               int s = 0;
               if(l%2 == 0)
               {
                  s = 1;
               }
               SparseSM::Worland::Id qid(nNr, nNr+1, o.a, o.b, l, 0, s);
               bMat = (qid.mat()*spasm.mat());
            }
            else
            {
               SparseSM::Worland::Id qid(nNr, nNc, o.a, o.b, l);
               bMat = qid.mat();
            }

            return bMat;
         };

         // Create block diagonal operator
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
#ifdef QUICC_INVISCID_2ND
         d.realOp = realOp2;
         d.imagOp = nullptr;
#else
         d.realOp = realOp;
         d.imagOp = nullptr;
#endif
      }
      else if(fieldId == std::make_pair(PhysicalNames::Velocity::id(),FieldComponents::Spectral::POL))
      {
         // Real part of operator
         auto realOp = [](const int nNr, const int nNc, const int l, std::shared_ptr<internal::BlockOptions> opts, const NonDimensional::NdMap& nds)
         {
            assert(nNr == nNc);

            SparseMatrix bMat;
            auto& o = *std::dynamic_pointer_cast<internal::BlockOptionsImpl>(opts);

            if(l > 0)
            {
               SparseSM::Worland::I4Lapl spasm(nNr, nNc, o.a, o.b, l, 2*o.truncateQI);

               // Correct Laplacian for 4th order system according to:
               // McFadden,Murray,Boisvert,
               // Elimination of Spurious Eigenvalues in the
               // Chebyshev Tau Spectral Method,
               // JCP 91, 228-239 (1990)
               // We simply drop the last column
               if(o.bcId == Bc::Name::NoSlip::id())
               {
                  SparseSM::Worland::Id qid(nNr, nNc, o.a, o.b, l, -1);
                  bMat = spasm.mat()*qid.mat();
               }
               else
               {
                  bMat = spasm.mat();
               }
            }
            else
            {
               SparseSM::Worland::Id qid(nNr, nNc, o.a, o.b, l);
               bMat = qid.mat();
            }

            return bMat;
         };

         // Real part of operator
         auto realOp2 = [](const int nNr, const int nNc, const int l, std::shared_ptr<internal::BlockOptions> opts, const NonDimensional::NdMap& nds)
         {
            assert(nNr == nNc);

            SparseMatrix bMat;
            auto& o = *std::dynamic_pointer_cast<internal::BlockOptionsImpl>(opts);

            if(l > 0)
            {
               SparseSM::Worland::I2Lapl spasm(nNr, nNc, o.a, o.b, l, 1*o.truncateQI);
               bMat = spasm.mat();
            }
            else
            {
               SparseSM::Worland::Id qid(nNr, nNc, o.a, o.b, l);
               bMat = qid.mat();
            }

            return bMat;
         };

         // Create block diagonal operator
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
#ifdef QUICC_INVISCID_2ND
         d.realOp = realOp2;
         d.imagOp = nullptr;
#else
         d.realOp = realOp;
         d.imagOp = nullptr;
#endif
      }
      else if(fieldId == std::make_pair(PhysicalNames::Temperature::id(), FieldComponents::Spectral::SCALAR))
      {
         // Real part of operator
         auto realOp = [](const int nNr, const int nNc, const int l, std::shared_ptr<internal::BlockOptions> opts, const NonDimensional::NdMap& nds)
         {
            auto& o = *std::dynamic_pointer_cast<internal::BlockOptionsImpl>(opts);

            SparseSM::Worland::I2 spasm(nNr, nNc, o.a, o.b, l, 1*o.truncateQI);
            SparseMatrix bMat = spasm.mat();

            return bMat;
         };

         // Create block diagonal operator
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = nullptr;
      }

      return descr;
   }

   void ModelBackend::addBlock(SparseMatrix& mat, const SparseMatrix& block, const int rowShift, const int colShift, const MHDFloat coeff) const
   {
      std::vector<Eigen::Triplet<MHDFloat> > triplets;
      triplets.reserve(block.nonZeros());
      for(int k = 0; k < block.outerSize(); ++k)
      {
         for(SparseMatrix::InnerIterator it(block, k); it; ++it)
         {
            triplets.emplace_back(Eigen::Triplet<MHDFloat>(it.row() + rowShift, it.col() + colShift, coeff*it.value()));
         }
      }
      SparseMatrix full(mat.rows(), mat.cols());
      full.setFromTriplets(triplets.begin(), triplets.end());
      mat += full;
   }

   int ModelBackend::blockSize(const SpectralFieldId& fId, const int m, const Resolution& res, const BcMap& bcs, const bool isGalerkin) const
   {
      auto nL = res.counter().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL, m);

      // Compute size
      auto s = 0;
      for(int l = m; l < nL; l++)
      {
         int tN, gN, rhs;
         ArrayI shift(3);
         this->blockInfo(tN, gN, shift, rhs, fId, res, l, bcs);
         if(isGalerkin)
         {
            s += gN;
         }
         else
         {
            s += tN;
         }
      }

      return s;
   }

   std::pair<int, int> ModelBackend::blockShape(const SpectralFieldId& rowId, const SpectralFieldId& colId, const int m, const Resolution& res, const BcMap& bcs, const bool isGalerkin, const bool dropRows) const
   {
      // Compute number of rows
      auto rows = this->blockSize(rowId, m, res, bcs, isGalerkin || dropRows);

      // Compute number of cols
      int cols = this->blockSize(colId, m, res, bcs, isGalerkin);

      return std::make_pair(rows, cols);
   }

   internal::SystemInfo ModelBackend::systemInfo(const SpectralFieldId& rowId, const SpectralFieldId& colId, const int m, const Resolution& res, const BcMap& bcs, const bool isGalerkin, const bool dropRows) const
   {
      auto shape = this->blockShape(rowId, colId, m, res, bcs, isGalerkin, dropRows);

      int sysN = 0;
      bool rowCount = true;
      bool colCount = true;
      int rowIdx = 0;
      int colIdx = 0;
      const auto& fields = this->implicitFields(rowId);
      for(auto it = fields.begin(); it != fields.end(); ++it)
      {
         int s = this->blockSize(*it, m, res, bcs, isGalerkin);
         sysN += s;

         // Get block index of rowId
         if(rowCount && rowId != *it)
         {
            rowIdx += s;
         }
         else if(rowId == *it)
         {
            rowCount = false;
         }

         // Get block index of colId
         if(colCount && colId != *it)
         {
            colIdx += s;
         }
         else if(colId == *it)
         {
            colCount = false;
         }
      }

      internal::SystemInfo info(sysN, shape.first, shape.second, rowIdx, colIdx);
      return info;
   }

   void ModelBackend::modelMatrix(DecoupledZSparse& rModelMatrix, const std::size_t opId, const Equations::CouplingInformation::FieldId_range imRange, const int matIdx, const std::size_t bcType, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      // Time operator
      if(opId == ModelOperator::Time::id())
      {
         for(auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
         {
            auto rowId = *pRowId;
            auto colId = rowId;
            auto descr = timeBlockBuilder(rowId, colId, res, eigs, bcs, nds);
            buildBlock(rModelMatrix, descr, rowId, colId, matIdx, bcType, res, eigs, bcs, nds, false);
         }
      }
      // Linear operator
      else if(opId == ModelOperator::ImplicitLinear::id() || opId == ModelOperator::SplitImplicitLinear::id())
      {
         bool isSplit = (opId == ModelOperator::SplitImplicitLinear::id());

         for(auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
         {
            for(auto pColId = imRange.first; pColId != imRange.second; pColId++)
            {
               auto rowId = *pRowId;
               auto colId = *pColId;
               auto descr = implicitBlockBuilder(rowId, colId, res, eigs, bcs, nds, isSplit);
               buildBlock(rModelMatrix, descr, rowId, colId, matIdx, bcType, res, eigs, bcs, nds, isSplit);
            }
         }
      }
      // Boundary operator
      else if(opId == ModelOperator::Boundary::id() || opId == ModelOperator::SplitBoundary::id())
      {
         bool isSplit = (opId == ModelOperator::SplitBoundary::id());

         for(auto pRowId = imRange.first; pRowId != imRange.second; pRowId++)
         {
            for(auto pColId = imRange.first; pColId != imRange.second; pColId++)
            {
               auto rowId = *pRowId;
               auto colId = *pColId;
               auto descr = boundaryBlockBuilder(rowId, colId, res, eigs, bcs, nds, isSplit);
               buildBlock(rModelMatrix, descr, rowId, colId, matIdx, bcType, res, eigs, bcs, nds, isSplit);
            }
         }
      }
      else
      {
         throw std::logic_error("Requested operator type is not implemented");
      }
   }

   std::vector<internal::BlockDescription> ModelBackend::boundaryBlockBuilder(const SpectralFieldId& rowId, const SpectralFieldId& colId, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds, const bool isSplit) const
   {
      std::vector<internal::BlockDescription> descr;

      // Create description with common options
      auto getDescription = [&]()->internal::BlockDescription&
      {
            descr.push_back({});
            auto& d = descr.back();
            auto opts = std::make_shared<internal::BlockOptionsImpl>();
            opts->a = Polynomial::Worland::WorlandBase::ALPHA_CHEBYSHEV;
            opts->b = Polynomial::Worland::WorlandBase::DBETA_CHEBYSHEV;
            opts->m = eigs.at(0);
            opts->bcId = bcs.find(colId.first)->second;
            opts->truncateQI = this->mcTruncateQI;
            opts->isSplitOperator = isSplit;
            d.opts = opts;

            return d;
      };

      if(rowId == colId)
      {
         // Real part of operator
         auto realOp = [](const int nNr, const int nNc, const int l, std::shared_ptr<internal::BlockOptions> opts, const NonDimensional::NdMap& nds)
         {
            SparseMatrix bMat(nNr, nNc);

            return bMat;
         };

         // Create block diagonal operator
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = realOp;
         d.imagOp = nullptr;
      }
      else
      {
         // Create block diagonal operator
         auto& d = getDescription();
         d.nRowShift = 0;
         d.nColShift = 0;
         d.realOp = nullptr;
         d.imagOp = nullptr;
      }

      return descr;
   }

   void ModelBackend::galerkinStencil(SparseMatrix& mat, const SpectralFieldId& fieldId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const bool makeSquare, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      assert(this->useGalerkin());
      assert(eigs.size() == 1);
      int m = eigs.at(0);

      // Compute system size
      const auto sysRows = systemInfo(fieldId, fieldId, m, res, bcs, makeSquare, false).blockRows;
      const auto sysCols = systemInfo(fieldId, fieldId, m, res, bcs, true, false).blockCols;

      auto nL = res.counter().dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL, m);

      if(mat.size() == 0)
      {
         mat.resize(sysRows,sysCols);
      }
      assert(mat.rows() == sysRows);
      assert(mat.cols() == sysCols);

      assert(eigs.size() == 1);

      int rowShift = 0;
      int colShift = 0;
      for(int l = m; l < nL; l++)
      {
         SparseMatrix S;
         this->stencil(S, fieldId, l, res, makeSquare, bcs, nds);
         this->addBlock(mat, S, rowShift, colShift);

         rowShift += S.rows();
         colShift += S.cols();
      }
   }

   void ModelBackend::explicitBlock(DecoupledZSparse& mat, const SpectralFieldId& fId, const std::size_t opId,  const SpectralFieldId fieldId, const int matIdx, const Resolution& res, const std::vector<MHDFloat>& eigs, const BcMap& bcs, const NonDimensional::NdMap& nds) const
   {
      // Explicit linear operator
      if(opId == ModelOperator::ExplicitLinear::id())
      {
         // Nothing to be done
         throw std::logic_error("There are no explicit linear operators");
      }
      // Explicit nonlinear operator
      else if(opId == ModelOperator::ExplicitNonlinear::id())
      {
         throw std::logic_error("There are no explicit nonlinear operators");
      }
      // Explicit nextstep operator
      else if(opId == ModelOperator::ExplicitNextstep::id())
      {
         throw std::logic_error("There are no explicit nextstep operators");
      }
   }

} // Implicit
} // RTC
} // Sphere
} // Boussinesq
} // Model
} // QuICC
