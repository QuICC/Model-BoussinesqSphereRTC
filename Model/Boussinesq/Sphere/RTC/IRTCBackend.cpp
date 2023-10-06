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
#include "QuICC/Polynomial/Worland/WorlandBase.hpp"
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
        nBc = 1;
    }
    else if (fId == std::make_pair(PhysicalNames::Velocity::id(),
                       FieldComponents::Spectral::POL))
    {
        nBc = 2;
    }
    else
    {
        nBc = 0;
    }

    return nBc;
}

void IRTCBackend::blockInfo(int& tN, int& gN, ArrayI& shift, int& rhs,
   const SpectralFieldId& fId, const Resolution& res, const MHDFloat l,
   const BcMap& bcs) const
{
    auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
    tN = nN;

    int shiftR = this->nBc(fId);
    if (this->useGalerkin())
    {
        gN = (nN - shiftR);
    }
    else
    {
        shiftR = 0;
        gN = nN;
    }

    // Set galerkin shifts
    shift(0) = shiftR;
    shift(1) = 0;
    shift(2) = 0;

    rhs = 1;
}

void IRTCBackend::applyTau(SparseMatrix& mat, const SpectralFieldId& rowId,
   const SpectralFieldId& colId, const int l, const Resolution& res,
   const BcMap& bcs, const NonDimensional::NdMap& nds,
   const bool isSplitOperator) const
{
    auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);

    auto a = Polynomial::Worland::WorlandBase::ALPHA_CHEBYSHEV;
    auto b = Polynomial::Worland::WorlandBase::DBETA_CHEBYSHEV;

    auto bcId = bcs.find(rowId.first)->second;

    SparseSM::Worland::Boundary::Operator bcOp(nN, nN, a, b, l);

    if (rowId == std::make_pair(PhysicalNames::Velocity::id(),
                    FieldComponents::Spectral::TOR) &&
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
    else if (rowId == std::make_pair(PhysicalNames::Velocity::id(),
                         FieldComponents::Spectral::POL) &&
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
    else if (rowId == std::make_pair(PhysicalNames::Temperature::id(),
                         FieldComponents::Spectral::SCALAR) &&
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
   const int l, const Resolution& res, const bool makeSquare, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
    auto nN = res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);

    auto a = Polynomial::Worland::WorlandBase::ALPHA_CHEBYSHEV;
    auto b = Polynomial::Worland::WorlandBase::DBETA_CHEBYSHEV;

    auto bcId = bcs.find(fieldId.first)->second;

    int s = this->nBc(fieldId);
    if (fieldId == std::make_pair(PhysicalNames::Velocity::id(),
                      FieldComponents::Spectral::TOR))
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
    else if (fieldId == std::make_pair(PhysicalNames::Velocity::id(),
                           FieldComponents::Spectral::POL))
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
    else if (fieldId == std::make_pair(PhysicalNames::Temperature::id(),
                           FieldComponents::Spectral::SCALAR))
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

    if (makeSquare)
    {
        SparseSM::Worland::Id qId(nN - s, nN, a, b, l);
        mat = qId.mat() * mat;
    }
}

void IRTCBackend::applyGalerkinStencil(SparseMatrix& mat,
   const SpectralFieldId& rowId, const SpectralFieldId& colId, const int lr,
   const int lc, const Resolution& res, const BcMap& bcs,
   const NonDimensional::NdMap& nds) const
{
    auto nNr = res.counter().dimensions(Dimensions::Space::SPECTRAL, lr)(0);

    auto a = Polynomial::Worland::WorlandBase::ALPHA_CHEBYSHEV;
    auto b = Polynomial::Worland::WorlandBase::DBETA_CHEBYSHEV;

    auto S = mat;
    this->stencil(S, colId, lc, res, false, bcs, nds);

    auto s = this->nBc(rowId);
    SparseSM::Worland::Id qId(nNr - s, nNr, a, b, lr, 0, s);
    mat = qId.mat() * (mat * S);
}

void IRTCBackend::buildBlock(DecoupledZSparse& decMat,
   const std::vector<details::BlockDescription>& descr,
   const SpectralFieldId& rowId, const SpectralFieldId& colId, const int matIdx,
   const std::size_t bcType, const Resolution& res, const int l0,
   const int maxL, const BcMap& bcs, const NonDimensional::NdMap& nds,
   const bool isSplitOperator) const
{
    // Compute system size
    const auto sysInfo =
       systemInfo(rowId, colId, l0, maxL, res, bcs, this->useGalerkin(), false);
    const auto& sysN = sysInfo.systemSize;
    const auto& baseRowShift = sysInfo.startRow;
    const auto& baseColShift = sysInfo.startCol;

    // Resize matrices the first time
    if (decMat.real().size() == 0)
    {
        decMat.real().resize(sysN, sysN);
        if (this->isComplex(rowId))
        {
            decMat.imag().resize(sysN, sysN);
        }
    }
    assert(decMat.real().rows() == sysN);
    assert(decMat.real().cols() == sysN);
    if (this->isComplex(rowId))
    {
        assert(decMat.imag().rows() == sysN);
        assert(decMat.imag().cols() == sysN);
    }

    int tN, gN, rhs;
    ArrayI shift(3);

    bool needStencil = (this->useGalerkin() &&
                        bcType == ModelOperatorBoundary::SolverNoTau::id());
    bool needTau = (bcType == ModelOperatorBoundary::SolverHasBc::id());

    for (auto&& d: descr)
    {
        assert(d.nRowShift == 0 || d.nColShift == 0);

        // Shift starting row
        int rowShift = baseRowShift;
        for (int s = 0; s < d.nRowShift; s++)
        {
            this->blockInfo(tN, gN, shift, rhs, rowId, res, l0 + s, bcs);
            rowShift += gN;
        }

        // Shift starting col
        int colShift = baseColShift;
        for (int s = 0; s < d.nColShift; s++)
        {
            this->blockInfo(tN, gN, shift, rhs, colId, res, l0 + s, bcs);
            colShift += gN;
        }

        int lShift = -d.nRowShift + d.nColShift;

        for (int l = l0 + d.nRowShift; l <= maxL - d.nColShift; l++)
        {
            auto nNr =
               res.counter().dimensions(Dimensions::Space::SPECTRAL, l)(0);
            auto nNc = res.counter().dimensions(Dimensions::Space::SPECTRAL,
               l + lShift)(0);

            //
            // Build real part of block
            if (d.realOp)
            {
                auto bMat = d.realOp(nNr, nNc, l, d.opts, nds);

                if (needStencil)
                {
                    this->applyGalerkinStencil(bMat, rowId, colId, l,
                       l + lShift, res, bcs, nds);
                }
                else if (needTau)
                {
                    this->applyTau(bMat, rowId, colId, l + lShift, res, bcs,
                       nds, isSplitOperator);
                }
                this->addBlock(decMat.real(), bMat, rowShift, colShift);
            }

            //
            // Build imaginary part of block
            if (d.imagOp)
            {
                auto bMat = d.imagOp(nNr, nNc, l, d.opts, nds);

                if (needStencil)
                {
                    this->applyGalerkinStencil(bMat, rowId, colId, l,
                       l + lShift, res, bcs, nds);
                }
                else if (needTau)
                {
                    this->applyTau(bMat, rowId, colId, l + lShift, res, bcs,
                       nds, isSplitOperator);
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

void IRTCBackend::addBlock(SparseMatrix& mat, const SparseMatrix& block,
   const int rowShift, const int colShift, const MHDFloat coeff) const
{
    std::vector<Eigen::Triplet<MHDFloat>> triplets;
    triplets.reserve(block.nonZeros());
    for (int k = 0; k < block.outerSize(); ++k)
    {
        for (SparseMatrix::InnerIterator it(block, k); it; ++it)
        {
            triplets.emplace_back(Eigen::Triplet<MHDFloat>(it.row() + rowShift,
               it.col() + colShift, coeff * it.value()));
        }
    }
    SparseMatrix full(mat.rows(), mat.cols());
    full.setFromTriplets(triplets.begin(), triplets.end());
    mat += full;
}

int IRTCBackend::blockSize(const SpectralFieldId& fId, const int l0,
   const int maxL, const Resolution& res, const BcMap& bcs,
   const bool isGalerkin) const
{
    // Compute size
    auto s = 0;
    for (int l = l0; l <= maxL; l++)
    {
        int tN, gN, rhs;
        ArrayI shift(3);
        this->blockInfo(tN, gN, shift, rhs, fId, res, l, bcs);
        if (isGalerkin)
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

std::pair<int, int> IRTCBackend::blockShape(const SpectralFieldId& rowId,
   const SpectralFieldId& colId, const int l0, const int maxL,
   const Resolution& res, const BcMap& bcs, const bool isGalerkin,
   const bool dropRows) const
{
    // Compute number of rows
    auto rows =
       this->blockSize(rowId, l0, maxL, res, bcs, isGalerkin || dropRows);

    // Compute number of cols
    int cols = this->blockSize(colId, l0, maxL, res, bcs, isGalerkin);

    return std::make_pair(rows, cols);
}

details::SystemInfo IRTCBackend::systemInfo(const SpectralFieldId& rowId,
   const SpectralFieldId& colId, const int l0, const int maxL,
   const Resolution& res, const BcMap& bcs, const bool isGalerkin,
   const bool dropRows) const
{
    auto shape =
       this->blockShape(rowId, colId, l0, maxL, res, bcs, isGalerkin, dropRows);

    int sysN = 0;
    bool rowCount = true;
    bool colCount = true;
    int rowIdx = 0;
    int colIdx = 0;
    const auto& fields = this->implicitFields(rowId);
    for (auto it = fields.begin(); it != fields.end(); ++it)
    {
        int s = this->blockSize(*it, l0, maxL, res, bcs, isGalerkin);
        sysN += s;

        // Get block index of rowId
        if (rowCount && rowId != *it)
        {
            rowIdx += s;
        }
        else if (rowId == *it)
        {
            rowCount = false;
        }

        // Get block index of colId
        if (colCount && colId != *it)
        {
            colIdx += s;
        }
        else if (colId == *it)
        {
            colCount = false;
        }
    }

    details::SystemInfo info(sysN, shape.first, shape.second, rowIdx, colIdx);
    return info;
}

} // namespace RTC
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
