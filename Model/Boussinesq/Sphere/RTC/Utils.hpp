/**
 * @file Utils.hpp
 * @brief Implementation utilities
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_UTILS_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_UTILS_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Simulation/SimulationBoundary.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/SparseSM/Bessel/BesselKind.hpp"
#include "Polynomial/SphericalBessel/Types.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace RTC {
   /**
    * @brief Get transform path
    *
    * @param spBcs   Simulation boundary conditions
    * @param fieldId Field ID
    */
   std::size_t getPathId(std::shared_ptr<SimulationBoundary> spBcs, const std::size_t fieldId,  const bool isNl, const FieldComponents::Spectral::Id comp = FieldComponents::Spectral::NOTUSED);

   /**
    * @brief Kind of Bessel basis
    *
    * @fId  Field ID
    */
   SparseSM::Bessel::BesselKind bKind(const SpectralFieldId& fId);

   /**
    * @brief Convert between SparseSM Bessel kind and Bessel basis
    *
    * @para kind  Bessel kind
    */
   Polynomial::SphericalBessel::sphjnl_basis_t bBasis(const SparseSM::Bessel::BesselKind kind);

} // namespace RTC
} // namespace Sphere
} // namespace Boussinesq
} // namespace Equations
} // namespace QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_UTILS_HPP
