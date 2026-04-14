/**
 * @file IExponentialRTCModel.hpp
 * @brief Implementation of the Boussinesq rotating thermal convection in a
 * sphere (Toroidal/Poloidal formulation)
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_EXPONENTIAL_IEXPONENTIALRTCMODEL_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_EXPONENTIAL_IEXPONENTIALRTCMODEL_HPP

// System includes
//
#include <string>

// Project includes
//
#include "QuICC/Simulation/Simulation.hpp"
#include "Model/Boussinesq/Sphere/RTC/IRTCModel.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace RTC {

namespace Exponential {

/**
 * @brief Implementation of the Boussinesq rotating thermal convection sphere
 * model (Toroidal/Poloidal formulation)
 */
class IExponentialRTCModel : public IRTCModel
{
public:
   /**
    * @brief Constructor
    */
   IExponentialRTCModel() = default;

   /**
    * @brief Destructor
    */
   virtual ~IExponentialRTCModel() = default;

   /**
    * @brief Exclude fields from initial state
    */
   virtual std::vector<std::size_t> excludedFieldIds() const override;

   /**
    * @brief Add the required equations
    *
    * @param spSim   Shared simulation object
    */
   virtual void addEquations(SharedSimulation spSim) override;

protected:
private:
};

} // namespace Exponential
} // namespace RTC
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_EXPONENTIAL_IEXPONENTIALRTCMODEL_HPP
