/**
 * @file Momentum.hpp
 * @brief Implementation of the vector Navier-Stokes equation for the Boussinesq
 * rotating thermal convection in a sphere
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_MOMENTUM_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_MOMENTUM_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/Equations/IVectorEquation.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Sphere {

namespace RTC {

/**
 * @brief Implementation of the vector Navier-Stokes equation for the Boussinesq
 * rotating thermal convection in a sphere
 */
class Momentum : public IVectorEquation
{
public:
   /**
    * @brief Simple constructor
    *
    * @param spEqParams  Shared equation parameters
    */
   Momentum(SharedEquationParameters spEqParams,
      SpatialScheme::SharedCISpatialScheme spScheme,
      std::shared_ptr<Model::IModelBackend> spBackend);

   /**
    * @brief Simple empty destructor
    */
   ~Momentum() = default;

   /**
    * @brief Initialize constraint kernel
    *
    * @param spMesh  Physical space mesh
    */
   void initConstraintKernel(
      const std::shared_ptr<std::vector<Array>> spMesh) final;

   /**
    * @brief Initialize nonlinear interaction kernel
    */
   void initNLKernel(const bool force = false) final;

protected:
   /**
    * @brief Set variable requirements
    */
   void setRequirements() final;

   /**
    * @brief Set the equation coupling information
    */
   void setCoupling() final;

   /**
    * @brief Set the nonlinear integration components
    */
   void setNLComponents() final;

private:
};

} // namespace RTC
} // namespace Sphere
} // namespace Boussinesq
} // namespace Equations
} // namespace QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_MOMENTUM_HPP
