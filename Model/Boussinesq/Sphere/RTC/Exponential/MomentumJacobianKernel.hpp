/**
 * @file MomentumJacobianKernel.hpp
 * @brief Physical kernel for the Momentum nonlinear kernel
 */

#ifndef QUICC_PHYSICAL_MOMENTUMJACOBIANKERNEL_HPP
#define QUICC_PHYSICAL_MOMENTUMJACOBIANKERNEL_HPP

// System includes
//
#include <memory>

// Project includes
//
#include "QuICC/PhysicalKernels/IPhysicalKernel.hpp"

namespace QuICC {

namespace Physical {

namespace Kernel {

/**
 * @brief Physical kernel for the Momentum nonlinear kernel
 */
class MomentumJacobianKernel : public IPhysicalKernel
{
public:
   /**
    * @brief Simple constructor
    */
   MomentumJacobianKernel() = default;

   /**
    * @brief Simple empty destructor
    */
   ~MomentumJacobianKernel() = default;

   /**
    * @brief Set the physical mesh on which kernel is working
    */
   virtual void setMesh(std::shared_ptr<std::vector<Array>> spMesh) override;

   /**
    * @brief Set the smart pointer to the velocity field
    *
    * \param name Name of the field
    * \param spField Shared pointer to the vector field
    */
   void setJacobianVelocity(std::size_t name,
      Framework::Selector::VariantSharedVectorVariable spField);

   /**
    * @brief Set the smart pointer to the velocity field
    *
    * \param name Name of the field
    * \param spField Shared pointer to the vector field
    */
   void setVelocity(std::size_t name,
      Framework::Selector::VariantSharedVectorVariable spField);

   /**
    * @brief Set the smart pointer to the temperature field
    *
    * \param name Name of the field
    * \param spField Shared pointer to the scalar field
    */
   void setJacobianTemperature(std::size_t name,
      Framework::Selector::VariantSharedScalarVariable spField);

   /**
    * @brief Set the smart pointer to the temperature field
    *
    * \param name Name of the field
    * \param spField Shared pointer to the scalar field
    */
   void setTemperature(std::size_t name,
      Framework::Selector::VariantSharedScalarVariable spField);

   /**
    * @brief Initialize kernel
    */
   void init(const MHDFloat inertia, const MHDFloat coriolis,
      const MHDFloat buoyancy);

   /**
    * @brief Compute the physical kernel
    *
    * @param rNLComp Nonlinear term component
    * @param id      ID of the component (allows for a more general
    * implementation)
    */
   virtual void compute(Framework::Selector::PhysicalScalarField& rNLComp,
      FieldComponents::Physical::Id id) const override;

protected:
   /**
    * @brief Get name ID of the velocity
    */
   std::size_t name() const;

private:
   /**
    * @brief Name ID of the jacobian velocity
    */
   std::size_t mJVName;

   /**
    * @brief Name ID of the velocity
    */
   std::size_t mVName;

   /**
    * @brief Name ID of the Jacobian temperature
    */
   std::size_t mJTName;

   /**
    * @brief Name ID of the temperature
    */
   std::size_t mTName;

   /**
    * @brief Scaling constant for inertial term
    */
   MHDFloat mInertia;

   /**
    * @brief Scaling constant for Coriolis term
    */
   MHDFloat mCoriolis;

   /**
    * @brief Scaling constant for Buoyancy term
    */
   MHDFloat mBuoyancy;

   /**
    * @brief Storage for the radial(theta) grid values (if required)
    */
   Array mRadius;

   /**
    * @brief Storage for the cos(theta) grid values (if required)
    */
   Array mCosTheta;

   /**
    * @brief Storage for the sin(theta) grid values (if required)
    */
   Array mSinTheta;
};

/// Typedef for a smart MomentumJacobianKernel
typedef std::shared_ptr<MomentumJacobianKernel> SharedMomentumJacobianKernel;

} // namespace Kernel
} // namespace Physical
} // namespace QuICC

#endif // QUICC_PHYSICAL_MOMENTUMJACOBIANKERNEL_HPP
