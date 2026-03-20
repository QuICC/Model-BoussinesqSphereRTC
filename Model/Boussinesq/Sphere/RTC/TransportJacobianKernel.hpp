/**
 * @file TransportJacobianKernel.hpp
 * @brief Physical kernel for the Transport nonlinear Jacobian
 */

#ifndef QUICC_PHYSICAL_TRANSPORTJACOBIANKERNEL_HPP
#define QUICC_PHYSICAL_TRANSPORTJACOBIANKERNEL_HPP

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
 * @brief Physical kernel for the Momentum nonlinear Jacobian
 */
class TransportJacobianKernel : public IPhysicalKernel
{
public:
   /**
    * @brief Simple constructor
    */
   TransportJacobianKernel() = default;

   /**
    * @brief Simple empty destructor
    */
   ~TransportJacobianKernel() = default;

   /**
    * @brief Set the physical mesh on which kernel is working
    */
   virtual void setMesh(std::shared_ptr<std::vector<Array>> spMesh) override;

   /**
    * @brief Set the smart pointer to the scalar field
    *
    * \param name Name of the field
    * \param spField Shared pointer to the scalar field
    */
   void setJacobianScalar(std::size_t name,
      Framework::Selector::VariantSharedScalarVariable spField);

   /**
    * @brief Set the smart pointer to the scalar field
    *
    * \param name Name of the field
    * \param spField Shared pointer to the scalar field
    */
   void setScalar(std::size_t name,
      Framework::Selector::VariantSharedScalarVariable spField);

   /**
    * @brief Set the smart pointer to the vector field
    *
    * \param name Name of the field
    * \param spField Shared pointer to the vector field
    */
   void setJacobianVector(std::size_t name,
      Framework::Selector::VariantSharedVectorVariable spField);

   /**
    * @brief Set the smart pointer to the vector field
    *
    * \param name Name of the field
    * \param spField Shared pointer to the vector field
    */
   void setVector(std::size_t name,
      Framework::Selector::VariantSharedVectorVariable spField);

   /**
    * @brief Initialize kernel
    */
   void init(const MHDFloat transport);

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
    * @brief Get name ID of the unknown
    */
   std::size_t name() const;

private:
   /**
    * @brief Name ID of the Jacobian scalar field
    */
   std::size_t mJName;

   /**
    * @brief Name ID of the scalar field
    */
   std::size_t mName;

   /**
    * @brief Name ID of the Jacobian vector field
    */
   std::size_t mJVName;

   /**
    * @brief Name ID of the vector field
    */
   std::size_t mVName;

   /**
    * @brief Scaling constant for transport term
    */
   MHDFloat mTransport;

   /**
    * @brief Storage for the radial grid values
    */
   Array mRadius;
};

/// Typedef for a smart TransportJacobianKernel
typedef std::shared_ptr<TransportJacobianKernel> SharedTransportJacobian;

} // namespace Kernel
} // namespace Physical
} // namespace QuICC

#endif // QUICC_PHYSICAL_TRANSPORTJACOBIANKERNEL_HPP
