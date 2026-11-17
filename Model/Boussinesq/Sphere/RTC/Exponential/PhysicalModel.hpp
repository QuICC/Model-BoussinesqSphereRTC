/**
 * @file PhysicalModel.hpp
 * @brief Implementation of the Boussinesq rotating thermal convection in a
 * sphere (Toroidal/Poloidal formulation) without coupled solve (standard
 * implementation)
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_EXPONENTIAL_PHYSICALMODEL_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_EXPONENTIAL_PHYSICALMODEL_HPP

// System includes
//
#include <string>

// Project includes
//
#include "QuICC/SpatialScheme/3D/WLFl.hpp"
#include "Model/Boussinesq/Sphere/RTC/Exponential/ModelBackend.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace RTC {

namespace Exponential {

/**
 * @brief Implementation of the Boussinesq rotating thermal convection sphere
 * model (Toroidal/Poloidal formulation) without coupled solve (standard
 * implementation)
 */
template <typename TBuilder> class PhysicalModel : public TBuilder
{
public:
   /// Typedef for the spatial scheme used
   typedef SpatialScheme::WLFl SchemeType;

   /**
    * @brief Constructor
    */
   PhysicalModel() = default;

   /**
    * @brief Destructor
    */
   virtual ~PhysicalModel() = default;

   /**
    * @brief Initialize specialized backend
    */
   void init() final;

protected:
private:
};

template <typename TBuilder> void PhysicalModel<TBuilder>::init()
{
   TBuilder::init();
#ifdef QUICC_MODEL_BOUSSINESQSPHERERTC_EXPONENTIAL_BACKEND_CPP

   this->mpBackend = std::make_shared<ModelBackend>();
#else
   throw std::logic_error("Python backend is not implemented");
#endif
}

} // namespace Exponential
} // namespace RTC
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_EXPONENTIAL_PHYSICALMODEL_HPP
