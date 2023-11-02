/**
 * @file PhysicalModel.hpp
 * @brief Implementation of the Boussinesq rotating thermal convection in a
 * sphere (Toroidal/Poloidal formulation) without coupled solve (standard
 * implementation)
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_EXPLICIT_PHYSICALMODEL_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_EXPLICIT_PHYSICALMODEL_HPP

// System includes
//
#include <string>

// Project includes
//
#include "Model/Boussinesq/Sphere/RTC/IRTCModel.hpp"
#include "QuICC/SpatialScheme/3D/WLFl.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace RTC {

namespace Explicit {

/**
 * @brief Implementation of the Boussinesq rotating thermal convection sphere
 * model (Toroidal/Poloidal formulation) without coupled solve (standard
 * implementation)
 */
class PhysicalModel : public IRTCModel
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

   /// Python script/module name
   std::string PYMODULE() final;

   /**
    * @brief Initialize specialized backend
    */
   void init() final;

protected:
private:
};

} // namespace Explicit
} // namespace RTC
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_EXPLICIT_PHYSICALMODEL_HPP
