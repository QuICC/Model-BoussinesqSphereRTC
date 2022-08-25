/**
 * @file PhysicalModel.hpp
 * @brief Implementation of the Boussinesq rotating thermal convection in a sphere (Toroidal/Poloidal formulation)
 */

#ifndef QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_IMPLICIT_PHYSICALMODEL_HPP
#define QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_IMPLICIT_PHYSICALMODEL_HPP

// Model version
#define QUICC_VERSION_MODEL_MAJOR 1
#define QUICC_VERSION_MODEL_MINOR 0
#define QUICC_VERSION_MODEL_PATCH 0

// Configuration includes
//

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "QuICC/Model/Boussinesq/Sphere/RTC/IRTCModel.hpp"
#include "QuICC/SpatialScheme/3D/WLFm.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace RTC {

namespace Implicit {

   /**
    * @brief Implementation of the Boussinesq rotating thermal convection sphere model (Toroidal/Poloidal formulation)
    */
   class PhysicalModel: public IRTCModel
   {
      public:
         /// Typedef for the spatial scheme used
         typedef SpatialScheme::WLFm SchemeType;

         /**
          * @brief Constructor
          */
         PhysicalModel() = default;

         /**
          * @brief Destructor
          */
         virtual ~PhysicalModel() = default;

         /// Python script/module name
         virtual std::string PYMODULE() override;

      protected:

      private:
   };

}
}
}
}
}
}

#endif // QUICC_MODEL_BOUSSINESQ_SPHERE_RTC_IMPLICIT_PHYSICALMODEL_HPP
