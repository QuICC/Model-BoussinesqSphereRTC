/**
 * @file PhysicalModel.cpp
 * @brief Source of the Boussinesq rotating thermal convection in a sphere (Toroidal/Poloidal formulation)
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Model/Boussinesq/Sphere/RTC/Implicit/PhysicalModel.hpp"

// Project includes
//

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace RTC {

namespace Implicit {

   std::string PhysicalModel::PYMODULE()
   {
      return "boussinesq.sphere.rtc.implicit.physical_model";
   }

}
}
}
}
}
}
