/**
 * @file IRTCBackend.cpp
 * @brief Source of the interface for model backend
 */

// System includes
//
#include <stdexcept>

// Project includes
//
#include "QuICC/Model/Boussinesq/Sphere/RTC/IRTCBackend.hpp"
#include "QuICC/ModelOperator/Time.hpp"
#include "QuICC/ModelOperator/ImplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitLinear.hpp"
#include "QuICC/ModelOperator/ExplicitNonlinear.hpp"
#include "QuICC/ModelOperator/ExplicitNextstep.hpp"
#include "QuICC/ModelOperator/Stencil.hpp"
#include "QuICC/ModelOperator/Boundary.hpp"
#include "QuICC/ModelOperator/SplitImplicitLinear.hpp"
#include "QuICC/ModelOperator/SplitBoundary.hpp"
#include "QuICC/ModelOperatorBoundary/FieldToRhs.hpp"
#include "QuICC/ModelOperatorBoundary/SolverHasBc.hpp"
#include "QuICC/ModelOperatorBoundary/SolverNoTau.hpp"
#include "QuICC/ModelOperatorBoundary/Stencil.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Bc/Name/FixedTemperature.hpp"
#include "QuICC/Bc/Name/FixedFlux.hpp"
#include "QuICC/Bc/Name/StressFree.hpp"
#include "QuICC/Bc/Name/NoSlip.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/NonDimensional/Prandtl.hpp"
#include "QuICC/NonDimensional/Rayleigh.hpp"
#include "QuICC/NonDimensional/Ekman.hpp"
#include "QuICC/NonDimensional/CflInertial.hpp"
#include "QuICC/Tools/IdToHuman.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace RTC {

   IRTCBackend::IRTCBackend()
      : IModelBackend(), mUseGalerkin(false), mUseSplitEquation(false)
   {
   }

   std::vector<std::string> IRTCBackend::fieldNames() const
   {
      std::vector<std::string> names = {
         PhysicalNames::Velocity().tag(),
         PhysicalNames::Temperature().tag()
      };

      return names;
   }

   std::vector<std::string> IRTCBackend::paramNames() const
   {
      std::vector<std::string> names = {
         NonDimensional::Prandtl().tag(),
         NonDimensional::Rayleigh().tag(),
         NonDimensional::Ekman().tag()
      };

      return names;
   }

   std::vector<bool> IRTCBackend::isPeriodicBox() const
   {
      std::vector<bool> periodic = {false, false, false};

      return periodic;
   }

   bool IRTCBackend::useGalerkin() const
   {
      return this->mUseGalerkin;
   }

   void IRTCBackend::enableGalerkin(const bool flag)
   {
      this->mUseGalerkin = flag;
   }

   bool IRTCBackend::useSplitEquation() const
   {
      return this->mUseSplitEquation;
   }

   void IRTCBackend::enableSplitEquation(const bool tag)
   {
      this->mUseSplitEquation = tag;
   }

   std::map<std::string,MHDFloat> IRTCBackend::automaticParameters(const std::map<std::string,MHDFloat>& cfg) const
   {
      auto E = cfg.find(NonDimensional::Ekman().tag())->second;

      std::map<std::string,MHDFloat> params = {
         {NonDimensional::CflInertial().tag(), 0.1*E}
      };

      return params;
   }

} // RTC
} // Sphere
} // Boussinesq
} // Model
} // QuICC
