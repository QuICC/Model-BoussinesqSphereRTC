/**
 * @file IRTCState.cpp
 * @brief Source of the Boussinesq rotating thermal convection in a sphere
 * (Toroidal/Poloidal formulation)
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Sphere/RTC/IRTCState.hpp"
#include "Model/Boussinesq/Sphere/RTC/gitHash.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Generator/States/Kernels/Sphere/BenchmarkTempC1.hpp"
#include "QuICC/Generator/States/Kernels/Sphere/ScalarYllPerturbation.hpp"
#include "QuICC/Generator/States/SphereExactScalarState.hpp"
#include "QuICC/Generator/States/SphereExactVectorState.hpp"
#include "QuICC/Io/Variable/StateFileWriter.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/SpectralKernels/MakeRandom.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace RTC {

VectorFormulation::Id IRTCState::SchemeFormulation()
{
   return VectorFormulation::TORPOL;
}

std::string IRTCState::version() const
{
   return std::string(gitHash);
}

void IRTCState::addStates(SharedStateGenerator spGen)
{
   // Shared pointer to equation
   Equations::SharedSphereExactScalarState spScalar;
   Equations::SharedSphereExactVectorState spVector;

   Spectral::Kernel::Complex3DMapType tSH;
   std::pair<Spectral::Kernel::Complex3DMapType::iterator, bool> ptSH;

   // Add temperature initial state generator
   spScalar =
      spGen->addEquation<Equations::SphereExactScalarState>(this->spBackend());
   spScalar->setIdentity(PhysicalNames::Temperature::id());
   switch (3)
   {
   case 0: {
      spScalar->setPhysicalNoise(1e-15);
   }
   break;

   case 1: {
      spScalar->setPhysicalConstant(1.0);
   }
   break;

   case 2: {
      tSH.clear();
      ptSH = tSH.insert(
         std::make_pair(std::make_pair(3, 3), std::map<int, MHDComplex>()));
      ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0, 2.0)));
      spScalar->setSpectralModes(tSH);
   }
   break;

   case 3: {
      auto spKernel =
         std::make_shared<Physical::Kernel::Sphere::BenchmarkTempC1>();
      spKernel->init(0.0, 1e-5);
      spScalar->setPhysicalKernel(spKernel);
   }
   break;

   case 4: {
      auto spKernel = std::make_shared<Spectral::Kernel::MakeRandom>(
         spGen->ss().has(SpatialScheme::Feature::ComplexSpectrum));
      std::vector<MHDFloat> ratios = {1e4, 1e4, 1e4};
      spKernel->setRatio(ratios);
      spKernel->init(-1e-4, 1e-4);
      spScalar->setSrcKernel(spKernel);
   }
   break;

   case 5: {
      auto spKernel =
         std::make_shared<Physical::Kernel::Sphere::ScalarYllPerturbation>();
      const MHDFloat amplitude_bg = 0.0;
      const MHDFloat eps = 1e-5;
      const int m = 3;
      spKernel->init(amplitude_bg, eps, m);
      spScalar->setPhysicalKernel(spKernel);
   }
   break;
   }

   // Add velocity initial state generator
   spVector =
      spGen->addEquation<Equations::SphereExactVectorState>(this->spBackend());
   spVector->setIdentity(PhysicalNames::Velocity::id());
   switch (3)
   {
   // Toroidal only
   case 0: {
      // Toroidal
      tSH.clear();
      ptSH = tSH.insert(
         std::make_pair(std::make_pair(1, 1), std::map<int, MHDComplex>()));
      ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
      spVector->setSpectralModes(FieldComponents::Spectral::TOR, tSH);
      // Poloidal
      tSH.clear();
      spVector->setSpectralModes(FieldComponents::Spectral::POL, tSH);
   }
   break;

   // Poloidal only
   case 1: {
      // Toroidal
      tSH.clear();
      spVector->setSpectralModes(FieldComponents::Spectral::TOR, tSH);
      // Poloidal
      tSH.clear();
      ptSH = tSH.insert(
         std::make_pair(std::make_pair(2, 0), std::map<int, MHDComplex>()));
      ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
      spVector->setSpectralModes(FieldComponents::Spectral::POL, tSH);
   }
   break;

   // Toroidal & Poloidal
   case 2: {
      // Toroidal
      tSH.clear();
      ptSH = tSH.insert(
         std::make_pair(std::make_pair(1, 1), std::map<int, MHDComplex>()));
      ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
      spVector->setSpectralModes(FieldComponents::Spectral::TOR, tSH);
      // Poloidal
      tSH.clear();
      ptSH = tSH.insert(
         std::make_pair(std::make_pair(2, 0), std::map<int, MHDComplex>()));
      ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
      spVector->setSpectralModes(FieldComponents::Spectral::POL, tSH);
   }
   break;

   case 3: {
      auto spKernel = std::make_shared<Spectral::Kernel::MakeRandom>(
         spGen->ss().has(SpatialScheme::Feature::ComplexSpectrum));
      std::vector<MHDFloat> ratios = {1e2, 1e2, 1e2};
      spKernel->setRatio(ratios);
      spKernel->init(-1e-15, 1e-15);
      spVector->setSrcKernel(FieldComponents::Spectral::TOR, spKernel);
      spVector->setSrcKernel(FieldComponents::Spectral::POL, spKernel);
   }
   break;

   case 4: {
      auto spKernel = std::make_shared<Spectral::Kernel::MakeRandom>(
         spGen->ss().has(SpatialScheme::Feature::ComplexSpectrum));
      std::vector<MHDFloat> ratios = {1e4, 1e4, 1e4};
      spKernel->setRatio(ratios);
      spKernel->init(-1e-4, 1e-4);
      spVector->setSrcKernel(FieldComponents::Spectral::TOR, spKernel);
      spVector->setSrcKernel(FieldComponents::Spectral::POL, spKernel);
   }
   break;
   }

   // Add output file
   auto spOut =
      std::make_shared<Io::Variable::StateFileWriter>(spGen->ss().tag(),
         spGen->ss().has(SpatialScheme::Feature::RegularSpectrum));
   spOut->expect(PhysicalNames::Temperature::id());
   spOut->expect(PhysicalNames::Velocity::id());
   spGen->addHdf5OutputFile(spOut);
}

} // namespace RTC
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
