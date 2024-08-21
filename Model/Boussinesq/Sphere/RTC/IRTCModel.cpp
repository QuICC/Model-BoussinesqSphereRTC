/**
 * @file IRTCModel.cpp
 * @brief Source of the Boussinesq rotating thermal convection in a sphere
 * (Toroidal/Poloidal formulation)
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Sphere/RTC/IRTCModel.hpp"
#include "Model/Boussinesq/Sphere/RTC/Momentum.hpp"
#include "Model/Boussinesq/Sphere/RTC/Transport.hpp"
#include "Model/Boussinesq/Sphere/RTC/gitHash.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/Generator/States/Kernels/Sphere/BenchmarkTempC1.hpp"
#include "QuICC/Generator/States/Kernels/Sphere/ScalarYllPerturbation.hpp"
#include "QuICC/Generator/States/RandomScalarState.hpp"
#include "QuICC/Generator/States/RandomVectorState.hpp"
#include "QuICC/Generator/States/SphereExactScalarState.hpp"
#include "QuICC/Generator/States/SphereExactVectorState.hpp"
#include "QuICC/Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "QuICC/Generator/Visualizers/VectorFieldVisualizer.hpp"
#include "QuICC/Io/Variable/FieldProbeWriter.hpp"
#include "QuICC/Io/Variable/SphereAngularMomentumWriter.hpp"
#include "QuICC/Io/Variable/SphereNusseltWriter.hpp"
#include "QuICC/Io/Variable/SphereScalarEnergyWriter.hpp"
#include "QuICC/Io/Variable/SphereScalarLSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereScalarMSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereScalarNSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolEnergyWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolLSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolMSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolNSpectrumWriter.hpp"
#include "QuICC/Io/Variable/StateFileReader.hpp"
#include "QuICC/Io/Variable/StateFileWriter.hpp"
#include "QuICC/Io/Variable/VisualizationFileWriter.hpp"
#include "QuICC/NonDimensional/CflInertial.hpp"
#include "QuICC/NonDimensional/Ekman.hpp"
#include "QuICC/NonDimensional/Prandtl.hpp"
#include "QuICC/NonDimensional/Rayleigh.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/SpectralKernels/MakeRandom.hpp"
#include "QuICC/Transform/Path/ValueScalar.hpp"
#include "QuICC/Transform/Path/ValueTorPol.hpp"
#include "QuICC/Transform/Path/NoSlipTorPol.hpp"
#include "QuICC/Transform/Path/InsulatingTorPol.hpp"
#include "QuICC/Transform/Path/StressFreeTorPol.hpp"
#include "QuICC/Bc/Name/FixedTemperature.hpp"
#include "QuICC/Bc/Name/NoSlip.hpp"
#include "QuICC/Bc/Name/StressFree.hpp"
#include "QuICC/Generator/States/Kernels/Sphere/TorPolHarmonic.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace RTC {

VectorFormulation::Id IRTCModel::SchemeFormulation()
{
   return VectorFormulation::TORPOL;
}

std::string IRTCModel::version() const
{
   return std::string(gitHash);
}

void IRTCModel::addEquations(SharedSimulation spSim)
{
   // Add transport equation
   spSim->addEquation<Equations::Boussinesq::Sphere::RTC::Transport>(
      this->spBackend());

   // Add Navier-Stokes equation
   spSim->addEquation<Equations::Boussinesq::Sphere::RTC::Momentum>(
      this->spBackend());
}

std::size_t IRTCModel::pathId(std::shared_ptr<SimulationBoundary> spBcs, const std::size_t fieldId) const
{
   std::size_t pathId;

   // Temperature
   if(fieldId == PhysicalNames::Temperature::id())
   {
      if(spBcs->bcId(fieldId) == Bc::Name::FixedTemperature::id())
      {
         pathId = Transform::Path::ValueScalar::id();
      }
      else
      {
         throw std::logic_error("Boundary condition for Temperature not implemented");
      }
   }
   // Velocity
   else if(fieldId == PhysicalNames::Velocity::id())
   {
      if(spBcs->bcId(fieldId) == Bc::Name::NoSlip::id())
      {
#if defined QUICC_BESSEL_VELOCITY_BC_VALUE_TOR_VALUE_POL
         pathId = Transform::Path::ValueTorPol::id();
#elif defined QUICC_BESSEL_VELOCITY_BC_VALUE_TOR_INSULATING_POL
         pathId = Transform::Path::InsulatingTorPol::id();
#elif defined QUICC_BESSEL_VELOCITY_BC_VALUE_TOR_NS_POL
         pathId = Transform::Path::NoSlipTorPol::id();
#else
#error "Unknown basis setup for Velocity field"
#endif
      }
      else if(spBcs->bcId(fieldId) == Bc::Name::StressFree::id())
      {
         pathId = Transform::Path::StressFreeTorPol::id();
      }
      else
      {
         throw std::logic_error("Boundary condition for Temperature not implemented");
      }
   }

   return pathId;
}

void IRTCModel::addStates(SharedStateGenerator spGen)
{
   // Create boundary object
   auto spBcs = spGen->createBoundary();
   std::size_t tempPathId = this->pathId(spBcs, PhysicalNames::Temperature::id());
   std::size_t velPathId = this->pathId(spBcs, PhysicalNames::Velocity::id());

   // Shared pointer to equation
   Equations::SharedSphereExactScalarState spScalar;
   Equations::SharedSphereExactVectorState spVector;

   Spectral::Kernel::Complex3DMapType tSH;
   std::pair<Spectral::Kernel::Complex3DMapType::iterator, bool> ptSH;

   // Add temperature initial state generator
   spScalar =
      spGen->addEquation<Equations::SphereExactScalarState>(this->spBackend());
   spScalar->setIdentity(PhysicalNames::Temperature::id());
   spScalar->setBackwardPath(tempPathId);
   spScalar->setForwardPath(tempPathId);
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
   spVector->setBackwardPath(velPathId);
   spVector->setForwardPath(velPathId);
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

   case 5: {
      auto spKernel = std::make_shared<Physical::Kernel::Sphere::TorPolHarmonic>();
      // Toroidal
      tSH.clear();
      ptSH = tSH.insert(
         std::make_pair(std::make_pair(1, 1), std::map<int, MHDComplex>()));
      ptSH.first->second.insert(std::make_pair(0, MHDComplex(1.0,1.0)));
      ptSH.first->second.insert(std::make_pair(1, MHDComplex(1.0,1.0)));
      ptSH.first->second.insert(std::make_pair(2, MHDComplex(1.0,1.0)));
      ptSH.first->second.insert(std::make_pair(3, MHDComplex(-3.0,-3.0)));
      spKernel->setModes(FieldComponents::Spectral::TOR, tSH);
      // Poloidal
      tSH.clear();
      ptSH = tSH.insert(
         std::make_pair(std::make_pair(1, 1), std::map<int, MHDComplex>()));
      ptSH.first->second.insert(std::make_pair(0, MHDComplex(1.0,1.0)));
      ptSH.first->second.insert(std::make_pair(1, MHDComplex(1.0,1.0)));
      ptSH.first->second.insert(std::make_pair(2, MHDComplex(-5.0,-5.0)));
      ptSH.first->second.insert(std::make_pair(3, MHDComplex(3.0,3.0)));
      spKernel->setModes(FieldComponents::Spectral::POL, tSH);
      spVector->setPhysicalKernel(spKernel);
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

void IRTCModel::addVisualizers(SharedVisualizationGenerator spVis)
{
   // Create boundary object
   auto spBcs = spVis->createBoundary();
   std::size_t tempPathId = this->pathId(spBcs, PhysicalNames::Temperature::id());
   std::size_t velPathId = this->pathId(spBcs, PhysicalNames::Velocity::id());

   // Shared pointer to basic field visualizer
   Equations::SharedScalarFieldVisualizer spScalar;
   Equations::SharedVectorFieldVisualizer spVector;

   // Add temperature field visualization
   spScalar =
      spVis->addEquation<Equations::ScalarFieldVisualizer>(this->spBackend());
   spScalar->setFields(true, true);
   spScalar->setIdentity(PhysicalNames::Temperature::id());
   spScalar->setBackwardPath(tempPathId);
   spScalar->setForwardPath(tempPathId);

   // Add velocity field visualization
   spVector =
      spVis->addEquation<Equations::VectorFieldVisualizer>(this->spBackend());
   spVector->setFields(true, false, true);
   spVector->setIdentity(PhysicalNames::Velocity::id());
   spVector->setBackwardPath(velPathId);
   spVector->setForwardPath(velPathId);

   // Add output file
   auto spOut = std::make_shared<Io::Variable::VisualizationFileWriter>(
      spVis->ss().tag());
   spOut->expect(PhysicalNames::Temperature::id());
   spOut->expect(PhysicalNames::Velocity::id());
   spVis->addHdf5OutputFile(spOut);
}

std::map<std::string, std::map<std::string, int>> IRTCModel::configTags() const
{
   std::map<std::string, int> onOff;
   onOff.emplace("enable", 1);

   std::map<std::string, int> options;
   options.emplace("enable", 0);
   options.emplace("numbered", 0);
   options.emplace("only_every", 1);

   std::map<std::string, std::map<std::string, int>> tags;
   // kinetic
   tags.emplace("kinetic_energy", onOff);
   tags.emplace("kinetic_l_spectrum", options);
   tags.emplace("kinetic_m_spectrum", options);
   tags.emplace("kinetic_n_spectrum", options);
   // temperature
   tags.emplace("temperature_energy", onOff);
   tags.emplace("temperature_l_spectrum", options);
   tags.emplace("temperature_m_spectrum", options);
   tags.emplace("temperature_n_spectrum", options);
   // diagnostic
   tags.emplace("angular_momentum", onOff);
   tags.emplace("nusselt", onOff);

   return tags;
}

void IRTCModel::addAsciiOutputFiles(SharedSimulation spSim)
{
   // Create boundary object
   auto spBcs = spSim->createBoundary();
   std::size_t tempPathId = this->pathId(spBcs, PhysicalNames::Temperature::id());
   std::size_t velPathId = this->pathId(spBcs, PhysicalNames::Velocity::id());

   // Create Nusselt writer
   this->enableAsciiFile<Io::Variable::SphereNusseltWriter>("nusselt", "",
      PhysicalNames::Temperature::id(), spSim);

   // Create temperature energy writer
   {
      auto spF = this->enableAsciiFile<Io::Variable::SphereScalarEnergyWriter>(
         "temperature_energy", "temperature", PhysicalNames::Temperature::id(),
         spSim);
      spF->setTransformPath(tempPathId);
   }

   // Create temperature L energy spectrum writer
   {
      auto spF = this->enableAsciiFile<Io::Variable::SphereScalarLSpectrumWriter>(
         "temperature_l_spectrum", "temperature", PhysicalNames::Temperature::id(),
         spSim);
         spF->setTransformPath(tempPathId);
   }

   // Create temperature M energy spectrum writer
   {   
      auto spF = this->enableAsciiFile<Io::Variable::SphereScalarMSpectrumWriter>(
         "temperature_m_spectrum", "temperature", PhysicalNames::Temperature::id(),
         spSim);
      spF->setTransformPath(tempPathId);
   }

   // Create temperature N power spectrum writer
   {   
      auto spF = this->enableAsciiFile<Io::Variable::SphereScalarNSpectrumWriter>(
            "temperature_n_spectrum", "temperature", PhysicalNames::Temperature::id(),
            spSim);
      spF->setTransformPath(tempPathId);
   }

   // Create kinetic energy writer
   {
      auto spF = this->enableAsciiFile<Io::Variable::SphereTorPolEnergyWriter>(
         "kinetic_energy", "kinetic", PhysicalNames::Velocity::id(), spSim);
      spF->setTransformPath(velPathId);
   }

   // Create kinetic L energy spectrum writer
   {
      auto spF = this->enableAsciiFile<Io::Variable::SphereTorPolLSpectrumWriter>(
         "kinetic_l_spectrum", "kinetic", PhysicalNames::Velocity::id(), spSim);
         spF->setTransformPath(velPathId);
   }

   // Create kinetic M energy spectrum writer
   {
      auto spF = this->enableAsciiFile<Io::Variable::SphereTorPolMSpectrumWriter>(
         "kinetic_m_spectrum", "kinetic", PhysicalNames::Velocity::id(), spSim);
         spF->setTransformPath(velPathId);
   }

   // Create kinetic N power spectrum writer
   {
      auto spF = this->enableAsciiFile<Io::Variable::SphereTorPolNSpectrumWriter>(
         "kinetic_n_spectrum", "kinetic", PhysicalNames::Velocity::id(), spSim);
         spF->setTransformPath(velPathId);
   }

   // Create angular momentum writer
   this->enableAsciiFile<Io::Variable::SphereAngularMomentumWriter>(
      "angular_momentum", "", PhysicalNames::Velocity::id(), spSim);

   // Examples of field physical space probes
   //
   const bool probeVelocity = false;
   const bool probeTemperature = false;

   // Add Velocity probe
   if (probeVelocity)
   {
      std::vector<MHDFloat> pos = {0.6, Math::PI / 2.0, 0};
      auto spFile = std::make_shared<Io::Variable::FieldProbeWriter>(
         "velocity_", spSim->ss().tag(), pos);
      spFile->expect(PhysicalNames::Velocity::id());
      spSim->addAsciiOutputFile(spFile);
   }

   // Add Temperature probe
   if (probeTemperature)
   {
      std::vector<MHDFloat> pos = {0.6, Math::PI / 2.0, 0};
      auto spFile = std::make_shared<Io::Variable::FieldProbeWriter>(
         "temperature_", spSim->ss().tag(), pos);
      spFile->expect(PhysicalNames::Temperature::id());
      spSim->addAsciiOutputFile(spFile);
   }
}

} // namespace RTC
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
