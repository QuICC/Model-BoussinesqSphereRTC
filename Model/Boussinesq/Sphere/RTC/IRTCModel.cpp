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
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"

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
   // Create Nusselt writer
   this->enableAsciiFile<Io::Variable::SphereNusseltWriter>("nusselt", "",
      PhysicalNames::Temperature::id(), spSim);

   // Create temperature energy writer
   this->enableAsciiFile<Io::Variable::SphereScalarEnergyWriter>(
      "temperature_energy", "temperature", PhysicalNames::Temperature::id(),
      spSim);

   // Create temperature L energy spectrum writer
   this->enableAsciiFile<Io::Variable::SphereScalarLSpectrumWriter>(
      "temperature_l_spectrum", "temperature", PhysicalNames::Temperature::id(),
      spSim);

   // Create temperature M energy spectrum writer
   this->enableAsciiFile<Io::Variable::SphereScalarMSpectrumWriter>(
      "temperature_m_spectrum", "temperature", PhysicalNames::Temperature::id(),
      spSim);

   // Create temperature N power spectrum writer
   this->enableAsciiFile<Io::Variable::SphereScalarNSpectrumWriter>(
      "temperature_n_spectrum", "temperature", PhysicalNames::Temperature::id(),
      spSim);

   // Create kinetic energy writer
   this->enableAsciiFile<Io::Variable::SphereTorPolEnergyWriter>(
      "kinetic_energy", "kinetic", PhysicalNames::Velocity::id(), spSim);

   // Create kinetic L energy spectrum writer
   this->enableAsciiFile<Io::Variable::SphereTorPolLSpectrumWriter>(
      "kinetic_l_spectrum", "kinetic", PhysicalNames::Velocity::id(), spSim);

   // Create kinetic M energy spectrum writer
   this->enableAsciiFile<Io::Variable::SphereTorPolMSpectrumWriter>(
      "kinetic_m_spectrum", "kinetic", PhysicalNames::Velocity::id(), spSim);

   // Create kinetic N power spectrum writer
   this->enableAsciiFile<Io::Variable::SphereTorPolNSpectrumWriter>(
      "kinetic_n_spectrum", "kinetic", PhysicalNames::Velocity::id(), spSim);

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
