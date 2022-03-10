/**
 * @file IRTCModel.cpp
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
#include "QuICC/Model/Boussinesq/Sphere/RTC/IRTCModel.hpp"

// Project includes
//
#include "QuICC/Model/Boussinesq/Sphere/RTC/Transport.hpp"
#include "QuICC/Model/Boussinesq/Sphere/RTC/Momentum.hpp"
#include "QuICC/Enums/FieldIds.hpp"
#include "QuICC/PhysicalNames/Temperature.hpp"
#include "QuICC/PhysicalNames/Velocity.hpp"
#include "QuICC/NonDimensional/Ekman.hpp"
#include "QuICC/NonDimensional/Rayleigh.hpp"
#include "QuICC/NonDimensional/Prandtl.hpp"
#include "QuICC/NonDimensional/CflInertial.hpp"
#include "QuICC/Io/Variable/StateFileReader.hpp"
#include "QuICC/Io/Variable/StateFileWriter.hpp"
#include "QuICC/Io/Variable/VisualizationFileWriter.hpp"
#include "QuICC/Io/Variable/SphereNusseltWriter.hpp"
#include "QuICC/Io/Variable/SphereScalarEnergyWriter.hpp"
#include "QuICC/Io/Variable/SphereScalarLSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereScalarMSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereScalarNSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolEnergyWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolLSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolMSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereTorPolNSpectrumWriter.hpp"
#include "QuICC/Io/Variable/SphereAngularMomentumWriter.hpp"
#include "QuICC/Generator/States/RandomScalarState.hpp"
#include "QuICC/Generator/States/RandomVectorState.hpp"
#include "QuICC/Generator/States/SphereExactScalarState.hpp"
#include "QuICC/Generator/States/SphereExactVectorState.hpp"
#include "QuICC/Generator/States/Kernels/Sphere/BenchmarkTempC1.hpp"
#include "QuICC/Generator/States/Kernels/Sphere/BenchmarkTempC1.hpp"
#include "QuICC/Generator/Visualizers/ScalarFieldVisualizer.hpp"
#include "QuICC/Generator/Visualizers/VectorFieldVisualizer.hpp"
#include "QuICC/SpectralKernels/MakeRandom.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace RTC {

   VectorFormulation::Id IRTCModel::SchemeFormulation()
   {
      return VectorFormulation::TORPOL;
   }

   void IRTCModel::addEquations(SharedSimulation spSim)
   {
      // Add transport equation
      spSim->addEquation<Equations::Boussinesq::Sphere::RTC::Transport>(this->spBackend());

      // Add Navier-Stokes equation
      spSim->addEquation<Equations::Boussinesq::Sphere::RTC::Momentum>(this->spBackend());
   }

   void IRTCModel::addStates(SharedStateGenerator spGen)
   {
      // Shared pointer to equation
      Equations::SharedSphereExactScalarState spScalar;
      Equations::SharedSphereExactVectorState spVector;

      Spectral::Kernel::Complex3DMapType tSH;
      std::pair<Spectral::Kernel::Complex3DMapType::iterator,bool> ptSH;

      // Add temperature initial state generator
      spScalar = spGen->addEquation<Equations::SphereExactScalarState>(this->spBackend());
      spScalar->setIdentity(PhysicalNames::Temperature::id());
      switch(3)
      {
         case 0:
            {
               spScalar->setPhysicalNoise(1e-15);
            }
            break;

         case 1:
            {
               spScalar->setPhysicalConstant(1.0);
            }
            break;

         case 2:
            {
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(3,3), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0,2.0)));
               spScalar->setSpectralModes(tSH);
            }
            break;

         case 3:
            {
               auto spKernel = std::make_shared<Physical::Kernel::Sphere::BenchmarkTempC1>();
               spKernel->init(0.0, 1e-5);
               spScalar->setPhysicalKernel(spKernel);
            }
            break;

         case 4:
            {
               auto spKernel = std::make_shared<Spectral::Kernel::MakeRandom>(spGen->ss().has(SpatialScheme::Feature::ComplexSpectrum));
               std::vector<MHDFloat> ratios = {1e4, 1e4, 1e4};
               spKernel->setRatio(ratios);
               spKernel->init(-1e-4, 1e-4);
               spScalar->setSrcKernel(spKernel);
            }
            break;
      }

      // Add velocity initial state generator
      spVector = spGen->addEquation<Equations::SphereExactVectorState>(this->spBackend());
      spVector->setIdentity(PhysicalNames::Velocity::id());
      switch(3)
      {
         // Toroidal only
         case 0:
            {
               // Toroidal
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(1,1), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setSpectralModes(FieldComponents::Spectral::TOR, tSH);
               // Poloidal
               tSH.clear();
               spVector->setSpectralModes(FieldComponents::Spectral::POL, tSH);
            }
            break;

         // Poloidal only
         case 1:
            {
               // Toroidal
               tSH.clear();
               spVector->setSpectralModes(FieldComponents::Spectral::TOR, tSH);
               // Poloidal
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(2,0), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setSpectralModes(FieldComponents::Spectral::POL, tSH);
            }
            break;

         // Toroidal & Poloidal
         case 2:
            {
               // Toroidal
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(1,1), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setSpectralModes(FieldComponents::Spectral::TOR, tSH);
               // Poloidal
               tSH.clear();
               ptSH = tSH.insert(std::make_pair(std::make_pair(2,0), std::map<int,MHDComplex>()));
               ptSH.first->second.insert(std::make_pair(7, MHDComplex(1.0)));
               spVector->setSpectralModes(FieldComponents::Spectral::POL, tSH);
            }
            break;

         case 3:
            {
               auto spKernel = std::make_shared<Spectral::Kernel::MakeRandom>(spGen->ss().has(SpatialScheme::Feature::ComplexSpectrum));
               std::vector<MHDFloat> ratios = {1e2, 1e2, 1e2};
               spKernel->setRatio(ratios);
               spKernel->init(-1e-15, 1e-15);
               spVector->setSrcKernel(FieldComponents::Spectral::TOR, spKernel);
               spVector->setSrcKernel(FieldComponents::Spectral::POL, spKernel);
            }
            break;

         case 4:
            {
               auto spKernel = std::make_shared<Spectral::Kernel::MakeRandom>(spGen->ss().has(SpatialScheme::Feature::ComplexSpectrum));
               std::vector<MHDFloat> ratios = {1e4, 1e4, 1e4};
               spKernel->setRatio(ratios);
               spKernel->init(-1e-4, 1e-4);
               spVector->setSrcKernel(FieldComponents::Spectral::TOR, spKernel);
               spVector->setSrcKernel(FieldComponents::Spectral::POL, spKernel);
            }
            break;
      }

      // Add output file
      auto spOut = std::make_shared<Io::Variable::StateFileWriter>(spGen->ss().tag(), spGen->ss().has(SpatialScheme::Feature::RegularSpectrum));
      spOut->expect(PhysicalNames::Temperature::id());
      spOut->expect(PhysicalNames::Velocity::id());
      spGen->addHdf5OutputFile(spOut);
   }

   void IRTCModel::addVisualizers(SharedVisualizationGenerator spVis)
   {
      // Shared pointer to basic field visualizer
      Equations::SharedScalarFieldVisualizer spScalar;
      Equations::SharedVectorFieldVisualizer spVector;

      // Add temperature field visualization
      spScalar = spVis->addEquation<Equations::ScalarFieldVisualizer>(this->spBackend());
      spScalar->setFields(true, true);
      spScalar->setIdentity(PhysicalNames::Temperature::id());

      // Add velocity field visualization
      spVector = spVis->addEquation<Equations::VectorFieldVisualizer>(this->spBackend());
      spVector->setFields(true, false, true);
      spVector->setIdentity(PhysicalNames::Velocity::id());

      // Add output file
      auto spOut = std::make_shared<Io::Variable::VisualizationFileWriter>(spVis->ss().tag());
      spOut->expect(PhysicalNames::Temperature::id());
      spOut->expect(PhysicalNames::Velocity::id());
      spVis->addHdf5OutputFile(spOut);
   }

   std::map<std::string, std::map<std::string,int> > IRTCModel::configTags() const
   {
      std::map<std::string,int> onOff;
      onOff.emplace("enable", 1);

      std::map<std::string,int> options;
      options.emplace("enable", 0);
      options.emplace("numbered", 0);
      options.emplace("only_every", 1);

      std::map<std::string,std::map<std::string,int> > tags;
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
      std::string tag;

      // Create Nusselt writer
      tag = "nusselt";
      if(spSim->config().model(tag).at("enable"))
      {
         auto spFile = std::make_shared<Io::Variable::SphereNusseltWriter>("", spSim->ss().tag());
         spFile->expect(PhysicalNames::Temperature::id());
         spSim->addAsciiOutputFile(spFile);
      }

      // Create temperature energy writer
      tag = "temperature_energy";
      if(spSim->config().model(tag).at("enable"))
      {
         auto spFile = std::make_shared<Io::Variable::SphereScalarEnergyWriter>("temperature", spSim->ss().tag());
         spFile->expect(PhysicalNames::Temperature::id());
         spSim->addAsciiOutputFile(spFile);
      }

      // Create temperature L energy spectrum writer
      tag = "temperature_l_spectrum";
      if(spSim->config().model(tag).at("enable"))
      {
         auto spFile = std::make_shared<Io::Variable::SphereScalarLSpectrumWriter>("temperature", spSim->ss().tag());
         spFile->expect(PhysicalNames::Temperature::id());
         if(spSim->config().model(tag).at("numbered"))
         {
            spFile->numberOutput();
         }
         spFile->onlyEvery(spSim->config().model(tag).at("only_every"));
         spSim->addAsciiOutputFile(spFile);
      }

      // Create temperature M energy spectrum writer
      tag = "temperature_m_spectrum";
      if(spSim->config().model(tag).at("enable"))
      {
         auto spspFile = std::make_shared<Io::Variable::SphereScalarMSpectrumWriter>("temperature", spSim->ss().tag());
         spspFile->expect(PhysicalNames::Temperature::id());
         if(spSim->config().model(tag).at("numbered"))
         {
            spspFile->numberOutput();
         }
         spspFile->onlyEvery(spSim->config().model(tag).at("only_every"));
         spSim->addAsciiOutputFile(spspFile);
      }

      // Create temperature N power spectrum writer
      tag = "temperature_n_spectrum";
      if(spSim->config().model(tag).at("enable"))
      {
         auto spFile = std::make_shared<Io::Variable::SphereScalarNSpectrumWriter>("temperature", spSim->ss().tag());
         spFile->expect(PhysicalNames::Temperature::id());
         if(spSim->config().model(tag).at("numbered"))
         {
            spFile->numberOutput();
         }
         spFile->onlyEvery(spSim->config().model(tag).at("only_every"));
         spSim->addAsciiOutputFile(spFile);
      }

      // Create kinetic energy writer
      tag = "kinetic_energy";
      if(spSim->config().model(tag).at("enable"))
      {
         auto spFile = std::make_shared<Io::Variable::SphereTorPolEnergyWriter>("kinetic", spSim->ss().tag());
         spFile->expect(PhysicalNames::Velocity::id());
         spSim->addAsciiOutputFile(spFile);
      }

      // Create kinetic L energy spectrum writer
      tag = "kinetic_l_spectrum";
      if(spSim->config().model(tag).at("enable"))
      {
         auto spFile = std::make_shared<Io::Variable::SphereTorPolLSpectrumWriter>("kinetic", spSim->ss().tag());
         spFile->expect(PhysicalNames::Velocity::id());
         if(spSim->config().model(tag).at("numbered"))
         {
            spFile->numberOutput();
         }
         spFile->onlyEvery(spSim->config().model(tag).at("only_every"));
         spSim->addAsciiOutputFile(spFile);
      }

      // Create kinetic M energy spectrum writer
      tag = "kinetic_m_spectrum";
      if(spSim->config().model(tag).at("enable"))
      {
         auto spFile = std::make_shared<Io::Variable::SphereTorPolMSpectrumWriter>("kinetic", spSim->ss().tag());
         spFile->expect(PhysicalNames::Velocity::id());
         if(spSim->config().model(tag).at("numbered"))
         {
            spFile->numberOutput();
         }
         spFile->onlyEvery(spSim->config().model(tag).at("only_every"));
         spSim->addAsciiOutputFile(spFile);
      }

      // Create kinetic N power spectrum writer
      tag = "kinetic_n_spectrum";
      if(spSim->config().model(tag).at("enable"))
      {
         auto spFile = std::make_shared<Io::Variable::SphereTorPolNSpectrumWriter>("kinetic", spSim->ss().tag());
         spFile->expect(PhysicalNames::Velocity::id());
         if(spSim->config().model(tag).at("numbered"))
         {
            spFile->numberOutput();
         }
         spFile->onlyEvery(spSim->config().model(tag).at("only_every"));
         spSim->addAsciiOutputFile(spFile);
      }

      // Create angular momentum writer
      tag = "angular_momentum";
      if(spSim->config().model(tag).at("enable"))
      {
         auto spFile = std::make_shared<Io::Variable::SphereAngularMomentumWriter>("", spSim->ss().tag());
         spFile->expect(PhysicalNames::Velocity::id());
         spSim->addAsciiOutputFile(spFile);
      }
   }

}
}
}
}
}
