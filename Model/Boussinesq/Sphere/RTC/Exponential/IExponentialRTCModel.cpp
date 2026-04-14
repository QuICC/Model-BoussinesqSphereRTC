/**
 * @file IExponentialRTCModel.cpp
 * @brief Source of the Boussinesq rotating thermal convection in a sphere
 * (Toroidal/Poloidal formulation)
 */

// System includes
//

// Project includes
//
#include "Model/Boussinesq/Sphere/RTC/Exponential/IExponentialRTCModel.hpp"
#include "Model/Boussinesq/Sphere/RTC/Momentum.hpp"
#include "Model/Boussinesq/Sphere/RTC/Transport.hpp"
#include "Model/Boussinesq/Sphere/RTC/Exponential/MomentumJacobian.hpp"
#include "Model/Boussinesq/Sphere/RTC/Exponential/TransportJacobian.hpp"
#include "Model/Boussinesq/Sphere/RTC/gitHash.hpp"
#include "QuICC/PhysicalNames/JacobianTemperature.hpp"
#include "QuICC/PhysicalNames/JacobianVelocity.hpp"

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Sphere {

namespace RTC {

namespace Exponential {

std::vector<std::size_t> IExponentialRTCModel::excludedFieldIds() const
{
   std::vector<std::size_t> fields = {
      PhysicalNames::JacobianVelocity::id(),
      PhysicalNames::JacobianTemperature::id(),
   };

   return fields;
}

void IExponentialRTCModel::addEquations(SharedSimulation spSim)
{
   auto optZero = std::make_shared<Equations::EquationOptions>(0, false);

   // Add transport equation
   spSim->addEquation<Equations::Boussinesq::Sphere::RTC::Transport>(
      this->spBackend(), optZero);

   // Add Navier-Stokes equation
   spSim->addEquation<Equations::Boussinesq::Sphere::RTC::Momentum>(
      this->spBackend(), optZero);

#if 1
   auto optOne = std::make_shared<Equations::EquationOptions>(1, false);

   // Add transport jacobian equation
   spSim->addEquation<Equations::Boussinesq::Sphere::RTC::Exponential::TransportJacobian>(
      this->spBackend(), optOne);

   // Add Navier-Stokes jacobian equation
   spSim->addEquation<Equations::Boussinesq::Sphere::RTC::Exponential::MomentumJacobian>(
      this->spBackend(), optOne);
#endif

   #ifdef QUICC_USE_MLIR_GRAPH
   // Add Graph

   // At this point we only support the explicit model
   // If the jitter is not set up, we default to the old implementation
   if (spSim->ss().has(SpatialScheme::Feature::SpectralOrdering123))
   {
      return;
   }

   std::string graphStr = R"mlir(
// type aliases
!real = tensor<?x?x?xf64>
!complex = tensor<?x?x?xcomplex<f64>>

func.func private @bwdScalar(%S: !complex) -> !real {
    // backward scalar path
    %S1 = quiccir.jw.prj %S : !complex -> !complex attributes{kind = "P"}
    %S1T = quiccir.transpose %S1 permutation = [1, 2, 0] : !complex -> !complex
    %S2 = quiccir.al.prj %S1T : !complex -> !complex attributes{kind = "P"}
    %S2T = quiccir.transpose %S2 permutation = [1, 2, 0] : !complex -> !complex
    %S3 = quiccir.fr.prj %S2T : !complex -> !real attributes{kind = "P"}
    return %S3 : !real
}

func.func private @bwdGradScalar(%S: !complex) -> (!real, !real, !real) {
    // grad R
    %SdR1 = quiccir.jw.prj %S : !complex -> !complex attributes{kind = "D1"}
    %SdR1T = quiccir.transpose %SdR1 permutation = [1, 2, 0] : !complex -> !complex
    %SdR2 = quiccir.al.prj %SdR1T : !complex -> !complex attributes{kind = "P"}
    %SdR2T = quiccir.transpose %SdR2 permutation = [1, 2, 0] : !complex -> !complex
    %SdR = quiccir.fr.prj %SdR2T : !complex -> !real attributes{kind = "P"}
    // grad Theta
    %SdTh1 = quiccir.jw.prj %S : !complex -> !complex attributes{kind = "DivR1_Zero"}
    %SdTh1T = quiccir.transpose %SdTh1 permutation = [1, 2, 0] : !complex -> !complex
    %SdTh2 = quiccir.al.prj %SdTh1T : !complex -> !complex attributes{kind = "D1"}
    %SdTh2T = quiccir.transpose %SdTh2 permutation = [1, 2, 0] : !complex -> !complex
    %SdTh = quiccir.fr.prj %SdTh2T : !complex -> !real attributes{kind = "P"}
    // grad Phi
    %SdPh1 = quiccir.jw.prj %S : !complex -> !complex attributes{kind = "DivR1_Zero"}
    %SdPh1T = quiccir.transpose %SdPh1 permutation = [1, 2, 0] : !complex -> !complex
    %SdPh2 = quiccir.al.prj %SdPh1T : !complex -> !complex attributes{kind = "DivS1Dp"}
    %SdPh2T = quiccir.transpose %SdPh2 permutation = [1, 2, 0] : !complex -> !complex
    %SdPh = quiccir.fr.prj %SdPh2T : !complex -> !real attributes{kind = "P"}
    return %SdR, %SdTh, %SdPh : !real, !real, !real
}

func.func private @bwdVector(%Tor: !complex, %Pol: !complex) -> (!real, !real, !real) {
    // R
    %PolR1 = quiccir.jw.prj %Pol : !complex -> !complex attributes{kind = "DivR1_Zero"}
    %PolR1T = quiccir.transpose %PolR1 permutation = [1, 2, 0] : !complex -> !complex
    %PolR2 = quiccir.al.prj %PolR1T : !complex -> !complex attributes{kind = "Ll"}
    %PolR2T = quiccir.transpose %PolR2 permutation = [1, 2, 0] : !complex -> !complex
    %R = quiccir.fr.prj %PolR2T : !complex -> !real attributes{kind = "P"}
    // Theta
    %TorTh1 = quiccir.jw.prj %Tor : !complex -> !complex attributes{kind = "P"}
    %TorTh1T = quiccir.transpose %TorTh1 permutation = [1, 2, 0] : !complex -> !complex
    %TorTh2 = quiccir.al.prj %TorTh1T : !complex -> !complex attributes{kind = "DivS1Dp"}
    %TorTh2T = quiccir.transpose %TorTh2 permutation = [1, 2, 0] : !complex -> !complex
    %TorTh3 = quiccir.fr.prj %TorTh2T : !complex -> !real attributes{kind = "P"}
    //
    %PolTh1 = quiccir.jw.prj %Pol : !complex -> !complex attributes{kind = "DivR1D1R1_Zero"}
    %PolTh1T = quiccir.transpose %PolTh1 permutation = [1, 2, 0] : !complex -> !complex
    %PolTh2 = quiccir.al.prj %PolTh1T : !complex -> !complex attributes{kind = "D1"}
    %PolTh2T = quiccir.transpose %PolTh2 permutation = [1, 2, 0] : !complex -> !complex
    %PolTh3 = quiccir.fr.prj %PolTh2T : !complex -> !real attributes{kind = "P"}
    //
    %Theta = quiccir.add %TorTh3, %PolTh3 : !real, !real -> !real
    // Phi
    %TorPh1 = quiccir.jw.prj %Tor : !complex -> !complex attributes{kind = "P"}
    %TorPh1T = quiccir.transpose %TorPh1 permutation = [1, 2, 0] : !complex -> !complex
    %TorPh2 = quiccir.al.prj %TorPh1T : !complex -> !complex attributes{kind = "D1"}
    %TorPh2T = quiccir.transpose %TorPh2 permutation = [1, 2, 0] : !complex -> !complex
    %TorPh3 = quiccir.fr.prj %TorPh2T : !complex -> !real attributes{kind = "P"}
    //
    %PolPh1 = quiccir.jw.prj %Pol : !complex -> !complex attributes{kind = "DivR1D1R1_Zero"}
    %PolPh1T = quiccir.transpose %PolPh1 permutation = [1, 2, 0] : !complex -> !complex
    %PolPh2 = quiccir.al.prj %PolPh1T : !complex -> !complex attributes{kind = "DivS1Dp"}
    %PolPh2T = quiccir.transpose %PolPh2 permutation = [1, 2, 0] : !complex -> !complex
    %PolPh3 = quiccir.fr.prj %PolPh2T : !complex -> !real attributes{kind = "P"}
    //
    %Phi = quiccir.sub %PolPh3, %TorPh3 : !real, !real -> !real
    return %R, %Theta, %Phi : !real, !real, !real
}

func.func private @bwdCurl(%Tor: !complex, %Pol: !complex) -> (!real, !real, !real) {
    // R
    %TorR1 = quiccir.jw.prj %Tor : !complex -> !complex attributes{kind = "DivR1_Zero"}
    %TorR1T = quiccir.transpose %TorR1 permutation = [1, 2, 0] : !complex -> !complex
    %TorR2 = quiccir.al.prj %TorR1T : !complex -> !complex attributes{kind = "Ll"}
    %TorR2T = quiccir.transpose %TorR2 permutation = [1, 2, 0] : !complex -> !complex
    %R = quiccir.fr.prj %TorR2T : !complex -> !real attributes{kind = "P"}
    // Theta
    %TorTh1 = quiccir.jw.prj %Tor : !complex -> !complex attributes{kind = "DivR1D1R1_Zero"}
    %TorTh1T = quiccir.transpose %TorTh1 permutation = [1, 2, 0] : !complex -> !complex
    %TorTh2 = quiccir.al.prj %TorTh1T : !complex -> !complex attributes{kind = "D1"}
    %TorTh2T = quiccir.transpose %TorTh2 permutation = [1, 2, 0] : !complex -> !complex
    %TorTh3 = quiccir.fr.prj %TorTh2T : !complex -> !real attributes{kind = "P"}
    //
    %PolTh1 = quiccir.jw.prj %Pol : !complex -> !complex attributes{kind = "SphLapl"}
    %PolTh1T = quiccir.transpose %PolTh1 permutation = [1, 2, 0] : !complex -> !complex
    %PolTh2 = quiccir.al.prj %PolTh1T : !complex -> !complex attributes{kind = "DivS1Dp"}
    %PolTh2T = quiccir.transpose %PolTh2 permutation = [1, 2, 0] : !complex -> !complex
    %PolTh3 = quiccir.fr.prj %PolTh2T : !complex -> !real attributes{kind = "P"}
    //
    %Theta = quiccir.sub %TorTh3, %PolTh3 : !real, !real -> !real
    // Phi
    %TorPh1 = quiccir.jw.prj %Tor : !complex -> !complex attributes{kind = "DivR1D1R1_Zero"}
    %TorPh1T = quiccir.transpose %TorPh1 permutation = [1, 2, 0] : !complex -> !complex
    %TorPh2 = quiccir.al.prj %TorPh1T : !complex -> !complex attributes{kind = "DivS1Dp"}
    %TorPh2T = quiccir.transpose %TorPh2 permutation = [1, 2, 0] : !complex -> !complex
    %TorPh3 = quiccir.fr.prj %TorPh2T : !complex -> !real attributes{kind = "P"}
    //
    %PolPh1 = quiccir.jw.prj %Pol : !complex -> !complex attributes{kind = "SphLapl"}
    %PolPh1T = quiccir.transpose %PolPh1 permutation = [1, 2, 0] : !complex -> !complex
    %PolPh2 = quiccir.al.prj %PolPh1T : !complex -> !complex attributes{kind = "D1"}
    %PolPh2T = quiccir.transpose %PolPh2 permutation = [1, 2, 0] : !complex -> !complex
    %PolPh3 = quiccir.fr.prj %PolPh2T : !complex -> !real attributes{kind = "P"}
    //
    %Phi = quiccir.add %PolPh3, %TorPh3 : !real, !real -> !real
    return %R, %Theta, %Phi : !real, !real, !real
}

func.func private @fwdScalar(%S: !real) -> !complex {
    // forward scalar Nl path
    %S1 = quiccir.fr.int %S : !real -> !complex attributes{kind = "P"}
    %S1T = quiccir.transpose %S1 permutation = [2, 0, 1] : !complex -> !complex
    %S2 = quiccir.al.int %S1T : !complex -> !complex attributes{kind = "P"}
    %S2T = quiccir.transpose %S2 permutation = [2, 0, 1] : !complex -> !complex
    %S3 = quiccir.jw.int %S2T : !complex -> !complex attributes{kind = "P"}

    return %S3 : !complex
}

func.func private @fwdVector(%R: !real, %Theta: !real, %Phi: !real) -> (!complex, !complex){
    // Tor
    %ThetaTor1 = quiccir.fr.int %Theta : !real -> !complex attributes{kind = "P"}
    %ThetaTor1T = quiccir.transpose %ThetaTor1 permutation = [2, 0, 1] : !complex -> !complex
    %ThetaTor2 = quiccir.al.int %ThetaTor1T : !complex -> !complex attributes{kind = "DivLlDivS1Dp"}
    %ThetaTor2T = quiccir.transpose %ThetaTor2 permutation = [2, 0, 1] : !complex -> !complex
    %ThetaTor3 = quiccir.jw.int %ThetaTor2T : !complex -> !complex attributes{kind = "P_Zero"}
    //
    %PhiTor1 = quiccir.fr.int %Phi : !real -> !complex attributes{kind = "P"}
    %PhiTor1T = quiccir.transpose %PhiTor1 permutation = [2, 0, 1] : !complex -> !complex
    %PhiTor2 = quiccir.al.int %PhiTor1T : !complex -> !complex attributes{kind = "DivLlD1"}
    %PhiTor2T = quiccir.transpose %PhiTor2 permutation = [2, 0, 1] : !complex -> !complex
    %PhiTor3 = quiccir.jw.int %PhiTor2T : !complex -> !complex attributes{kind = "P_Zero"}
    //
    %Tor = quiccir.sub %ThetaTor3, %PhiTor3 : !complex, !complex -> !complex

    // Pol
    %RPol1 = quiccir.fr.int %R : !real -> !complex attributes{kind = "P"}
    %RPol1T = quiccir.transpose %RPol1 permutation = [2, 0, 1] : !complex -> !complex
    %RPol2 = quiccir.al.int %RPol1T : !complex -> !complex attributes{kind = "P"}
    %RPol2T = quiccir.transpose %RPol2 permutation = [2, 0, 1] : !complex -> !complex
    %RPol3 = quiccir.jw.int %RPol2T : !complex -> !complex attributes{kind = "DivR1_Zero"}
    //
    %ThetaPol1 = quiccir.fr.int %Theta : !real -> !complex attributes{kind = "P"}
    %ThetaPol1T = quiccir.transpose %ThetaPol1 permutation = [2, 0, 1] : !complex -> !complex
    %ThetaPol2 = quiccir.al.int %ThetaPol1T : !complex -> !complex attributes{kind = "DivLlD1"}
    %ThetaPol2T = quiccir.transpose %ThetaPol2 permutation = [2, 0, 1] : !complex -> !complex
    %ThetaPol3 = quiccir.jw.int %ThetaPol2T : !complex -> !complex attributes{kind = "DivR1D1R1_Zero"}
    //
    %PhiPol1 = quiccir.fr.int %Phi : !real -> !complex attributes{kind = "P"}
    %PhiPol1T = quiccir.transpose %PhiPol1 permutation = [2, 0, 1] : !complex -> !complex
    %PhiPol2 = quiccir.al.int %PhiPol1T : !complex -> !complex attributes{kind = "DivLlDivS1Dp"}
    %PhiPol2T = quiccir.transpose %PhiPol2 permutation = [2, 0, 1] : !complex -> !complex
    %PhiPol3 = quiccir.jw.int %PhiPol2T : !complex -> !complex attributes{kind = "DivR1D1R1_Zero"}
    //
    %tmp = quiccir.add %ThetaPol3, %PhiPol3 : !complex, !complex -> !complex
    %Pol = quiccir.sub %tmp, %RPol3 : !complex, !complex -> !complex
    return %Tor, %Pol : !complex, !complex
}

func.func private @nlScalar(%UR: !real, %UTheta: !real, %UPhi: !real,
    %TdR: !real, %TdTheta: !real, %TdPhi: !real) -> !real {
    // U dot grad T
    %DotT = quiccir.dot(%UR, %UTheta, %UPhi), (%TdR, %TdTheta, %TdPhi) :
        (!real, !real, !real), (!real, !real, !real) ->
        !real
        attributes{kind = "transport"}
    // U dot R
    %DotR = quiccir.mul.const %UR : !real -> !real attributes{kind = "transport"}
    %TPhysNl = quiccir.sub %DotT, %DotR : !real, !real -> !real
    return %TPhysNl : !real
}

func.func private @nlVector(%UR: !real, %UTheta: !real, %UPhi: !real,
    %CurlR: !real, %CurlTheta: !real, %CurlPhi: !real, %T: !real) -> (!real, !real, !real) {
    // Cross
    %Cross:3 = quiccir.cross(%CurlR, %CurlTheta, %CurlPhi), (%UR, %UTheta, %UPhi) :
        (!real, !real, !real), (!real, !real, !real) ->
        (!real, !real, !real)
        attributes{kind = "inertia"}
    // Add buoyancy
    %Buoy = quiccir.mul.const %T : !real -> !real attributes{kind = "buoyancy"}
    %RTmp = quiccir.sub %Cross#0, %Buoy : !real, !real -> !real
    // Add coriolis
    %Cor0 = quiccir.mul.const %UPhi : !real -> !real attributes{kind = "coriolis_sin"}
    %RNl = quiccir.sub %RTmp, %Cor0 : !real, !real -> !real
    %Cor1 = quiccir.mul.const %UPhi : !real -> !real attributes{kind = "coriolis_cos"}
    %ThetaNl = quiccir.sub %Cross#1, %Cor1 : !real, !real -> !real
    %Cor2a = quiccir.mul.const %UR : !real -> !real attributes{kind = "coriolis_sin"}
    %PhiTmp = quiccir.add %Cross#2, %Cor2a : !real, !real -> !real
    %Cor2b = quiccir.mul.const %UTheta : !real -> !real attributes{kind = "coriolis_cos"}
    %PhiNl = quiccir.add %PhiTmp, %Cor2b : !real, !real -> !real
    return %RNl, %ThetaNl, %PhiNl : !real, !real, !real
}


func.func @entry(%T: !complex, %Tor: !complex, %Pol: !complex) -> (!complex, !complex, !complex, !real, !real, !real) {
    %TPhys = call @bwdScalar(%T) : (!complex) -> !real
    %TGrad:3 = call @bwdGradScalar(%T) : (!complex) -> (!real, !real, !real)
    %Vel:3 = call @bwdVector(%Tor, %Pol) : (!complex, !complex) -> (!real, !real, !real)
    %Curl:3 = call @bwdCurl(%Tor, %Pol) : (!complex, !complex) -> (!real, !real, !real)
    %TPhysNl = call @nlScalar(%Vel#0, %Vel#1, %Vel#2, %TGrad#0, %TGrad#1, %TGrad#2) : (!real, !real, !real, !real, !real, !real) -> !real
    %VelNl:3 = call @nlVector(%Vel#0, %Vel#1, %Vel#2, %Curl#0, %Curl#1, %Curl#2, %TPhys) : (!real, !real, !real, !real, !real, !real, !real) -> (!real, !real, !real)
    %TNl = call @fwdScalar(%TPhysNl) : (!real) -> !complex
    %TorNl, %PolNl = call @fwdVector(%VelNl#0, %VelNl#1, %VelNl#2) : (!real, !real, !real) -> (!complex, !complex)
    return %TNl, %TorNl, %PolNl, %Vel#0, %Vel#1, %Vel#2: !complex, !complex, !complex, !real, !real, !real
}
   )mlir";

   MHDFloat T = 1.0 / spSim->eqParams()->nd(NonDimensional::Ekman::id());
   MHDFloat Ra = spSim->eqParams()->nd(NonDimensional::Rayleigh::id());
   Graph::PhysicalParameters<MHDFloat> physParams;
   physParams.transport = 1.0;
   physParams.inertia = 1.0;
   physParams.coriolis = T;
   physParams.buoyancy = Ra * T;
   spSim->addGraph(graphStr, physParams);
   #endif
}

} // namespace Exponential
} // namespace RTC
} // namespace Sphere
} // namespace Boussinesq
} // namespace Model
} // namespace QuICC
