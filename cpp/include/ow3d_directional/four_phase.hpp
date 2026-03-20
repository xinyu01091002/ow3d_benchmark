#pragma once

#include "ow3d_directional/config.hpp"
#include "ow3d_directional/mf12_runner.hpp"

#include <filesystem>
#include <vector>

namespace ow3d_directional {

struct FourPhaseSeparation {
  mf12_cpp::Matrix eta1;
  mf12_cpp::Matrix eta2;
  mf12_cpp::Matrix eta3;
  mf12_cpp::Matrix eta4;
  mf12_cpp::Matrix phi1;
  mf12_cpp::Matrix phi2;
  mf12_cpp::Matrix phi3;
  mf12_cpp::Matrix phi4;
};

bool supports_standard_four_phase_separation(const std::vector<double>& phases_deg);

FourPhaseSeparation separate_four_phase_fields(
    const std::vector<double>& phases_deg,
    const std::vector<mf12_cpp::Matrix>& eta_phases,
    const std::vector<mf12_cpp::Matrix>& phi_phases);

}  // namespace ow3d_directional
