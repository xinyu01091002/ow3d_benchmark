#pragma once

#include "mf12_cpp/compare.hpp"
#include "mf12_cpp/spectral.hpp"
#include "ow3d_directional/config.hpp"

namespace ow3d_directional {

struct Mf12RunOutputs {
  mf12_cpp::ResultBundle nonlinear;
  mf12_cpp::ResultBundle linear;
};

Mf12RunOutputs run_directional_case(const GeneratorConfig& config, const CaseParameters& params);
mf12_cpp::LoadedCase build_loaded_case(
    int order,
    const GeneratorConfig& config,
    const CaseParameters& params);

}  // namespace ow3d_directional
