#pragma once

#include "ow3d_directional/config.hpp"
#include "ow3d_directional/mf12_runner.hpp"

#include <filesystem>

namespace ow3d_directional {

void write_visualizations(
    const std::filesystem::path& output_dir,
    const GeneratorConfig& config,
    const CaseParameters& params,
    const Mf12RunOutputs& outputs);

}  // namespace ow3d_directional
