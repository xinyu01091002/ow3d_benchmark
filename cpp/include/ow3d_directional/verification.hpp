#pragma once

#include "mf12_cpp/case_io.hpp"
#include "ow3d_directional/config.hpp"
#include "ow3d_directional/mf12_runner.hpp"

#include <map>

namespace ow3d_directional {

MatlabReference load_matlab_reference(const std::filesystem::path& reference_dir);
std::map<std::string, double> compare_to_reference(
    const mf12_cpp::Matrix& eta,
    const mf12_cpp::Matrix& phi,
    const MatlabReference& reference);
bool metrics_within_default_tolerance(const std::map<std::string, double>& metrics);

}  // namespace ow3d_directional
