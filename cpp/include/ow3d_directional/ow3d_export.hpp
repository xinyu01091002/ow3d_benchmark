#pragma once

#include "ow3d_directional/config.hpp"
#include "ow3d_directional/mf12_runner.hpp"

#include <filesystem>

namespace ow3d_directional {

std::filesystem::path make_case_directory_name(
    const GeneratorConfig& config,
    const CaseParameters& params);
std::filesystem::path make_case_group_directory_name(
    const GeneratorConfig& config,
    const CaseParameters& params);

void ensure_directory(const std::filesystem::path& path);
void write_batch_readme(
    const std::filesystem::path& file_name,
    const GeneratorConfig& config,
    double lx,
    double ly,
    double dx,
    double dy);

void write_ow3d_init(
    const std::filesystem::path& file_name,
    const mf12_cpp::Matrix& eta,
    const mf12_cpp::Matrix& phi,
    const CaseParameters& params);

void write_ow3d_inp(
    const std::filesystem::path& file_name,
    const GeneratorConfig& config,
    const CaseParameters& params);

void write_case_readme(
    const std::filesystem::path& file_name,
    const GeneratorConfig& config,
    const CaseParameters& params);

void write_field_csv(
    const std::filesystem::path& file_name,
    const mf12_cpp::Matrix& matrix);

}  // namespace ow3d_directional
