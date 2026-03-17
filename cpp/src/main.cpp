#include "ow3d_directional/config.hpp"
#include "ow3d_directional/four_phase.hpp"
#include "ow3d_directional/mf12_runner.hpp"
#include "ow3d_directional/ow3d_export.hpp"
#include "ow3d_directional/spectrum.hpp"
#include "ow3d_directional/verification.hpp"
#include "ow3d_directional/visualization.hpp"

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using ow3d_directional::CaseParameters;
using ow3d_directional::GeneratorConfig;

std::size_t pair_count_for(std::size_t n_comp) {
  return n_comp < 2U ? 0U : n_comp * (n_comp - 1U) / 2U;
}

std::size_t triplet_count_for(std::size_t n_comp) {
  return n_comp < 3U ? 0U : n_comp * (n_comp - 1U) * (n_comp - 2U) / 6U;
}

void print_spectrum_summary(const CaseParameters& params) {
  const auto& spec = params.spectrum_definition;
  std::cout << "Spectrum retention: "
            << spec.candidate_components << " candidate bins -> "
            << spec.threshold_filtered_components << " after tiny-weight filter -> "
            << spec.energy_target_components << " to reach requested energy -> "
            << spec.n_components << " kept for reconstruction";
  if (spec.n_components < spec.energy_target_components) {
    std::cout << " (capped by max_components)";
  }
  std::cout << "\n";
  std::cout << "Retained weight fraction: " << std::setprecision(6) << spec.retained_weight_frac
            << " | MF12 reconstruction counts: n=" << spec.n_components
            << ", pairs=" << pair_count_for(static_cast<std::size_t>(spec.n_components))
            << ", triplets=" << triplet_count_for(static_cast<std::size_t>(spec.n_components))
            << "\n";
}

void print_usage() {
  std::cout
      << "ow3d_directional_generator commands:\n"
      << "  inspect-config <config.json>\n"
      << "  generate <config.json>\n"
      << "  verify-matlab <config.json> <reference_dir>\n";
}

CaseParameters first_case_from_config(const GeneratorConfig& config) {
  return ow3d_directional::build_case_parameters(
      config,
      config.physics.kd_list.at(0),
      config.spectrum.akp_list.at(0),
      config.spectrum.alpha_list.at(0),
      config.spectrum.phases_deg.at(0),
      config.timing.t_init_periods_list.at(0));
}

int run_inspect_config(const std::filesystem::path& config_path) {
  const auto config = ow3d_directional::load_config(config_path);
  std::cout << ow3d_directional::config_summary_json(config) << "\n";
  return 0;
}

int run_generate(const std::filesystem::path& config_path) {
  const auto config = ow3d_directional::load_config(config_path);
  const auto preview = first_case_from_config(config);
  ow3d_directional::ensure_directory(config.output.output_root);
  ow3d_directional::write_batch_readme(
      config.output.output_root / "readme",
      config,
      preview.lx,
      preview.ly,
      preview.dx,
      preview.dy);

  int written = 0;
  for (double kd : config.physics.kd_list) {
    for (double akp : config.spectrum.akp_list) {
      for (double alpha : config.spectrum.alpha_list) {
        for (double t_init_periods : config.timing.t_init_periods_list) {
          const auto summary_params = ow3d_directional::build_case_parameters(
              config, kd, akp, alpha, config.spectrum.phases_deg.at(0), t_init_periods);
          std::cout << "\nGroup: "
                    << ow3d_directional::make_case_group_directory_name(config, summary_params).string()
                    << "\n";
          print_spectrum_summary(summary_params);

          std::vector<CaseParameters> phase_params;
          std::vector<mf12_cpp::Matrix> nonlinear_eta_phases;
          std::vector<mf12_cpp::Matrix> nonlinear_phi_phases;
          for (double phase_deg : config.spectrum.phases_deg) {
            const auto params = ow3d_directional::build_case_parameters(
                config, kd, akp, alpha, phase_deg, t_init_periods);
            const auto outputs = ow3d_directional::run_directional_case(config, params);
            const auto write_dir = config.output.output_root /
                                   ow3d_directional::make_case_directory_name(config, params);
            ow3d_directional::ensure_directory(write_dir);
            ow3d_directional::write_ow3d_init(write_dir / "OceanWave3D.init", outputs.nonlinear.eta, outputs.nonlinear.phi, params);
            ow3d_directional::write_ow3d_inp(write_dir / "OceanWave3D.inp", config, params);
            ow3d_directional::write_case_readme(write_dir / "OW_readme.txt", config, params);
            ow3d_directional::write_field_csv(write_dir / "eta_linear.csv", outputs.linear.eta);
            ow3d_directional::write_field_csv(write_dir / "phi_linear.csv", outputs.linear.phi);
            ow3d_directional::write_field_csv(write_dir / "eta_nonlinear.csv", outputs.nonlinear.eta);
            ow3d_directional::write_field_csv(write_dir / "phi_nonlinear.csv", outputs.nonlinear.phi);
            ow3d_directional::write_field_csv(write_dir / "x.csv", outputs.nonlinear.x);
            ow3d_directional::write_field_csv(write_dir / "y.csv", outputs.nonlinear.y);
            if (config.visualization.enabled) {
              ow3d_directional::write_visualizations(write_dir, config, params, outputs);
            }
            phase_params.push_back(params);
            nonlinear_eta_phases.push_back(outputs.nonlinear.eta);
            nonlinear_phi_phases.push_back(outputs.nonlinear.phi);
            ++written;
            std::cout << "wrote: " << write_dir.string() << "\n";
          }

          if (ow3d_directional::supports_standard_four_phase_separation(config.spectrum.phases_deg)) {
            const auto separated = ow3d_directional::separate_four_phase_fields(
                config.spectrum.phases_deg, nonlinear_eta_phases, nonlinear_phi_phases);
            const auto separation_dir = config.output.output_root / "four_phase_separation" /
                                        ow3d_directional::make_case_group_directory_name(config, summary_params);
            ow3d_directional::write_four_phase_separation(separation_dir, config, phase_params, separated);
            std::cout << "wrote four-phase separation: " << separation_dir.string() << "\n";
          } else if (config.spectrum.phases_deg.size() == 4U) {
            std::cout << "skipped four-phase separation because phases are not the standard [0, 90, 180, 270] set\n";
          }
        }
      }
    }
  }

  std::cout << "{\n"
            << "  \"output_root\": \"" << config.output.output_root.string() << "\",\n"
            << "  \"cases_written\": " << written << "\n"
            << "}\n";
  return 0;
}

int run_verify_matlab(const std::filesystem::path& config_path, const std::filesystem::path& reference_dir) {
  const auto config = ow3d_directional::load_config(config_path);
  const auto params = first_case_from_config(config);
  const auto outputs = ow3d_directional::run_directional_case(config, params);
  const auto reference = ow3d_directional::load_matlab_reference(reference_dir);
  const auto metrics = ow3d_directional::compare_to_reference(outputs.nonlinear.eta, outputs.nonlinear.phi, reference);
  const bool pass = ow3d_directional::metrics_within_default_tolerance(metrics);
  std::cout << "{\n"
            << "  \"pass\": " << (pass ? "true" : "false") << ",\n"
            << "  \"metrics\": {\n";
  bool first = true;
  for (const auto& [key, value] : metrics) {
    if (!first) {
      std::cout << ",\n";
    }
    std::cout << "    \"" << key << "\": " << std::setprecision(17) << value;
    first = false;
  }
  std::cout << "\n  }\n}\n";
  return pass ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 3) {
      print_usage();
      return 1;
    }

    const std::string command = argv[1];
    if (command == "inspect-config") {
      return run_inspect_config(argv[2]);
    }
    if (command == "generate") {
      return run_generate(argv[2]);
    }
    if (command == "verify-matlab") {
      if (argc < 4) {
        print_usage();
        return 1;
      }
      return run_verify_matlab(argv[2], argv[3]);
    }

    print_usage();
    return 1;
  } catch (const std::exception& exc) {
    std::cerr << "ow3d_directional_generator error: " << exc.what() << "\n";
    return 1;
  }
}
