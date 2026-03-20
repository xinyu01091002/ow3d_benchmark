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
using ow3d_directional::CombinedCaseConfig;
using ow3d_directional::GeneratorConfig;

struct PlannedCaseGroup {
  double kd = 0.0;
  double akp = 0.0;
  double alpha = 0.0;
  double t_init_periods = 0.0;
  std::vector<double> phases_deg;
};

std::size_t pair_count_for(std::size_t n_comp) {
  return n_comp < 2U ? 0U : n_comp * (n_comp - 1U) / 2U;
}

std::size_t triplet_count_for(std::size_t n_comp) {
  return n_comp < 3U ? 0U : n_comp * (n_comp - 1U) * (n_comp - 2U) / 6U;
}

std::vector<PlannedCaseGroup> planned_case_groups(const GeneratorConfig& config) {
  std::vector<PlannedCaseGroup> groups;
  if (!config.combined_cases.empty()) {
    groups.reserve(config.combined_cases.size());
    for (const CombinedCaseConfig& combined : config.combined_cases) {
      PlannedCaseGroup group;
      group.kd = combined.kd;
      group.akp = combined.akp;
      group.alpha = combined.alpha;
      group.t_init_periods = combined.t_init_periods;
      group.phases_deg = combined.phases_deg.empty() ? config.spectrum.phases_deg : combined.phases_deg;
      groups.push_back(std::move(group));
    }
    return groups;
  }

  for (double kd : config.physics.kd_list) {
    for (double akp : config.spectrum.akp_list) {
      for (double alpha : config.spectrum.alpha_list) {
        for (double t_init_periods : config.timing.t_init_periods_list) {
          PlannedCaseGroup group;
          group.kd = kd;
          group.akp = akp;
          group.alpha = alpha;
          group.t_init_periods = t_init_periods;
          group.phases_deg = config.spectrum.phases_deg;
          groups.push_back(std::move(group));
        }
      }
    }
  }
  return groups;
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
  const auto groups = planned_case_groups(config);
  if (groups.empty() || groups.front().phases_deg.empty()) {
    throw std::runtime_error("No parameter sets available to generate.");
  }
  return ow3d_directional::build_case_parameters(
      config,
      groups.front().kd,
      groups.front().akp,
      groups.front().alpha,
      groups.front().phases_deg.front(),
      groups.front().t_init_periods);
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
  const auto groups = planned_case_groups(config);
  for (const PlannedCaseGroup& group : groups) {
    if (group.phases_deg.empty()) {
      throw std::runtime_error("Each generated case group must contain at least one phase.");
    }
    const auto summary_params = ow3d_directional::build_case_parameters(
        config, group.kd, group.akp, group.alpha, group.phases_deg.front(), group.t_init_periods);
          std::cout << "\nGroup: "
                    << ow3d_directional::make_case_group_directory_name(config, summary_params).string()
                    << "\n";
          print_spectrum_summary(summary_params);

          std::vector<CaseParameters> phase_params;
          std::vector<mf12_cpp::Matrix> nonlinear_eta_phases;
          std::vector<mf12_cpp::Matrix> nonlinear_phi_phases;
          for (double phase_deg : group.phases_deg) {
            const auto params = ow3d_directional::build_case_parameters(
                config, group.kd, group.akp, group.alpha, phase_deg, group.t_init_periods);
            const auto outputs = ow3d_directional::run_directional_case(config, params);
            const auto write_dir = config.output.output_root /
                                   ow3d_directional::make_case_directory_name(config, params);
            ow3d_directional::ensure_directory(write_dir);
            ow3d_directional::write_ow3d_init(write_dir / "OceanWave3D.init", outputs.nonlinear.eta, outputs.nonlinear.phi, params);
            ow3d_directional::write_ow3d_inp(write_dir / "OceanWave3D.inp", config, params);
            ow3d_directional::write_case_readme(write_dir / "OW_readme.txt", config, params);
            if (config.visualization.enabled) {
              ow3d_directional::write_visualizations(write_dir, config, params, outputs);
            }
            phase_params.push_back(params);
            nonlinear_eta_phases.push_back(outputs.nonlinear.eta);
            nonlinear_phi_phases.push_back(outputs.nonlinear.phi);
            ++written;
            std::cout << "wrote: " << write_dir.string() << "\n";
          }

          if (!ow3d_directional::supports_standard_four_phase_separation(group.phases_deg) &&
              group.phases_deg.size() == 4U) {
            std::cout << "skipped four-phase separation because phases are not the standard [0, 90, 180, 270] set\n";
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
