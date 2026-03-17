#pragma once

#include "mf12_cpp/case_io.hpp"

#include <filesystem>
#include <string>
#include <vector>

namespace ow3d_directional {

struct PhysicsConfig {
  double g = 9.81;
  double kp = 0.00279;
  double tp = 12.0;
  std::vector<double> kd_list{1.0};
};

struct NonlinearConfig {
  int order = 3;
  std::string subharmonic_mode = "auto";
};

struct SpectrumConfig {
  std::string preset = "semi_gaussian_directional";
  std::vector<double> akp_list{0.12};
  std::vector<double> alpha_list{8.0};
  std::vector<double> phases_deg{0.0, 90.0, 180.0, 270.0};
  double heading_deg = 0.0;
  double spread_deg = 15.0;
  double energy_keep_frac = 0.999;
  int max_components = 800;
};

struct DomainConfig {
  double lx_lambda = 40.0;
  double ly_lambda = 15.0;
  int nx = 1025;
  int ny = 257;
};

struct TimingConfig {
  double t_focus = 0.0;
  std::vector<double> t_init_periods_list{-40.0, -30.0, -20.0};
  double t_end_periods = 5.0;
  int steps_per_period = 30;
  double focus_edge_padding_fraction = 0.05;
};

struct VisualizationConfig {
  bool enabled = true;
  bool write_linear_field = true;
  bool write_nonlinear_field = true;
  bool write_surface_style = true;
  std::string centerline_section = "mid";
  std::string format = "png";
  int width = 1400;
  int height = 900;
};

struct OutputConfig {
  std::filesystem::path output_root = std::filesystem::path("directional initial condition") / "cpp_generator";
  int store_surface_stride = 4;
  int surface_format = 1;
};

struct GeneratorConfig {
  PhysicsConfig physics;
  NonlinearConfig nonlinear;
  SpectrumConfig spectrum;
  DomainConfig domain;
  TimingConfig timing;
  VisualizationConfig visualization;
  OutputConfig output;
};

struct SpectrumDefinition {
  std::vector<double> kx;
  std::vector<double> ky;
  std::vector<double> amp;
  int n_components = 0;
  int candidate_components = 0;
  int threshold_filtered_components = 0;
  int energy_target_components = 0;
  double kmin = 0.0;
  double kmax = 0.0;
  double heading_deg = 0.0;
  double spread_deg = 0.0;
  double depth = 0.0;
  double energy_keep_frac = 0.0;
  double retained_weight_frac = 0.0;
};

struct CaseParameters {
  double kd = 0.0;
  double akp = 0.0;
  double alpha = 0.0;
  double phase_deg = 0.0;
  double tp = 0.0;
  double h = 0.0;
  double lx = 0.0;
  double ly = 0.0;
  double dx = 0.0;
  double dy = 0.0;
  double dt = 0.0;
  int n_steps = 0;
  double t_eval = 0.0;
  double t_init_periods = 0.0;
  double t_end_periods = 0.0;
  double duration_periods = 0.0;
  double focus_x = 0.0;
  double focus_y = 0.0;
  double focus_x_fraction = 0.0;
  SpectrumDefinition spectrum_definition;
  std::vector<double> a;
  std::vector<double> b;
};

struct MatlabReference {
  mf12_cpp::Matrix eta;
  mf12_cpp::Matrix phi;
  mf12_cpp::Matrix x;
  mf12_cpp::Matrix y;
};

GeneratorConfig load_config(const std::filesystem::path& path);
std::string config_summary_json(const GeneratorConfig& config);

}  // namespace ow3d_directional
