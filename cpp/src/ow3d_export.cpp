#include "ow3d_directional/ow3d_export.hpp"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace ow3d_directional {
namespace {

std::string now_stamp() {
  const auto now = std::chrono::system_clock::now();
  const std::time_t tt = std::chrono::system_clock::to_time_t(now);
  std::tm tm{};
#if defined(_WIN32)
  localtime_s(&tm, &tt);
#else
  localtime_r(&tt, &tm);
#endif
  std::ostringstream out;
  out << std::put_time(&tm, "%d-%b-%Y %H:%M:%S");
  return out.str();
}

mf12_cpp::Matrix transpose(const mf12_cpp::Matrix& input) {
  mf12_cpp::Matrix output;
  output.rows = input.cols;
  output.cols = input.rows;
  output.values.assign(input.values.size(), 0.0);
  for (std::size_t r = 0; r < input.rows; ++r) {
    for (std::size_t c = 0; c < input.cols; ++c) {
      output.values[c * output.cols + r] = input.values[r * input.cols + c];
    }
  }
  return output;
}

double max_value(const mf12_cpp::Matrix& matrix) {
  double value = matrix.values.empty() ? 0.0 : matrix.values.front();
  for (double candidate : matrix.values) {
    if (candidate > value) {
      value = candidate;
    }
  }
  return value;
}

}  // namespace

void ensure_directory(const std::filesystem::path& path) {
  std::filesystem::create_directories(path);
}

std::filesystem::path make_case_directory_name(
    const GeneratorConfig& config,
    const CaseParameters& params) {
  std::ostringstream name;
  name << make_case_group_directory_name(config, params).string()
       << "_phi_" << std::llround(params.phase_deg);
  return name.str();
}

std::filesystem::path make_case_group_directory_name(
    const GeneratorConfig& config,
    const CaseParameters& params) {
  std::ostringstream name;
  name << "T_init" << std::llround(params.t_init_periods)
       << "_Tend" << std::llround(params.t_end_periods)
       << "_Tp_kd" << std::fixed << std::setprecision(1) << params.kd
       << "_spread_" << std::llround(config.spectrum.spread_deg)
       << "_heading_" << std::llround(config.spectrum.heading_deg)
       << "_Akp_" << std::setw(3) << std::setfill('0') << std::llround(params.akp * 100.0)
       << "_alpha_" << std::fixed << std::setprecision(1) << params.alpha;
  return name.str();
}

void write_batch_readme(
    const std::filesystem::path& file_name,
    const GeneratorConfig& config,
    double lx,
    double ly,
    double dx,
    double dy) {
  std::ofstream out(file_name);
  if (!out) {
    throw std::runtime_error("Unable to write batch readme: " + file_name.string());
  }

  out << "Directional OW3D initial-condition batch\n";
  out << "Updated: " << now_stamp() << "\n\n";
  out << "Purpose\n";
  out << "Standalone C++ directional OW3D generation batch using the MF12 spectral implementation.\n\n";
  out << "Workflow\n";
  out << "- Full direct spectral reconstruction for each requested phase.\n";
  out << "- Surface-only OW3D export.\n";
  out << "- Quick-look visualization and field CSV dumps enabled by default.\n\n";
  out << "Current settings snapshot\n";
  out << "- output directory = " << config.output.output_root.string() << "\n";
  out << "- phases = [";
  for (std::size_t i = 0; i < config.spectrum.phases_deg.size(); ++i) {
    if (i > 0) {
      out << ", ";
    }
    out << config.spectrum.phases_deg[i];
  }
  out << "] deg\n";
  out << "- heading = " << config.spectrum.heading_deg << " deg\n";
  out << "- spread = " << config.spectrum.spread_deg << " deg\n";
  out << "- domain = [" << lx << ", " << ly << "] m\n";
  out << "- grid = [" << config.domain.nx << ", " << config.domain.ny << "], dx = " << dx << " m, dy = " << dy << " m\n";
  out << "- surface stride = " << config.output.store_surface_stride << "\n";
}

void write_ow3d_init(
    const std::filesystem::path& file_name,
    const mf12_cpp::Matrix& eta_in,
    const mf12_cpp::Matrix& phi_in,
    const CaseParameters& params) {
  const auto eta = transpose(eta_in);
  const auto phi = transpose(phi_in);

  std::ofstream out(file_name);
  if (!out) {
    throw std::runtime_error("Unable to write OceanWave3D.init: " + file_name.string());
  }

  out << std::fixed << std::setprecision(6);
  out << " H=" << max_value(eta)
      << " nx=" << eta.rows
      << " ny=" << eta.cols
      << " dx=" << params.dx
      << " dy=" << params.dy
      << " akp=" << params.akp
      << " shift=" << params.phase_deg << "\n";
  out << std::scientific << std::setprecision(6)
      << params.lx << " " << params.ly << " "
      << eta.rows << " " << eta.cols << " " << params.dt;

  for (std::size_t ry = 0; ry < eta.cols; ++ry) {
    for (std::size_t rx = 0; rx < eta.rows; ++rx) {
      const std::size_t index = rx * eta.cols + ry;
      out << "\n" << std::scientific << std::setprecision(6)
          << eta.values[index] << " " << phi.values[index];
    }
  }
}

void write_ow3d_inp(
    const std::filesystem::path& file_name,
    const GeneratorConfig& config,
    const CaseParameters& params) {
  std::ofstream out(file_name);
  if (!out) {
    throw std::runtime_error("Unable to write OceanWave3D.inp: " + file_name.string());
  }

  out << "Data for MF12 directional wave-group initialization " << now_stamp() << " <-\n";
  out << "-1 2 <-\n";
  out << std::llround(params.lx) << " "
      << std::llround(params.ly) << " "
      << std::llround(params.h) << " "
      << config.domain.nx << " "
      << config.domain.ny << " 17 0 0 1 1 1 1 <-\n";
  out << "4 4 4 1 1 1 <-\n";
  out << params.n_steps + 1 << " " << std::fixed << std::setprecision(6) << params.dt << " 1 0. 1 <-\n";
  out << "9.81 <-\n";
  out << "1 3 0 55 1e-6 1e-6 1 V 1 1 20 <-\n";
  out << "0.05 1.00 1.84 2 0 0 1 6 32 <-\n";
  out << config.output.store_surface_stride << " " << config.output.surface_format << " <-\n";
  out << "1 0 <-\n";
  out << "0 6 10 0.08 0.08 0.4 <-\n";
  out << "0 8. 3 X 0.0 <-\n";
  out << "0 0 <-\n";
  out << "0 2.0 2 0 0 1 0 <-\n";
  out << "0 <-\n";
  out << "33  8. 2. 80. 20. -1 -11 100. 50. run06.el 22.5 1.0 3.3 <-\n";
}

void write_case_readme(
    const std::filesystem::path& file_name,
    const GeneratorConfig& config,
    const CaseParameters& params) {
  std::ofstream out(file_name);
  if (!out) {
    throw std::runtime_error("Unable to write case readme: " + file_name.string());
  }

  out << "MF12 third-order directional wave-group initialization\n";
  out << "Generated: " << now_stamp() << "\n";
  out << std::fixed << std::setprecision(4)
      << "Akp=" << params.akp << ", Alpha=" << params.alpha
      << ", kp=" << std::setprecision(5) << config.physics.kp
      << ", kd=" << std::setprecision(2) << params.kd
      << ", h=" << std::setprecision(4) << params.h << "\n";
  out << std::setprecision(6)
      << "Tp=" << params.tp << " s, t_init=" << params.t_eval << " s, phase shift=" << params.phase_deg << " deg\n";
  out << "Heading=" << params.spectrum_definition.heading_deg
      << " deg, spreading sigma=" << params.spectrum_definition.spread_deg << " deg\n";
  out << "Components kept=" << params.spectrum_definition.n_components
      << ", candidate bins=" << params.spectrum_definition.candidate_components
      << ", after weight filter=" << params.spectrum_definition.threshold_filtered_components
      << ", energy-target count=" << params.spectrum_definition.energy_target_components << "\n";
  out << "Energy keep request=" << std::setprecision(5) << params.spectrum_definition.energy_keep_frac
      << ", actual retained weight fraction=" << params.spectrum_definition.retained_weight_frac << "\n";
  out << "k-range=[" << std::setprecision(6) << params.spectrum_definition.kmin
      << ", " << params.spectrum_definition.kmax << "] 1/m\n";
  out << "Focus point: x=" << params.focus_x
      << " m (" << std::setprecision(4) << params.focus_x_fraction
      << " Lx), y=" << params.focus_y << " m\n";
  out << "Grid: Nx=" << config.domain.nx << ", Ny=" << config.domain.ny
      << ", dx=" << std::setprecision(6) << params.dx
      << " m, dy=" << params.dy << " m\n";
  out << "Model: MF12 spectral coefficients/order-3, full direct reconstruction for each phase\n";
  out << "Initial condition time relative to focus: " << std::setprecision(2) << params.t_init_periods << " Tp\n";
  out << "Target end time relative to focus: " << params.t_end_periods << " Tp\n";
  out << "Total OW3D duration after initialization: " << params.duration_periods << " Tp\n";
  out << "Surface output stride: every " << std::abs(config.output.store_surface_stride) << " time step(s)\n";
  out << "Kinematic output: disabled\n";
}

void write_field_csv(
    const std::filesystem::path& file_name,
    const mf12_cpp::Matrix& matrix) {
  std::ofstream out(file_name);
  if (!out) {
    throw std::runtime_error("Unable to write CSV: " + file_name.string());
  }
  out << std::setprecision(17);
  for (std::size_t r = 0; r < matrix.rows; ++r) {
    for (std::size_t c = 0; c < matrix.cols; ++c) {
      if (c > 0) {
        out << ",";
      }
      out << matrix.values[r * matrix.cols + c];
    }
    out << "\n";
  }
}

}  // namespace ow3d_directional
