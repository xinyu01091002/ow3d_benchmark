#include "ow3d_directional/verification.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace ow3d_directional {
namespace {

double vector_norm(const std::vector<double>& values) {
  double sum = 0.0;
  for (double value : values) {
    sum += value * value;
  }
  return std::sqrt(sum);
}

std::map<std::string, double> compare_fields(
    const mf12_cpp::Matrix& candidate,
    const mf12_cpp::Matrix& reference,
    const std::string& prefix) {
  if (candidate.rows != reference.rows || candidate.cols != reference.cols) {
    throw std::runtime_error("Field dimensions do not match for comparison.");
  }
  std::vector<double> diff(candidate.values.size(), 0.0);
  double max_abs = 0.0;
  double sum_sq = 0.0;
  for (std::size_t i = 0; i < candidate.values.size(); ++i) {
    diff[i] = candidate.values[i] - reference.values[i];
    const double abs_val = std::abs(diff[i]);
    if (abs_val > max_abs) {
      max_abs = abs_val;
    }
    sum_sq += diff[i] * diff[i];
  }
  const double denom = std::max(vector_norm(reference.values), std::numeric_limits<double>::epsilon());
  return {
      {prefix + "_max_abs_err", max_abs},
      {prefix + "_rms_err", std::sqrt(sum_sq / static_cast<double>(diff.size()))},
      {prefix + "_relative_l2_err", vector_norm(diff) / denom},
  };
}

}  // namespace

MatlabReference load_matlab_reference(const std::filesystem::path& reference_dir) {
  MatlabReference reference;
  reference.eta = mf12_cpp::load_csv_matrix(reference_dir / "eta.csv");
  reference.phi = mf12_cpp::load_csv_matrix(reference_dir / "phi.csv");
  reference.x = mf12_cpp::load_csv_matrix(reference_dir / "x.csv");
  reference.y = mf12_cpp::load_csv_matrix(reference_dir / "y.csv");
  return reference;
}

std::map<std::string, double> compare_to_reference(
    const mf12_cpp::Matrix& eta,
    const mf12_cpp::Matrix& phi,
    const MatlabReference& reference) {
  std::map<std::string, double> metrics;
  const auto eta_metrics = compare_fields(eta, reference.eta, "eta");
  const auto phi_metrics = compare_fields(phi, reference.phi, "phi");
  metrics.insert(eta_metrics.begin(), eta_metrics.end());
  metrics.insert(phi_metrics.begin(), phi_metrics.end());
  return metrics;
}

bool metrics_within_default_tolerance(const std::map<std::string, double>& metrics) {
  constexpr double max_abs_tol = 1e-8;
  constexpr double rel_l2_tol = 1e-8;
  for (const auto& [key, value] : metrics) {
    if (key.find("max_abs_err") != std::string::npos && value > max_abs_tol) {
      return false;
    }
    if (key.find("relative_l2_err") != std::string::npos && value > rel_l2_tol) {
      return false;
    }
  }
  return true;
}

}  // namespace ow3d_directional
