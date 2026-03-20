#include "ow3d_directional/four_phase.hpp"

#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

#if defined(MF12_HAVE_FFTW)
#include <fftw3.h>
#endif

namespace ow3d_directional {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

double normalize_phase_deg(double phase_deg) {
  double value = std::fmod(phase_deg, 360.0);
  if (value < 0.0) {
    value += 360.0;
  }
  if (std::abs(value - 360.0) < 1e-9) {
    value = 0.0;
  }
  return value;
}

bool approx_equal(double lhs, double rhs) {
  return std::abs(lhs - rhs) < 1e-6;
}

std::array<int, 4> phase_index_map(const std::vector<double>& phases_deg) {
  if (phases_deg.size() != 4U) {
    throw std::runtime_error("four-phase separation requires exactly four phases.");
  }
  std::array<int, 4> indices{-1, -1, -1, -1};
  for (std::size_t i = 0; i < phases_deg.size(); ++i) {
    const double value = normalize_phase_deg(phases_deg[i]);
    if (approx_equal(value, 0.0)) {
      indices[0] = static_cast<int>(i);
    } else if (approx_equal(value, 90.0)) {
      indices[1] = static_cast<int>(i);
    } else if (approx_equal(value, 180.0)) {
      indices[2] = static_cast<int>(i);
    } else if (approx_equal(value, 270.0)) {
      indices[3] = static_cast<int>(i);
    }
  }
  for (int idx : indices) {
    if (idx < 0) {
      throw std::runtime_error("four-phase separation requires phases 0, 90, 180, 270 degrees.");
    }
  }
  return indices;
}

void assert_same_shape(const std::vector<mf12_cpp::Matrix>& fields) {
  if (fields.empty()) {
    throw std::runtime_error("four-phase separation received no fields.");
  }
  const std::size_t rows = fields.front().rows;
  const std::size_t cols = fields.front().cols;
  for (const auto& field : fields) {
    if (field.rows != rows || field.cols != cols) {
      throw std::runtime_error("four-phase separation requires all fields to share the same dimensions.");
    }
  }
}

std::vector<std::complex<double>> fft_1d(const std::vector<std::complex<double>>& in) {
  const std::size_t n = in.size();
#if defined(MF12_HAVE_FFTW)
  std::vector<std::complex<double>> out(n);
  auto* src = reinterpret_cast<fftw_complex*>(const_cast<std::complex<double>*>(in.data()));
  auto* dst = reinterpret_cast<fftw_complex*>(out.data());
  fftw_plan plan = fftw_plan_dft_1d(static_cast<int>(n), src, dst, FFTW_FORWARD, FFTW_ESTIMATE);
  if (plan == nullptr) {
    throw std::runtime_error("Failed to create FFTW forward plan for four-phase separation.");
  }
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  return out;
#else
  std::vector<std::complex<double>> out(n, std::complex<double>(0.0, 0.0));
  for (std::size_t k = 0; k < n; ++k) {
    for (std::size_t x = 0; x < n; ++x) {
      const double angle = -2.0 * kPi * static_cast<double>(k * x) / static_cast<double>(n);
      out[k] += in[x] * std::exp(std::complex<double>(0.0, angle));
    }
  }
  return out;
#endif
}

std::vector<std::complex<double>> ifft_1d(const std::vector<std::complex<double>>& in) {
  const std::size_t n = in.size();
#if defined(MF12_HAVE_FFTW)
  std::vector<std::complex<double>> out(n);
  auto* src = reinterpret_cast<fftw_complex*>(const_cast<std::complex<double>*>(in.data()));
  auto* dst = reinterpret_cast<fftw_complex*>(out.data());
  fftw_plan plan = fftw_plan_dft_1d(static_cast<int>(n), src, dst, FFTW_BACKWARD, FFTW_ESTIMATE);
  if (plan == nullptr) {
    throw std::runtime_error("Failed to create FFTW inverse plan for four-phase separation.");
  }
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  const double inv_n = 1.0 / static_cast<double>(n);
  for (auto& value : out) {
    value *= inv_n;
  }
  return out;
#else
  std::vector<std::complex<double>> out(n, std::complex<double>(0.0, 0.0));
  for (std::size_t x = 0; x < n; ++x) {
    for (std::size_t k = 0; k < n; ++k) {
      const double angle = 2.0 * kPi * static_cast<double>(k * x) / static_cast<double>(n);
      out[x] += in[k] * std::exp(std::complex<double>(0.0, angle));
    }
    out[x] /= static_cast<double>(n);
  }
  return out;
#endif
}

mf12_cpp::Matrix negative_imag_hilbert_like(const mf12_cpp::Matrix& field) {
  mf12_cpp::Matrix out;
  out.rows = field.rows;
  out.cols = field.cols;
  out.values.assign(field.values.size(), 0.0);

  const std::size_t rows = field.rows;
  const std::size_t cols = field.cols;
  const std::size_t floor_half = rows / 2U;
  const std::size_t ceil_half = (rows + 1U) / 2U;

  for (std::size_t c = 0; c < cols; ++c) {
    std::vector<std::complex<double>> column(rows);
    for (std::size_t r = 0; r < rows; ++r) {
      column[r] = std::complex<double>(field.values[r * cols + c], 0.0);
    }
    auto spectrum = fft_1d(column);
    for (std::size_t k = 0; k < rows; ++k) {
      const std::size_t matlab_idx = k + 1U;
      double mask = 0.0;
      if (matlab_idx >= ceil_half) {
        mask = -1.0;
      } else if (matlab_idx <= floor_half) {
        mask = 1.0;
      }
      spectrum[k] *= mask;
    }
    const auto transformed = ifft_1d(spectrum);
    for (std::size_t r = 0; r < rows; ++r) {
      out.values[r * cols + c] = -std::imag(transformed[r]);
    }
  }

  return out;
}

mf12_cpp::Matrix linear_combine(
    const std::vector<mf12_cpp::Matrix>& fields,
    const std::array<double, 8>& weights) {
  mf12_cpp::Matrix out;
  out.rows = fields.front().rows;
  out.cols = fields.front().cols;
  out.values.assign(fields.front().values.size(), 0.0);
  for (std::size_t i = 0; i < fields.size(); ++i) {
    const double w = weights[i];
    for (std::size_t j = 0; j < out.values.size(); ++j) {
      out.values[j] += w * fields[i].values[j];
    }
  }
  return out;
}

}  // namespace

bool supports_standard_four_phase_separation(const std::vector<double>& phases_deg) {
  try {
    (void) phase_index_map(phases_deg);
    return true;
  } catch (const std::exception&) {
    return false;
  }
}

FourPhaseSeparation separate_four_phase_fields(
    const std::vector<double>& phases_deg,
    const std::vector<mf12_cpp::Matrix>& eta_phases,
    const std::vector<mf12_cpp::Matrix>& phi_phases) {
  const auto indices = phase_index_map(phases_deg);
  assert_same_shape(eta_phases);
  assert_same_shape(phi_phases);

  std::vector<mf12_cpp::Matrix> eta_ordered(4);
  std::vector<mf12_cpp::Matrix> phi_ordered(4);
  for (int i = 0; i < 4; ++i) {
    eta_ordered[static_cast<std::size_t>(i)] = eta_phases[static_cast<std::size_t>(indices[static_cast<std::size_t>(i)])];
    phi_ordered[static_cast<std::size_t>(i)] = phi_phases[static_cast<std::size_t>(indices[static_cast<std::size_t>(i)])];
  }

  std::vector<mf12_cpp::Matrix> eta_all = eta_ordered;
  std::vector<mf12_cpp::Matrix> phi_all = phi_ordered;
  for (const auto& field : eta_ordered) {
    eta_all.push_back(negative_imag_hilbert_like(field));
  }
  for (const auto& field : phi_ordered) {
    phi_all.push_back(negative_imag_hilbert_like(field));
  }

  const std::array<std::array<double, 8>, 4> coef{{
      {{0.25, 0.0, -0.25, 0.0, 0.0, -0.25, 0.0, 0.25}},
      {{0.25, -0.25, 0.25, -0.25, 0.0, 0.0, 0.0, 0.0}},
      {{0.25, 0.0, -0.25, 0.0, 0.0, 0.25, 0.0, -0.25}},
      {{0.25, 0.25, 0.25, 0.25, 0.0, 0.0, 0.0, 0.0}},
  }};

  FourPhaseSeparation separated;
  separated.eta1 = linear_combine(eta_all, coef[0]);
  separated.eta2 = linear_combine(eta_all, coef[1]);
  separated.eta3 = linear_combine(eta_all, coef[2]);
  separated.eta4 = linear_combine(eta_all, coef[3]);
  separated.phi1 = linear_combine(phi_all, coef[0]);
  separated.phi2 = linear_combine(phi_all, coef[1]);
  separated.phi3 = linear_combine(phi_all, coef[2]);
  separated.phi4 = linear_combine(phi_all, coef[3]);
  return separated;
}

}  // namespace ow3d_directional
