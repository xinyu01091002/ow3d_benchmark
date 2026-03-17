#include "ow3d_directional/config.hpp"
#include "ow3d_directional/mf12_runner.hpp"
#include "ow3d_directional/spectrum.hpp"

#include <cmath>
#include <stdexcept>

namespace {

void expect(bool condition, const char* message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void test_focus_point() {
  ow3d_directional::GeneratorConfig config;
  double xf = 0.0;
  double yf = 0.0;
  double frac = 0.0;
  ow3d_directional::resolve_focus_point(config, 100.0, 20.0, -20.0, &xf, &yf, &frac);
  expect(xf > 0.0 && xf < 100.0, "focus x should lie inside the domain");
  expect(std::abs(yf - 10.0) < 1e-12, "focus y should be at the centerline");
  expect(frac > 0.0 && frac < 1.0, "focus fraction should lie inside the domain");
}

void test_spectrum_builder() {
  ow3d_directional::GeneratorConfig config;
  config.domain.nx = 65;
  config.domain.ny = 33;
  const auto spectrum = ow3d_directional::build_directional_group_spectrum(config, 1000.0, 500.0, 100.0, 0.12, 8.0);
  expect(!spectrum.kx.empty(), "spectrum should retain components");
  expect(spectrum.kx.size() == spectrum.ky.size(), "kx and ky should match");
  expect(spectrum.kx.size() == spectrum.amp.size(), "amplitudes should match component count");
  double amp_sum = 0.0;
  for (double value : spectrum.amp) {
    expect(std::isfinite(value), "amplitudes should be finite");
    expect(value > 0.0, "amplitudes should be positive");
    amp_sum += value;
  }
  expect(amp_sum > 0.0, "retained amplitude sum should be positive");
}

void test_case_parameters() {
  ow3d_directional::GeneratorConfig config;
  const auto params = ow3d_directional::build_case_parameters(config, 1.0, 0.12, 8.0, 0.0, -40.0);
  expect(params.n_steps > 0, "n_steps should be positive");
  expect(params.a.size() == params.b.size(), "a and b should match");
  expect(params.a.size() == params.spectrum_definition.kx.size(), "component arrays should match spectrum");
}

void test_linear_reconstruction_matches_direct_sum() {
  ow3d_directional::GeneratorConfig config;
  config.physics.kp = 0.0279;
  config.physics.tp = 12.0;
  config.physics.kd_list = {1.0};
  config.spectrum.akp_list = {0.06};
  config.spectrum.alpha_list = {8.0};
  config.spectrum.energy_keep_frac = 0.995;
  config.spectrum.max_components = 40;
  config.domain.nx = 65;
  config.domain.ny = 33;
  config.domain.lx_lambda = 5.0;
  config.domain.ly_lambda = 3.0;
  config.timing.t_init_periods_list = {-2.0};
  config.timing.t_end_periods = 1.0;
  config.timing.steps_per_period = 12;

  const auto params = ow3d_directional::build_case_parameters(config, 1.0, 0.06, 8.0, 0.0, -2.0);
  const auto outputs = ow3d_directional::run_directional_case(config, params);

  const std::size_t ix = 26;
  const std::size_t iy = 16;
  const double x = static_cast<double>(ix) * params.dx;
  const double y = static_cast<double>(iy) * params.dy;
  double direct = 0.0;
  for (std::size_t i = 0; i < params.a.size(); ++i) {
    const double kx = params.spectrum_definition.kx[i];
    const double ky = params.spectrum_definition.ky[i];
    const double k = std::hypot(kx, ky);
    const double omega = std::sqrt(config.physics.g * k * std::tanh(params.h * k));
    const double theta = omega * params.t_eval - kx * x - ky * y;
    direct += params.a[i] * std::cos(theta) + params.b[i] * std::sin(theta);
  }
  const double reconstructed =
      outputs.linear.eta.values[iy * outputs.linear.eta.cols + ix];
  expect(std::abs(reconstructed - direct) < 1e-6, "linear reconstruction should match direct sum");
}

}  // namespace

int main() {
  test_focus_point();
  test_spectrum_builder();
  test_case_parameters();
  test_linear_reconstruction_matches_direct_sum();
  return 0;
}
