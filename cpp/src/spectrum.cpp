#include "ow3d_directional/spectrum.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace ow3d_directional {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

double deg2rad(double degrees) {
  return degrees * kPi / 180.0;
}

double wrap_angle(double angle) {
  return std::atan2(std::sin(angle), std::cos(angle));
}

double gaussian_spreading(double theta, double spread_angle_deg) {
  const double sigma = std::max(deg2rad(spread_angle_deg), 1e-12);
  const double wrapped = wrap_angle(theta);
  return std::exp(-0.5 * (wrapped / sigma) * (wrapped / sigma));
}

struct WeightedComponent {
  double weight = 0.0;
  double kx = 0.0;
  double ky = 0.0;
  double kmag = 0.0;
};

}  // namespace

void resolve_focus_point(
    const GeneratorConfig& config,
    double lx,
    double ly,
    double t_init_periods,
    double* xf,
    double* yf,
    double* focus_x_fraction) {
  const double duration_periods = config.timing.t_end_periods - t_init_periods;
  if (duration_periods <= 0.0) {
    throw std::runtime_error("t_end_periods must be larger than t_init_periods.");
  }
  double focus_fraction = -t_init_periods / duration_periods;
  const double pad = config.timing.focus_edge_padding_fraction;
  focus_fraction = pad + (1.0 - 2.0 * pad) * focus_fraction;
  focus_fraction = std::min(std::max(focus_fraction, pad), 1.0 - pad);
  *focus_x_fraction = focus_fraction;
  *xf = focus_fraction * lx;
  *yf = 0.5 * ly;
}

SpectrumDefinition build_directional_group_spectrum(
    const GeneratorConfig& config,
    double lx,
    double ly,
    double h,
    double akp,
    double alpha) {
  if (config.spectrum.preset != "semi_gaussian_directional") {
    throw std::runtime_error("Unsupported spectrum preset: " + config.spectrum.preset);
  }

  const double dkx = 2.0 * kPi / lx;
  const double dky = 2.0 * kPi / ly;

  const int nx = config.domain.nx;
  const int ny = config.domain.ny;

  std::vector<WeightedComponent> components;
  components.reserve(static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny) / 2U);

  const double kw_left = 0.004606;
  const double kw_right = std::sqrt((config.physics.kp * config.physics.kp) /
                                    (2.0 * std::log(std::pow(10.0, alpha))));

  for (int jy = -ny / 2; jy < (ny + 1) / 2; ++jy) {
    for (int ix = -nx / 2; ix < (nx + 1) / 2; ++ix) {
      const double kx = static_cast<double>(ix) * dkx;
      const double ky = static_cast<double>(jy) * dky;
      if (!(kx > 0.0 || (kx == 0.0 && ky > 0.0))) {
        continue;
      }

      const double kmag = std::hypot(kx, ky);
      const double theta = std::atan2(ky, kx);
      const double kw = kmag <= config.physics.kp ? kw_left : kw_right;
      const double sk = std::exp(-((kmag - config.physics.kp) * (kmag - config.physics.kp)) / (2.0 * kw * kw));
      const double d = gaussian_spreading(theta - deg2rad(config.spectrum.heading_deg), config.spectrum.spread_deg);
      const double weight = sk * d;
      components.push_back({weight, kx, ky, kmag});
    }
  }

  const double max_weight = std::max_element(
      components.begin(), components.end(),
      [](const WeightedComponent& lhs, const WeightedComponent& rhs) {
        return lhs.weight < rhs.weight;
      })->weight;

  std::vector<WeightedComponent> filtered;
  filtered.reserve(components.size());
  for (const auto& component : components) {
    if (component.weight > 1e-10 * max_weight) {
      filtered.push_back(component);
    }
  }

  std::sort(filtered.begin(), filtered.end(), [](const WeightedComponent& lhs, const WeightedComponent& rhs) {
    return lhs.weight > rhs.weight;
  });

  double total_weight = 0.0;
  for (const auto& component : filtered) {
    total_weight += component.weight;
  }

  const double target_weight = config.spectrum.energy_keep_frac * total_weight;
  int n_keep = 0;
  double cumulative = 0.0;
  for (const auto& component : filtered) {
    cumulative += component.weight;
    ++n_keep;
    if (cumulative >= target_weight) {
      break;
    }
  }
  const int n_keep_energy = n_keep;
  n_keep = std::min(n_keep, config.spectrum.max_components);

  SpectrumDefinition definition;
  definition.kx.resize(n_keep);
  definition.ky.resize(n_keep);
  definition.amp.resize(n_keep);
  definition.n_components = n_keep;
  definition.candidate_components = static_cast<int>(components.size());
  definition.threshold_filtered_components = static_cast<int>(filtered.size());
  definition.energy_target_components = n_keep_energy;
  definition.heading_deg = config.spectrum.heading_deg;
  definition.spread_deg = config.spectrum.spread_deg;
  definition.depth = h;
  definition.energy_keep_frac = config.spectrum.energy_keep_frac;

  double amp_sum = 0.0;
  double retained_weight = 0.0;
  for (int i = 0; i < n_keep; ++i) {
    definition.kx[i] = filtered[i].kx;
    definition.ky[i] = filtered[i].ky;
    definition.amp[i] = filtered[i].weight;
    amp_sum += definition.amp[i];
    retained_weight += filtered[i].weight;
  }
  definition.retained_weight_frac = total_weight > 0.0 ? retained_weight / total_weight : 0.0;

  const double target_amp_sum = akp / config.physics.kp;
  const double scale = target_amp_sum / std::max(amp_sum, 1e-12);
  for (double& value : definition.amp) {
    value *= scale;
  }

  auto [min_it, max_it] = std::minmax_element(filtered.begin(), filtered.begin() + n_keep,
      [](const WeightedComponent& lhs, const WeightedComponent& rhs) {
        return lhs.kmag < rhs.kmag;
      });
  definition.kmin = min_it->kmag;
  definition.kmax = max_it->kmag;

  return definition;
}

CaseParameters build_case_parameters(
    const GeneratorConfig& config,
    double kd,
    double akp,
    double alpha,
    double phase_deg,
    double t_init_periods) {
  CaseParameters params;
  params.kd = kd;
  params.akp = akp;
  params.alpha = alpha;
  params.phase_deg = phase_deg;
  params.tp = config.physics.tp;
  params.h = kd / config.physics.kp;

  const double lambda_p = 2.0 * kPi / config.physics.kp;
  params.lx = config.domain.lx_lambda * lambda_p;
  params.ly = config.domain.ly_lambda * lambda_p;
  params.dx = params.lx / static_cast<double>(config.domain.nx);
  params.dy = params.ly / static_cast<double>(config.domain.ny);
  params.dt = params.tp / static_cast<double>(config.timing.steps_per_period);
  params.t_eval = t_init_periods * params.tp;
  params.t_init_periods = t_init_periods;
  params.t_end_periods = config.timing.t_end_periods;
  params.duration_periods = params.t_end_periods - params.t_init_periods;
  params.n_steps = static_cast<int>(std::llround(params.duration_periods * config.timing.steps_per_period));

  resolve_focus_point(
      config, params.lx, params.ly, t_init_periods,
      &params.focus_x, &params.focus_y, &params.focus_x_fraction);

  params.spectrum_definition = build_directional_group_spectrum(
      config, params.lx, params.ly, params.h, akp, alpha);

  const double phase_shift = deg2rad(phase_deg);
  const std::size_t n = params.spectrum_definition.kx.size();
  params.a.resize(n);
  params.b.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    const double kx = params.spectrum_definition.kx[i];
    const double ky = params.spectrum_definition.ky[i];
    const double k = std::hypot(kx, ky);
    const double omega_lin = std::sqrt(config.physics.g * k * std::tanh(params.h * k));
    const double phase_focus = -(kx * params.focus_x + ky * params.focus_y) +
                               omega_lin * config.timing.t_focus + phase_shift;
    params.a[i] = params.spectrum_definition.amp[i] * std::cos(phase_focus);
    params.b[i] = params.spectrum_definition.amp[i] * std::sin(phase_focus);
  }

  return params;
}

}  // namespace ow3d_directional
