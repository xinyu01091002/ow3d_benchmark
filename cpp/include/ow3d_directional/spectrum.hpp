#pragma once

#include "ow3d_directional/config.hpp"

namespace ow3d_directional {

SpectrumDefinition build_directional_group_spectrum(
    const GeneratorConfig& config,
    double lx,
    double ly,
    double h,
    double akp,
    double alpha);

void resolve_focus_point(
    const GeneratorConfig& config,
    double lx,
    double ly,
    double t_init_periods,
    double* xf,
    double* yf,
    double* focus_x_fraction);

CaseParameters build_case_parameters(
    const GeneratorConfig& config,
    double kd,
    double akp,
    double alpha,
    double phase_deg,
    double t_init_periods);

}  // namespace ow3d_directional
