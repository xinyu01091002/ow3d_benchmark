#include "ow3d_directional/mf12_runner.hpp"

#include <map>
#include <string>

namespace ow3d_directional {
namespace {

mf12_cpp::Matrix make_vector_matrix(const std::vector<double>& values) {
  mf12_cpp::Matrix matrix;
  matrix.rows = values.size();
  matrix.cols = 1;
  matrix.values = values;
  return matrix;
}

}  // namespace

mf12_cpp::LoadedCase build_loaded_case(
    int order,
    const GeneratorConfig& config,
    const CaseParameters& params) {
  mf12_cpp::LoadedCase loaded;
  loaded.manifest.case_id = "ow3d_directional_generated";
  loaded.manifest.description = "Directional OW3D generation case";
  loaded.manifest.purpose = "MF12 spectral reconstruction for OW3D export";
  loaded.manifest.inputs.order = order;
  loaded.manifest.inputs.g = config.physics.g;
  loaded.manifest.inputs.h = params.h;
  loaded.manifest.inputs.Ux = 0.0;
  loaded.manifest.inputs.Uy = 0.0;
  loaded.manifest.inputs.Lx = params.lx;
  loaded.manifest.inputs.Ly = params.ly;
  loaded.manifest.inputs.Nx = config.domain.nx;
  loaded.manifest.inputs.Ny = config.domain.ny;
  loaded.manifest.inputs.t = params.t_eval;
  loaded.manifest.inputs.subharmonic_mode = config.nonlinear.subharmonic_mode;

  loaded.arrays["a"] = make_vector_matrix(params.a);
  loaded.arrays["b"] = make_vector_matrix(params.b);
  loaded.arrays["kx"] = make_vector_matrix(params.spectrum_definition.kx);
  loaded.arrays["ky"] = make_vector_matrix(params.spectrum_definition.ky);
  return loaded;
}

Mf12RunOutputs run_directional_case(const GeneratorConfig& config, const CaseParameters& params) {
  Mf12RunOutputs outputs;
  outputs.linear = mf12_cpp::run_case(build_loaded_case(1, config, params), 1, false);
  outputs.nonlinear = mf12_cpp::run_case(build_loaded_case(config.nonlinear.order, config, params), 1, false);
  return outputs;
}

}  // namespace ow3d_directional
