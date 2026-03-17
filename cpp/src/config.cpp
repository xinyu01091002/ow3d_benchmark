#include "ow3d_directional/config.hpp"

#include <nlohmann/json.hpp>

#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace ow3d_directional {
namespace {

using json = nlohmann::json;

json read_json(const std::filesystem::path& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("Unable to open config file: " + path.string());
  }
  json parsed;
  in >> parsed;
  return parsed;
}

template <typename T>
std::vector<T> read_vector_or_scalar(const json& node, const char* key, const std::vector<T>& defaults) {
  if (!node.contains(key)) {
    return defaults;
  }
  const auto& value = node.at(key);
  if (value.is_array()) {
    return value.get<std::vector<T>>();
  }
  return {value.get<T>()};
}

template <typename T>
T read_with_default(const json& node, const char* key, const T& default_value) {
  if (!node.contains(key)) {
    return default_value;
  }
  return node.at(key).get<T>();
}

}  // namespace

GeneratorConfig load_config(const std::filesystem::path& path) {
  const json root = read_json(path);
  GeneratorConfig config;

  if (root.contains("physics")) {
    const auto& node = root.at("physics");
    config.physics.g = read_with_default(node, "g", config.physics.g);
    config.physics.kp = read_with_default(node, "kp", config.physics.kp);
    config.physics.tp = read_with_default(node, "Tp", config.physics.tp);
    config.physics.kd_list = read_vector_or_scalar<double>(node, "kd", config.physics.kd_list);
    if (node.contains("kd_list")) {
      config.physics.kd_list = node.at("kd_list").get<std::vector<double>>();
    }
  }

  if (root.contains("nonlinear")) {
    const auto& node = root.at("nonlinear");
    config.nonlinear.order = read_with_default(node, "order", config.nonlinear.order);
    config.nonlinear.subharmonic_mode = read_with_default(node, "subharmonic_mode", config.nonlinear.subharmonic_mode);
  }

  if (root.contains("spectrum")) {
    const auto& node = root.at("spectrum");
    config.spectrum.preset = read_with_default(node, "preset", config.spectrum.preset);
    config.spectrum.akp_list = read_vector_or_scalar<double>(node, "Akp", config.spectrum.akp_list);
    config.spectrum.alpha_list = read_vector_or_scalar<double>(node, "Alpha", config.spectrum.alpha_list);
    if (node.contains("Akp_list")) {
      config.spectrum.akp_list = node.at("Akp_list").get<std::vector<double>>();
    }
    if (node.contains("Alpha_list")) {
      config.spectrum.alpha_list = node.at("Alpha_list").get<std::vector<double>>();
    }
    if (node.contains("phases_deg")) {
      config.spectrum.phases_deg = node.at("phases_deg").get<std::vector<double>>();
    }
    config.spectrum.heading_deg = read_with_default(node, "heading_deg", config.spectrum.heading_deg);
    config.spectrum.spread_deg = read_with_default(node, "spread_deg", config.spectrum.spread_deg);
    config.spectrum.energy_keep_frac = read_with_default(node, "energy_keep_frac", config.spectrum.energy_keep_frac);
    config.spectrum.max_components = read_with_default(node, "max_components", config.spectrum.max_components);
  }

  if (root.contains("domain")) {
    const auto& node = root.at("domain");
    config.domain.lx_lambda = read_with_default(node, "Lx_lambda", config.domain.lx_lambda);
    config.domain.ly_lambda = read_with_default(node, "Ly_lambda", config.domain.ly_lambda);
    config.domain.nx = read_with_default(node, "Nx", config.domain.nx);
    config.domain.ny = read_with_default(node, "Ny", config.domain.ny);
  }

  if (root.contains("timing")) {
    const auto& node = root.at("timing");
    config.timing.t_focus = read_with_default(node, "t_focus", config.timing.t_focus);
    config.timing.t_init_periods_list = read_vector_or_scalar<double>(node, "t_init_periods_list", config.timing.t_init_periods_list);
    if (node.contains("t_init_periods")) {
      config.timing.t_init_periods_list = read_vector_or_scalar<double>(node, "t_init_periods", config.timing.t_init_periods_list);
    }
    config.timing.t_end_periods = read_with_default(node, "t_end_periods", config.timing.t_end_periods);
    config.timing.steps_per_period = read_with_default(node, "steps_per_period", config.timing.steps_per_period);
    config.timing.focus_edge_padding_fraction = read_with_default(
        node, "focus_edge_padding_fraction", config.timing.focus_edge_padding_fraction);
  }

  if (root.contains("output")) {
    const auto& node = root.at("output");
    config.output.output_root = std::filesystem::path(read_with_default(node, "output_root", config.output.output_root.string()));
    config.output.store_surface_stride = read_with_default(node, "store_surface_stride", config.output.store_surface_stride);
    config.output.surface_format = read_with_default(node, "surface_format", config.output.surface_format);
  }

  if (root.contains("visualization")) {
    const auto& node = root.at("visualization");
    config.visualization.enabled = read_with_default(node, "enabled", config.visualization.enabled);
    config.visualization.write_linear_field = read_with_default(
        node, "write_linear_field", config.visualization.write_linear_field);
    config.visualization.write_nonlinear_field = read_with_default(
        node, "write_nonlinear_field", config.visualization.write_nonlinear_field);
    config.visualization.write_surface_style = read_with_default(
        node, "write_surface_style", config.visualization.write_surface_style);
    config.visualization.centerline_section = read_with_default(
        node, "centerline_section", config.visualization.centerline_section);
    config.visualization.format = read_with_default(node, "format", config.visualization.format);
    config.visualization.width = read_with_default(node, "width", config.visualization.width);
    config.visualization.height = read_with_default(node, "height", config.visualization.height);
  }

  if (root.contains("phases_deg")) {
    config.spectrum.phases_deg = root.at("phases_deg").get<std::vector<double>>();
  }

  return config;
}

std::string config_summary_json(const GeneratorConfig& config) {
  json root;
  root["physics"] = {
      {"g", config.physics.g},
      {"kp", config.physics.kp},
      {"Tp", config.physics.tp},
      {"kd_list", config.physics.kd_list},
  };
  root["nonlinear"] = {
      {"order", config.nonlinear.order},
      {"subharmonic_mode", config.nonlinear.subharmonic_mode},
  };
  root["spectrum"] = {
      {"preset", config.spectrum.preset},
      {"Akp_list", config.spectrum.akp_list},
      {"Alpha_list", config.spectrum.alpha_list},
      {"phases_deg", config.spectrum.phases_deg},
      {"heading_deg", config.spectrum.heading_deg},
      {"spread_deg", config.spectrum.spread_deg},
      {"energy_keep_frac", config.spectrum.energy_keep_frac},
      {"max_components", config.spectrum.max_components},
  };
  root["domain"] = {
      {"Lx_lambda", config.domain.lx_lambda},
      {"Ly_lambda", config.domain.ly_lambda},
      {"Nx", config.domain.nx},
      {"Ny", config.domain.ny},
  };
  root["timing"] = {
      {"t_focus", config.timing.t_focus},
      {"t_init_periods_list", config.timing.t_init_periods_list},
      {"t_end_periods", config.timing.t_end_periods},
      {"steps_per_period", config.timing.steps_per_period},
      {"focus_edge_padding_fraction", config.timing.focus_edge_padding_fraction},
  };
  root["output"] = {
      {"output_root", config.output.output_root.string()},
      {"store_surface_stride", config.output.store_surface_stride},
      {"surface_format", config.output.surface_format},
  };
  root["visualization"] = {
      {"enabled", config.visualization.enabled},
      {"write_linear_field", config.visualization.write_linear_field},
      {"write_nonlinear_field", config.visualization.write_nonlinear_field},
      {"write_surface_style", config.visualization.write_surface_style},
      {"centerline_section", config.visualization.centerline_section},
      {"format", config.visualization.format},
      {"width", config.visualization.width},
      {"height", config.visualization.height},
  };
  return root.dump(2);
}

}  // namespace ow3d_directional
