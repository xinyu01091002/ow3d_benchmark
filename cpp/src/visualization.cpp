#include "ow3d_directional/visualization.hpp"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "ow3d_directional/ow3d_export.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace ow3d_directional {
namespace {

struct Image {
  int width = 0;
  int height = 0;
  std::vector<std::uint8_t> pixels;

  Image(int w, int h)
      : width(w),
        height(h),
        pixels(static_cast<std::size_t>(w) * static_cast<std::size_t>(h) * 3U, 255) {}

  void set_pixel(int x, int y, std::uint8_t r, std::uint8_t g, std::uint8_t b) {
    if (x < 0 || y < 0 || x >= width || y >= height) {
      return;
    }
    const std::size_t idx = (static_cast<std::size_t>(y) * static_cast<std::size_t>(width) +
                             static_cast<std::size_t>(x)) * 3U;
    pixels[idx + 0] = r;
    pixels[idx + 1] = g;
    pixels[idx + 2] = b;
  }
};

double clamp01(double value) {
  return std::min(1.0, std::max(0.0, value));
}

void save_png(const std::filesystem::path& path, const Image& image) {
  ensure_directory(path.parent_path());
  if (stbi_write_png(path.string().c_str(), image.width, image.height, 3, image.pixels.data(), image.width * 3) == 0) {
    throw std::runtime_error("Failed to write PNG: " + path.string());
  }
}

void draw_line(Image& image, int x0, int y0, int x1, int y1, std::uint8_t r, std::uint8_t g, std::uint8_t b) {
  const int dx = std::abs(x1 - x0);
  const int sx = x0 < x1 ? 1 : -1;
  const int dy = -std::abs(y1 - y0);
  const int sy = y0 < y1 ? 1 : -1;
  int err = dx + dy;

  while (true) {
    image.set_pixel(x0, y0, r, g, b);
    if (x0 == x1 && y0 == y1) {
      break;
    }
    const int e2 = 2 * err;
    if (e2 >= dy) {
      err += dy;
      x0 += sx;
    }
    if (e2 <= dx) {
      err += dx;
      y0 += sy;
    }
  }
}

void draw_rect_outline(Image& image, int x0, int y0, int x1, int y1, std::uint8_t r, std::uint8_t g, std::uint8_t b) {
  draw_line(image, x0, y0, x1, y0, r, g, b);
  draw_line(image, x1, y0, x1, y1, r, g, b);
  draw_line(image, x1, y1, x0, y1, r, g, b);
  draw_line(image, x0, y1, x0, y0, r, g, b);
}

void fill_rect(Image& image, int x0, int y0, int x1, int y1, std::uint8_t r, std::uint8_t g, std::uint8_t b) {
  const int xmin = std::max(0, std::min(x0, x1));
  const int xmax = std::min(image.width - 1, std::max(x0, x1));
  const int ymin = std::max(0, std::min(y0, y1));
  const int ymax = std::min(image.height - 1, std::max(y0, y1));
  for (int y = ymin; y <= ymax; ++y) {
    for (int x = xmin; x <= xmax; ++x) {
      image.set_pixel(x, y, r, g, b);
    }
  }
}

void colorize(double normalized, std::uint8_t* r, std::uint8_t* g, std::uint8_t* b) {
  const double t = clamp01(normalized);
  const double blue = clamp01(1.5 - 2.0 * t);
  const double red = clamp01(2.0 * t - 0.5);
  const double green = clamp01(1.0 - std::abs(2.0 * t - 1.0));
  *r = static_cast<std::uint8_t>(std::round(red * 255.0));
  *g = static_cast<std::uint8_t>(std::round(green * 255.0));
  *b = static_cast<std::uint8_t>(std::round(blue * 255.0));
}

double matrix_min(const mf12_cpp::Matrix& matrix) {
  return *std::min_element(matrix.values.begin(), matrix.values.end());
}

double matrix_max(const mf12_cpp::Matrix& matrix) {
  return *std::max_element(matrix.values.begin(), matrix.values.end());
}

int resolve_centerline_index(const GeneratorConfig& config, const mf12_cpp::Matrix& y) {
  if (config.visualization.centerline_section != "mid") {
    return static_cast<int>(y.values.size() / 2U);
  }
  int best_index = 0;
  double best_abs = std::abs(y.values.front());
  for (std::size_t i = 1; i < y.values.size(); ++i) {
    const double candidate = std::abs(y.values[i]);
    if (candidate < best_abs) {
      best_abs = candidate;
      best_index = static_cast<int>(i);
    }
  }
  return best_index;
}

Image make_heatmap(const mf12_cpp::Matrix& field, int width, int height) {
  Image image(width, height);
  const double min_v = matrix_min(field);
  const double max_v = matrix_max(field);
  const double span = std::max(max_v - min_v, 1e-12);

  for (int py = 0; py < height; ++py) {
    const double v = static_cast<double>(py) / static_cast<double>(std::max(height - 1, 1));
    const std::size_t src_y = static_cast<std::size_t>(std::llround((1.0 - v) * static_cast<double>(field.rows - 1)));
    for (int px = 0; px < width; ++px) {
      const double u = static_cast<double>(px) / static_cast<double>(std::max(width - 1, 1));
      const std::size_t src_x = static_cast<std::size_t>(std::llround(u * static_cast<double>(field.cols - 1)));
      const double value = field.values[src_y * field.cols + src_x];
      std::uint8_t r = 0;
      std::uint8_t g = 0;
      std::uint8_t b = 0;
      colorize((value - min_v) / span, &r, &g, &b);
      image.set_pixel(px, py, r, g, b);
    }
  }

  draw_rect_outline(image, 0, 0, width - 1, height - 1, 20, 20, 20);
  return image;
}

Image make_overview_plot(
    const mf12_cpp::Matrix& eta_linear,
    const mf12_cpp::Matrix& eta_nonlinear,
    const mf12_cpp::Matrix& phi_linear,
    const mf12_cpp::Matrix& phi_nonlinear,
    const GeneratorConfig& config) {
  Image image(config.visualization.width, config.visualization.height);
  fill_rect(image, 0, 0, image.width - 1, image.height - 1, 255, 255, 255);

  const int margin = 40;
  const int gap = 20;
  const int panel_w = (image.width - 2 * margin - gap) / 2;
  const int panel_h = (image.height - 2 * margin - gap) / 2;

  auto blit = [&](const Image& src, int x0, int y0) {
    for (int y = 0; y < src.height; ++y) {
      for (int x = 0; x < src.width; ++x) {
        const std::size_t idx = (static_cast<std::size_t>(y) * static_cast<std::size_t>(src.width) +
                                 static_cast<std::size_t>(x)) * 3U;
        image.set_pixel(x0 + x, y0 + y, src.pixels[idx], src.pixels[idx + 1], src.pixels[idx + 2]);
      }
    }
  };

  blit(make_heatmap(eta_linear, panel_w, panel_h), margin, margin);
  blit(make_heatmap(eta_nonlinear, panel_w, panel_h), margin + panel_w + gap, margin);

  auto draw_series_panel = [&](int x0, int y0, const mf12_cpp::Matrix& linear, const mf12_cpp::Matrix& nonlinear) {
    fill_rect(image, x0, y0, x0 + panel_w, y0 + panel_h, 250, 250, 250);
    draw_rect_outline(image, x0, y0, x0 + panel_w, y0 + panel_h, 30, 30, 30);
    const int line_index = static_cast<int>(linear.rows / 2U);
    double min_v = std::numeric_limits<double>::max();
    double max_v = -std::numeric_limits<double>::max();
    for (std::size_t i = 0; i < linear.cols; ++i) {
      min_v = std::min(min_v, linear.values[static_cast<std::size_t>(line_index) * linear.cols + i]);
      max_v = std::max(max_v, linear.values[static_cast<std::size_t>(line_index) * linear.cols + i]);
      min_v = std::min(min_v, nonlinear.values[static_cast<std::size_t>(line_index) * nonlinear.cols + i]);
      max_v = std::max(max_v, nonlinear.values[static_cast<std::size_t>(line_index) * nonlinear.cols + i]);
    }
    const double span = std::max(max_v - min_v, 1e-12);
    auto ypix = [&](double value) {
      const double t = (value - min_v) / span;
      return y0 + panel_h - 1 - static_cast<int>(std::llround(t * static_cast<double>(panel_h - 1)));
    };

    for (std::size_t i = 1; i < linear.cols; ++i) {
      const int px0 = x0 + static_cast<int>(std::llround((static_cast<double>(i - 1) / static_cast<double>(linear.cols - 1)) * (panel_w - 1)));
      const int px1 = x0 + static_cast<int>(std::llround((static_cast<double>(i) / static_cast<double>(linear.cols - 1)) * (panel_w - 1)));
      const double l0 = linear.values[static_cast<std::size_t>(line_index) * linear.cols + (i - 1)];
      const double l1 = linear.values[static_cast<std::size_t>(line_index) * linear.cols + i];
      const double n0 = nonlinear.values[static_cast<std::size_t>(line_index) * nonlinear.cols + (i - 1)];
      const double n1 = nonlinear.values[static_cast<std::size_t>(line_index) * nonlinear.cols + i];
      draw_line(image, px0, ypix(l0), px1, ypix(l1), 25, 90, 170);
      draw_line(image, px0, ypix(n0), px1, ypix(n1), 200, 50, 40);
    }
  };

  draw_series_panel(margin, margin + panel_h + gap, eta_linear, eta_nonlinear);
  draw_series_panel(margin + panel_w + gap, margin + panel_h + gap, phi_linear, phi_nonlinear);
  return image;
}

Image make_surface_style_plot(const mf12_cpp::Matrix& field, int width, int height) {
  Image image(width, height);
  fill_rect(image, 0, 0, width - 1, height - 1, 255, 255, 255);
  const double min_v = matrix_min(field);
  const double max_v = matrix_max(field);
  const double span = std::max(max_v - min_v, 1e-12);
  const int margin = 40;
  const int usable_w = width - 2 * margin;
  const int usable_h = height - 2 * margin;

  for (std::size_t row = 0; row < field.rows; row += std::max<std::size_t>(1, field.rows / 36U)) {
    int prev_x = 0;
    int prev_y = 0;
    bool has_prev = false;
    for (std::size_t col = 0; col < field.cols; ++col) {
      const double z = field.values[row * field.cols + col];
      const double zn = (z - min_v) / span - 0.5;
      const double x = static_cast<double>(col) / static_cast<double>(std::max<std::size_t>(field.cols - 1, 1));
      const double y = static_cast<double>(row) / static_cast<double>(std::max<std::size_t>(field.rows - 1, 1));
      const int px = margin + static_cast<int>(std::llround((x * 0.8 + y * 0.2) * usable_w));
      const int py = margin + usable_h - static_cast<int>(std::llround((y * 0.55 + (zn + 0.5) * 0.35) * usable_h));
      if (has_prev) {
        draw_line(image, prev_x, prev_y, px, py, 50, 80, 140);
      }
      prev_x = px;
      prev_y = py;
      has_prev = true;
    }
  }

  draw_rect_outline(image, 0, 0, width - 1, height - 1, 20, 20, 20);
  return image;
}

}  // namespace

void write_visualizations(
    const std::filesystem::path& output_dir,
    const GeneratorConfig& config,
    const CaseParameters& params,
    const Mf12RunOutputs& outputs) {
  (void) params;

  if (config.visualization.write_linear_field) {
    save_png(output_dir / "linear_field_eta.png",
             make_heatmap(outputs.linear.eta, config.visualization.width, config.visualization.height));
  }
  if (config.visualization.write_nonlinear_field) {
    save_png(output_dir / "nonlinear_field_eta.png",
             make_heatmap(outputs.nonlinear.eta, config.visualization.width, config.visualization.height));
  }

  save_png(output_dir / "mf12_directional_overview.png",
           make_overview_plot(outputs.linear.eta, outputs.nonlinear.eta, outputs.linear.phi, outputs.nonlinear.phi, config));

  if (config.visualization.write_surface_style) {
    save_png(output_dir / "linear_wave_group_surface.png",
             make_surface_style_plot(outputs.linear.eta, 900, 700));
  }
}

}  // namespace ow3d_directional
