#include "geometry.h"

#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>

namespace primitives {
namespace geometry {

const float kWidthRange = 1;
const float kLengthRange = 1;

Eigen::Vector2d GeneratePoint() {
  Eigen::Vector2d p(
      static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / kWidthRange)),
      static_cast<float>(rand()) /
          (static_cast<float>(RAND_MAX / kLengthRange)));
  return p;
}

std::vector<float> GenerateRandomVec(const float size, const float min_range,
                                     const float max_range) {
  std::random_device random_device;
  std::mt19937 mersenne_engine{random_device()};
  std::uniform_real_distribution<float> distribution{min_range, max_range};

  auto gen = [&distribution, &mersenne_engine]() {
    return distribution(mersenne_engine);
  };
  std::vector<float> vec(size);
  std::generate(std::begin(vec), std::end(vec), gen);
  return vec;
}

std::vector<double> ComputeDirectionalComponents(const size_t num_sides,
                                                 const size_t range) {
  std::vector<float> coordinates = GenerateRandomVec(num_sides, 0, range);
  std::sort(coordinates.begin(), coordinates.end());

  // Generate binary mask
  std::vector<uint32_t> binary_mask(num_sides - 2, 0);
  size_t split_size = std::ceil(binary_mask.size() / 2);
  for (size_t i = 0; i < split_size; ++i) {
    binary_mask[i] = 1;
  }

  // Shuffle the array.
  std::shuffle(
      binary_mask.begin(), binary_mask.end(),
      std::default_random_engine(
          std::chrono::system_clock::now().time_since_epoch().count()));

  std::vector<double> directional_vectors;
  directional_vectors.reserve(num_sides);
  size_t previous_index_up = 0;
  size_t previous_index_down = 0;
  double sum_up = 0;
  double sum_down = 0;
  for (size_t i = 0; i < binary_mask.size(); ++i) {
    if (binary_mask[i]) {
      directional_vectors.push_back(
          -1 * (coordinates[i + 1] - coordinates[previous_index_up]));
      sum_up += std::abs((coordinates[i + 1] - coordinates[previous_index_up]));
      previous_index_up = i + 1;
    } else {
      directional_vectors.push_back(coordinates[i + 1] -
                                    coordinates[previous_index_down]);
      sum_down += (coordinates[i + 1] - coordinates[previous_index_down]);
      previous_index_down = i + 1;
    }
  }
  // Add two vectors connecting to the last coordinate, one going up and the
  // other down.
  size_t index = coordinates.size() - 1;
  directional_vectors.push_back(coordinates[index] -
                                coordinates[previous_index_down]);
  directional_vectors.push_back(
      -1 * (coordinates[index] - coordinates[previous_index_up]));
  sum_up += (coordinates[index] - coordinates[previous_index_up]);
  sum_down += (coordinates[index] - coordinates[previous_index_down]);

  float coord_range = coordinates[index] - coordinates[0];
  float kErr = 0.00001;
  assert(std::abs(sum_up - coord_range) < kErr);
  assert(std::abs(sum_down - coord_range) < kErr);
  assert(directional_vectors.size() == num_sides);
  return directional_vectors;
}

Polygon ConstructConvexPolygon(const size_t num_sides) {
  std::vector<double> x_components =
      ComputeDirectionalComponents(num_sides, kWidthRange);
  std::vector<double> y_components =
      ComputeDirectionalComponents(num_sides, kLengthRange);
  assert(x_components.size() == num_sides);
  assert(y_components.size() == num_sides);
  int current_int = 0;
  auto increment_int = [&current_int]() { return current_int++; };
  std::vector<size_t> x_sequence(x_components.size());
  std::generate_n(x_sequence.begin(), x_sequence.size(), increment_int);
  std::vector<size_t> y_sequence = x_sequence;
  std::shuffle(
      x_sequence.begin(), x_sequence.end(),
      std::default_random_engine(
          std::chrono::system_clock::now().time_since_epoch().count()));
  std::shuffle(
      y_sequence.begin(), y_sequence.end(),
      std::default_random_engine(
          std::chrono::system_clock::now().time_since_epoch().count()));
  assert(x_sequence.size() == num_sides);
  assert(x_components.size() == num_sides);
  // Randomly pair each directional component and compute vectors.
  std::vector<Eigen::Vector2d> polygon_sides;
  polygon_sides.reserve(x_sequence.size());
  for (size_t i = 0; i < num_sides; ++i) {
    polygon_sides.emplace_back(x_components[x_sequence[i]],
                               y_components[y_sequence[i]]);
  }
  // Compute the angle of each vector with respect to the +y-axis.
  std::vector<std::pair<size_t, float>> index_angles;
  for (size_t i = 0; i < polygon_sides.size(); ++i) {
    float angle = std::atan2(polygon_sides[i].y(), polygon_sides[i].x());
    if (angle < 0) {
      angle += 2 * M_PI;
    }
    index_angles.push_back(std::make_pair(i, angle));
  }
  auto compare_angles = [](const std::pair<size_t, float>& a,
                           const std::pair<size_t, float>& b) {
    return a.second < b.second;
  };
  std::sort(index_angles.begin(), index_angles.end(), compare_angles);

  Polygon polygon(num_sides);
  Eigen::Vector2d current_start = Eigen::Vector2d::Zero();
  for (size_t i = 0; i < index_angles.size(); ++i) {
    polygon[i].first = current_start;
    polygon[i].second = current_start + polygon_sides[index_angles[i].first];
    current_start = polygon[i].second;
    std::cout << polygon[i].first.transpose() << " --> "
              << polygon[i].second.transpose() << std::endl;
  }
  return polygon;
}

}  // namespace geometry
}  // namespace primitives
