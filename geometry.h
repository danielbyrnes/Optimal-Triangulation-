#pragma once

#include <Eigen/Dense>
#include <utility>
#include <vector>

namespace primitives {
namespace geometry {

struct triplet {
  size_t i;
  size_t j;
  size_t k;
};
struct triangle {
  Eigen::Vector2d a;
  Eigen::Vector2d b;
  Eigen::Vector2d c;
};
typedef std::pair<Eigen::Vector2d, Eigen::Vector2d> LineSegment;
typedef std::vector<LineSegment> Polygon;

Eigen::Vector2d GeneratePoint();
Polygon ConstructConvexPolygon(const size_t num_sides);

}  // namespace geometry
}  // namespace primitives
