#include "geometry.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <limits>

// Use recursion to trace the minimum weight triangulation.
void TraceOptimalTriangleIndices(
    const std::vector<std::vector<double>>& triangulation_cost, size_t i,
    size_t k, std::vector<primitives::geometry::triplet>* triangulation) {
  if (i + 1 == k) {
    return;
  }
  double min_cost = std::numeric_limits<double>::max();
  size_t min_j = 0;
  // compute j that minimizes cost.
  for (size_t j = i + 1; j < k; ++j) {
    if (triangulation_cost[i][j] + triangulation_cost[j][k] < min_cost) {
      min_cost = triangulation_cost[i][j] + triangulation_cost[j][k];
      min_j = j;
    }
  }
  primitives::geometry::triplet triangle_indices;
  triangle_indices.i = i;
  triangle_indices.j = min_j;
  triangle_indices.k = k;
  triangulation->push_back(triangle_indices);

  TraceOptimalTriangleIndices(triangulation_cost, i, min_j, triangulation);
  TraceOptimalTriangleIndices(triangulation_cost, min_j, k, triangulation);
}

// Use dynamic programming to compute the minimal cost triangulation.
std::vector<primitives::geometry::triplet> TriangulatePolygon(
    const primitives::geometry::Polygon& polygon) {
  size_t num_sides = polygon.size();
  std::vector<primitives::geometry::triplet> triangulation;
  std::vector<std::vector<double>> triangulation_cost(num_sides);
  for (int i = 0; i < num_sides; ++i) {
    triangulation_cost[i].resize(num_sides, std::numeric_limits<double>::max());
  }
  for (int d = 1; d < num_sides; ++d) {
    for (int i = 0; i < num_sides - d; ++i) {
      int k = i + d;
      if (d == 1) {
        triangulation_cost[i][k] = 0;
      }
      for (int j = i + 1; j < k; ++j) {
        // compute distance from i to j and j to k.
        std::cout << triangulation_cost[i][j] << " " << triangulation_cost[j][k]
                  << std::endl;
        std::cout << "i->j: " << (polygon[i].first - polygon[j].first).norm()
                  << std::endl;
        std::cout << "j->k: " << (polygon[k].first - polygon[j].first).norm()
                  << std::endl;
        double cost = triangulation_cost[i][j] + triangulation_cost[j][k] +
                      (polygon[i].first - polygon[j].first).norm() +
                      (polygon[k].first - polygon[j].first).norm();
        if (cost < triangulation_cost[i][k]) {
          triangulation_cost[i][k] = cost;
        }
      }
      std::cout << std::endl;
    }
  }

  // TODO(dcbyrnes) Remove this once done debugging.
  for (size_t i = 0; i < num_sides; ++i) {
    for (size_t j = 0; j < num_sides; ++j) {
      if (j < i + 1) {
        std::cout << " ";
      } else {
        std::cout << " " << triangulation_cost[i][j] << " ";
      }
    }
    std::cout << std::endl;
  }

  std::cout << "Triangulation cost: " << triangulation_cost[0][num_sides - 1]
            << std::endl;
  // Compute the triangle indices of the minimum weight triangulation.
  TraceOptimalTriangleIndices(triangulation_cost, 0, num_sides - 1,
                              &triangulation);

  assert(triangulation.size() == num_sides - 2);
  return triangulation;
}

int main(int argc, char* argv[]) {
  int num_sides = 0;
  if (argc == 1 || atoi(argv[1]) < 3) {
    num_sides = 10;
  } else {
    num_sides = atoi(argv[1]);
  }
  printf("Constructing convex polygon with %d sides...\n", num_sides);
  auto polygon = primitives::geometry::ConstructConvexPolygon(num_sides);

  const std::string kPolygonFilename("polygon.txt");
  printf("Writing polygon vertices to file:\t %s \n", kPolygonFilename.c_str());
  std::ofstream polygon_output_file;
  polygon_output_file.open(kPolygonFilename);
  for (auto const& line_segment : polygon) {
    polygon_output_file << line_segment.first.x() << ","
                        << line_segment.first.y() << std::endl;
  }
  polygon_output_file.close();

  printf("Computing the minimum weight triangulation...\n");
  std::vector<primitives::geometry::triplet> triangulation =
      TriangulatePolygon(polygon);

  const std::string kTriangulationFilename("triangle_indices.txt");
  printf("Writing triangulation to file:\t %s \n",
         kTriangulationFilename.c_str());
  std::ofstream triangulation_output_file;
  triangulation_output_file.open(kTriangulationFilename);
  for (auto const& triangle_indices : triangulation) {
    printf("a: (%f,%f) b: (%f,%f) c: (%f,%f) \n",
           polygon[triangle_indices.i].first.x(),
           polygon[triangle_indices.i].first.y(),
           polygon[triangle_indices.j].first.x(),
           polygon[triangle_indices.j].first.y(),
           polygon[triangle_indices.k].first.x(),
           polygon[triangle_indices.k].first.y());
    triangulation_output_file << triangle_indices.i << "," << triangle_indices.j
                              << "," << triangle_indices.k << std::endl;
  }
  triangulation_output_file.close();
  return EXIT_SUCCESS;
}
