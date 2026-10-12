#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#if defined(_WIN32)
#if defined(INFLATOR_BUILDING_LIBRARY)
#define INFLATOR_EXPORT __declspec(dllexport)
#else
#define INFLATOR_EXPORT __declspec(dllimport)
#endif
#else
#define INFLATOR_EXPORT __attribute__((visibility("default")))
#endif

namespace inflator
{
    enum class ParameterType { Position, Thickness, Blending, Other };

    struct Request
    {
        std::string type = "2D_doubly_periodic";
        std::string wire_path;
        std::vector<double> parameters;
        std::string meshing_options = "{}";
        int tiles = 1;
        std::size_t graph_radius = 2;
        double default_thickness = 0.07;
    };

    struct Mesh
    {
        int dimension;
        std::vector<double> parameters;
        std::vector<ParameterType> parameter_types;
        std::vector<std::array<double, 3>> vertices;
        std::vector<std::array<std::size_t, 4>> elements;
        std::vector<std::vector<std::array<double, 3>>> shape_velocities;
    };

    /// Generate a stitched periodic mesh and parameter-major boundary normal velocities.
    /// Empty parameters select defaults. Triangles use the first three element indices.
    /// Throws std::exception on failure; no third-party headers or types cross this interface.
    INFLATOR_EXPORT Mesh inflate(const Request &request);
}
