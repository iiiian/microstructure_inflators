#include <inflator/Inflator.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace
{
    void require(bool condition, const char *message)
    {
        if (!condition)
        {
            throw std::runtime_error(message);
        }
    }

    void check(const inflator::Mesh &mesh)
    {
        require(!mesh.elements.empty(), "Empty mesh");
        require(mesh.shape_velocities.size() == mesh.parameters.size(), "Missing velocities");
        for (const auto &velocities : mesh.shape_velocities)
        {
            require(velocities.size() == mesh.vertices.size(), "Wrong velocity size");
            for (const auto &velocity : velocities)
            {
                for (const double component : velocity)
                {
                    require(std::isfinite(component), "Nonfinite velocity");
                }
            }
        }
        for (const auto &element : mesh.elements)
        {
            for (int corner = 0; corner <= mesh.dimension; ++corner)
            {
                require(element[corner] < mesh.vertices.size(), "Invalid element index");
            }
            std::array<std::array<double, 3>, 3> edges{};
            for (int corner = 0; corner < mesh.dimension; ++corner)
            {
                for (int axis = 0; axis < mesh.dimension; ++axis)
                {
                    edges[corner][axis] = mesh.vertices[element[corner + 1]][axis] - mesh.vertices[element[0]][axis];
                }
            }
            double determinant = edges[0][0] * edges[1][1] - edges[0][1] * edges[1][0];
            if (mesh.dimension == 3)
            {
                determinant = edges[0][0] * (edges[1][1] * edges[2][2] - edges[1][2] * edges[2][1])
                    - edges[0][1] * (edges[1][0] * edges[2][2] - edges[1][2] * edges[2][0])
                    + edges[0][2] * (edges[1][0] * edges[2][1] - edges[1][1] * edges[2][0]);
            }
            require(determinant > 0, "Inverted or degenerate element");
        }
        for (int axis = 0; axis < mesh.dimension; ++axis)
        {
            double minimum = mesh.vertices.front()[axis];
            double maximum = minimum;
            for (const auto &vertex : mesh.vertices)
            {
                minimum = std::min(minimum, vertex[axis]);
                maximum = std::max(maximum, vertex[axis]);
            }
            for (const auto &vertex : mesh.vertices)
            {
                if (std::abs(vertex[axis] - minimum) > 1e-10)
                {
                    continue;
                }
                bool found = false;
                for (const auto &candidate : mesh.vertices)
                {
                    double distance = 0;
                    for (int component = 0; component < mesh.dimension; ++component)
                    {
                        distance += std::abs(candidate[component] - vertex[component]
                            - (component == axis ? maximum - minimum : 0));
                    }
                    found |= distance < 1e-9;
                }
                require(found, "Unmatched periodic vertex");
            }
        }
    }
}

int main(int argc, char **argv)
{
    const auto directory = std::filesystem::temp_directory_path()
        / ("inflator-api-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::filesystem::create_directory(directory);
    try
    {
        for (const int dimension : {2, 3})
        {
            if (argc > 1 && dimension != std::stoi(argv[1]))
            {
                continue;
            }
            std::cout << "Testing " << dimension << "D inflation" << std::endl;
            const auto wire = directory / "wire.obj";
            std::ofstream output(wire);
            output << "v 0 0 0\nv 1 0 0\nv 0 1 0\nv -1 0 0\nv 0 -1 0\n";
            if (dimension == 3)
            {
                output << "v 0 0 1\nv 0 0 -1\n";
            }
            output << "l 1 2\nl 1 3\nl 1 4\nl 1 5\n";
            if (dimension == 3)
            {
                output << "l 1 6\nl 1 7\n";
            }
            output.close();
            inflator::Request request;
            request.type = dimension == 2 ? "2D_orthotropic" : "orthotropic";
            request.wire_path = wire.string();
            request.default_thickness = 0.15;
            request.meshing_options = R"({"maxArea":0.02,"marchingSquaresGridSize":24,"forceMSGridSize":true,"facetSize":0.2,"facetDistance":0.02,"cellSize":0.4,"edgeSize":0.2})";
            const auto cell = inflator::inflate(request);
            require(cell.dimension == dimension, "Wrong dimension");
            check(cell);
            request.parameters = cell.parameters;
            const auto repeated = inflator::inflate(request);
            require(cell.vertices == repeated.vertices && cell.elements == repeated.elements, "Nondeterministic meshing");
            request.tiles = 2;
            const auto tiled = inflator::inflate(request);
            check(tiled);
            require(tiled.elements.size() == cell.elements.size() * (1 << dimension), "Incorrect tile count");
            require(tiled.vertices.size() < cell.vertices.size() * (1 << dimension), "Tile seams were not stitched");
            std::cout << dimension << "D: " << cell.vertices.size() << " vertices, " << cell.elements.size() << " elements\n";
        }
        inflator::Request invalid;
        invalid.tiles = 0;
        bool caught = false;
        try
        {
            inflator::inflate(invalid);
        }
        catch (const std::invalid_argument &)
        {
            caught = true;
        }
        require(caught, "Exceptions do not cross the library boundary");
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        std::filesystem::remove_all(directory);
        return 1;
    }
    std::filesystem::remove_all(directory);
    return 0;
}
