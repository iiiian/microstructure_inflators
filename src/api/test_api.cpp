#include <inflator/Inflator.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
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
		require(mesh.dimension == 2 || mesh.dimension == 3, "Wrong dimension");
		require(mesh.vertices.size() % mesh.dimension == 0, "Wrong coordinate layout");
		require(mesh.elements.size() % (mesh.dimension + 1) == 0, "Wrong element layout");
		const std::size_t vertex_count = mesh.vertices.size() / mesh.dimension;
		require(mesh.parameter_kinds.size() == mesh.shape_parameters.size(), "Missing parameter kinds");
		require(mesh.boundary_shape_velocities.size() == mesh.shape_parameters.size() * mesh.vertices.size(), "Wrong velocity layout");
		for (const double velocity : mesh.boundary_shape_velocities)
		{
			require(std::isfinite(velocity), "Nonfinite velocity");
		}
		const auto coordinate = [&](std::size_t vertex, int axis) { return mesh.vertices[vertex * mesh.dimension + axis]; };
		for (std::size_t offset = 0; offset < mesh.elements.size(); offset += mesh.dimension + 1)
		{
			const auto *element = mesh.elements.data() + offset;
			for (int corner = 0; corner <= mesh.dimension; ++corner)
			{
				require(element[corner] < vertex_count, "Invalid element index");
			}
			std::array<std::array<double, 3>, 3> edges{};
			for (int corner = 0; corner < mesh.dimension; ++corner)
			{
				for (int axis = 0; axis < mesh.dimension; ++axis)
				{
					edges[corner][axis] = coordinate(element[corner + 1], axis) - coordinate(element[0], axis);
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
			double minimum = coordinate(0, axis);
			double maximum = minimum;
			for (std::size_t vertex = 0; vertex < vertex_count; ++vertex)
			{
				minimum = std::min(minimum, coordinate(vertex, axis));
				maximum = std::max(maximum, coordinate(vertex, axis));
			}
			for (std::size_t vertex = 0; vertex < vertex_count; ++vertex)
			{
				if (std::abs(coordinate(vertex, axis) - minimum) > 1e-10)
				{
					continue;
				}
				bool found = false;
				for (std::size_t candidate = 0; candidate < vertex_count; ++candidate)
				{
					double distance = 0;
					for (int component = 0; component < mesh.dimension; ++component)
					{
						distance += std::abs(coordinate(candidate, component) - coordinate(vertex, component)
											 - (component == axis ? maximum - minimum : 0));
					}
					found |= distance < 1e-9;
				}
				require(found, "Unmatched periodic vertex");
			}
		}
	}

	template <typename Callable>
	void require_invalid(Callable call)
	{
		bool caught = false;
		try
		{
			call();
		}
		catch (const std::invalid_argument &)
		{
			caught = true;
		}
		require(caught, "Expected invalid_argument");
	}
} // namespace

int main(int argc, char **argv)
{
	const auto directory = std::filesystem::temp_directory_path()
						   / ("inflator-api-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
	std::filesystem::create_directory(directory);
	try
	{
		const std::array<inflator::CellType, 3> cell_types = {inflator::CellType::Periodic2D, inflator::CellType::Orthotropic2D, inflator::CellType::Orthotropic3D};
		const int first = argc > 1 ? std::stoi(argv[1]) : 0;
		const int end = argc > 1 ? first + 1 : cell_types.size();
		require(first >= 0 && end <= int(cell_types.size()), "Invalid test mode");
		for (int mode = first; mode < end; ++mode)
		{
			const auto cell_type = cell_types[mode];
			const int dimension = cell_type == inflator::CellType::Orthotropic3D ? 3 : 2;
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
			inflator::Options options;
			options.cell_type = cell_type;
			options.default_thickness = 0.15;
			options.max_triangle_area = 0.02;
			options.contour_resolution = 24;
			options.force_contour_resolution = true;
			options.facet_size = 0.2;
			options.facet_distance = 0.02;
			options.cell_size = 0.4;
			options.edge_size = 0.2;
			const auto cell = inflator::inflate(wire, options);
			require(cell.dimension == dimension, "Wrong dimension");
			check(cell);
			auto parameters = cell.shape_parameters;
			const auto repeated = inflator::inflate(wire, options, parameters.data(), parameters.size());
			require(parameters == cell.shape_parameters, "Borrowed input was modified");
			require(cell.vertices == repeated.vertices && cell.elements == repeated.elements, "Nondeterministic meshing");
			require(cell.boundary_shape_velocities == repeated.boundary_shape_velocities, "Nondeterministic velocities");
			options.tiles = 2;
			const auto tiled = inflator::inflate(wire, options, parameters.data(), parameters.size());
			check(tiled);
			require(tiled.elements.size() == cell.elements.size() * (1 << dimension), "Incorrect tile count");
			require(tiled.vertices.size() < cell.vertices.size() * (1 << dimension), "Tile seams were not stitched");
			std::cout << dimension << "D: " << cell.vertices.size() / dimension << " vertices, " << cell.elements.size() / (dimension + 1) << " elements\n";

			require_invalid([&]() { (void)inflator::inflate(wire, options, nullptr, parameters.size()); });
			require_invalid([&]() { (void)inflator::inflate(wire, options, parameters.data(), parameters.size() - 1); });
			parameters[0] = std::numeric_limits<double>::quiet_NaN();
			require_invalid([&]() { (void)inflator::inflate(wire, options, parameters.data(), parameters.size()); });
			options.tiles = 0;
			require_invalid([&]() { (void)inflator::inflate(wire, options); });
			options.tiles = 1;
			options.default_thickness = -1;
			require_invalid([&]() { (void)inflator::inflate(wire, options); });
			options.default_thickness = 0.15;
			options.max_triangle_area = 0;
			options.cell_size = 0;
			require_invalid([&]() { (void)inflator::inflate(wire, options); });
			options.cell_type = static_cast<inflator::CellType>(-1);
			require_invalid([&]() { (void)inflator::inflate(wire, options); });
		}
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
