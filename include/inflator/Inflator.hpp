#pragma once

#include <cstddef>
#include <filesystem>
#include <vector>

#define INFLATOR_EXPORT __attribute__((visibility("default")))

namespace inflator
{
	/// Supported geometry and symmetry combinations.
	enum class CellType
	{
		Periodic2D,
		Orthotropic2D,
		Orthotropic3D
	};

	/// Meaning of a scalar shape parameter.
	enum class ParameterKind
	{
		Position,
		Thickness,
		Blending
	};

	struct Options
	{
		CellType cell_type = CellType::Periodic2D;
		/// Number of stitched repetitions along each active axis.
		std::size_t tiles = 1;
		/// Graph-edge neighborhood used to evaluate the implicit geometry.
		std::size_t graph_radius = 2;
		/// Initial bar radius when shape parameters are omitted.
		double default_thickness = 0.07;

		/// Sampling resolution for 2D contours and 3D boundary features.
		std::size_t contour_resolution = 256;
		/// In 2D, use contour_resolution instead of deriving it from triangle area.
		bool force_contour_resolution = false;
		/// Maximum triangle area for 2D meshing.
		double max_triangle_area = 0.001;

		/// Surface-facet size bound for 3D meshing.
		double facet_size = 0.025;
		/// Surface-approximation distance bound for 3D meshing.
		double facet_distance = 0.002;
		/// Tetrahedral cell-size bound for 3D meshing.
		double cell_size = 0.15;
		/// Feature-edge size bound for 3D meshing.
		double edge_size = 0.025;
	};

	/// Owns generated geometry, evaluated parameters, and boundary shape velocities.
	struct Mesh
	{
		/// Spatial dimension: 2 for triangles, 3 for tetrahedra.
		int dimension = 0;
		/// Vertex-major coordinates; dimension entries per vertex.
		std::vector<double> vertices;
		/// Element-major, zero-based vertex indices; dimension + 1 per element.
		std::vector<std::size_t> elements;
		/// Evaluated values, ordered as position DOFs, radii, then blending values.
		std::vector<double> shape_parameters;
		/// One classification per entry in shape_parameters.
		std::vector<ParameterKind> parameter_kinds;
		/// Layout [parameter][vertex][coordinate]; normal boundary velocities, zero in the interior.
		/// These describe geometric shape variation, not discrete remeshing.
		std::vector<double> boundary_shape_velocities;
	};

	/**
	 * @brief Inflate an OBJ line graph into a stitched periodic volume mesh.
	 * Coordinates and meshing lengths use normalized cell units, without application-specific scaling.
	 * @param wire_path Graph file used to construct the symmetry-reduced geometry.
	 * @param options Geometry, tiling, and meshing settings.
	 * @param shape_parameters Borrowed values in graph parameter order. Omit to use graph-derived
	 * positions, options.default_thickness, and default joint blending.
	 * @param parameter_count Number of values; zero selects defaults, otherwise must match the graph.
	 * @return Owning mesh, actual parameters, and corresponding boundary shape velocities.
	 * @throws std::invalid_argument Invalid settings or parameter input.
	 * @throws std::runtime_error Graph loading or mesh generation failed.
	 * Inputs are not retained. Calls are synchronous and internally serialized.
	 */
	[[nodiscard]] INFLATOR_EXPORT Mesh inflate(
		const std::filesystem::path &wire_path,
		const Options &options,
		const double *shape_parameters = nullptr,
		std::size_t parameter_count = 0);
} // namespace inflator
