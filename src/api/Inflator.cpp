#include <inflator/Inflator.hpp>

#include <IsosurfaceInflator.hh>
#include <MeshFEM/PeriodicBoundaryMatcher.hh>
#include <MeshFEM/SimplicialMesh.hh>

#include <cmath>
#include <array>
#include <mutex>
#include <numeric>
#include <stdexcept>

namespace inflator
{
	namespace
	{
		template <std::size_t Dimension>
		Mesh copy_and_tile(const IsosurfaceInflator &source, const Options &options)
		{
			Mesh result;
			result.dimension = Dimension;
			result.shape_parameters = source.inflatedParams();
			for (std::size_t parameter = 0; parameter < result.shape_parameters.size(); ++parameter)
			{
				if (source.isPositionParam(parameter))
				{
					result.parameter_kinds.push_back(ParameterKind::Position);
				}
				else if (source.isThicknessParam(parameter))
				{
					result.parameter_kinds.push_back(ParameterKind::Thickness);
				}
				else
				{
					result.parameter_kinds.push_back(ParameterKind::Blending);
				}
			}

			const auto &vertices = source.vertices();
			const BBox<PointND<Dimension>> bounds(vertices);
			const auto dimensions = bounds.dimensions();
			SimplicialMesh<Dimension> topology(source.elements(), vertices.size());
			std::vector<PointND<Dimension>> boundary_points;
			for (const auto &vertex : topology.boundaryVertices())
			{
				const auto &point = vertices[vertex.volumeVertex().index()].point;
				boundary_points.push_back(truncateFrom3D<PointND<Dimension>>(point));
			}

			std::vector<PeriodicBoundaryMatcher::FaceMembership<Dimension>> membership;
			std::vector<std::vector<std::size_t>> node_sets;
			std::vector<std::size_t> node_set_for_node;
			PeriodicBoundaryMatcher::determineCellBoundaryFaceMembership(boundary_points, bounds, membership);
			PeriodicBoundaryMatcher::match(boundary_points, bounds, membership, node_sets, node_set_for_node);

			std::vector<std::vector<std::size_t>> equivalence(node_sets.size());
			for (std::size_t group = 0; group < node_sets.size(); ++group)
			{
				for (const auto boundary_vertex : node_sets[group])
				{
					equivalence[group].push_back(topology.boundaryVertex(boundary_vertex).volumeVertex().index());
				}
			}

			std::vector<std::array<double, Dimension>> tiled_vertices;
			std::vector<std::array<std::size_t, Dimension + 1>> tiled_elements;
			for (const auto &vertex : vertices)
			{
				std::array<double, Dimension> point;
				for (std::size_t axis = 0; axis < Dimension; ++axis)
				{
					point[axis] = vertex.point[axis];
				}
				tiled_vertices.push_back(point);
			}

			for (const auto &element : source.elements())
			{
				std::array<std::size_t, Dimension + 1> indices;
				for (std::size_t corner = 0; corner <= Dimension; ++corner)
				{
					indices[corner] = element[corner];
				}
				tiled_elements.push_back(indices);
			}

			std::vector<std::size_t> original(vertices.size());
			std::iota(original.begin(), original.end(), 0);
			for (std::size_t axis = 0; axis < Dimension; ++axis)
			{
				const std::size_t vertex_count = tiled_vertices.size();
				const std::size_t element_count = tiled_elements.size();
				for (std::size_t tile = 1; tile < options.tiles; ++tile)
				{
					std::vector<std::size_t> mapping(vertex_count);
					for (std::size_t vertex = 0; vertex < vertex_count; ++vertex)
					{
						auto point = tiled_vertices[vertex];
						point[axis] += tile * dimensions[axis];
						std::size_t destination = tiled_vertices.size();
						std::size_t group = equivalence.size();
						const auto boundary = topology.vertex(original[vertex]).boundaryVertex();
						if (boundary)
						{
							group = node_set_for_node.at(boundary.index());
							for (const auto candidate : equivalence.at(group))
							{
								double squared_distance = 0;
								for (std::size_t component = 0; component < Dimension; ++component)
								{
									const double difference = point[component] - tiled_vertices[candidate][component];
									squared_distance += difference * difference;
								}

								if (squared_distance < 1e-16)
								{
									destination = candidate;
									break;
								}
							}
						}

						if (destination == tiled_vertices.size())
						{
							tiled_vertices.push_back(point);
							original.push_back(original[vertex]);
							if (group < equivalence.size())
							{
								equivalence[group].push_back(destination);
							}
						}
						mapping[vertex] = destination;
					}

					for (std::size_t element = 0; element < element_count; ++element)
					{
						auto indices = tiled_elements[element];
						for (std::size_t corner = 0; corner <= Dimension; ++corner)
						{
							indices[corner] = mapping[indices[corner]];
						}
						tiled_elements.push_back(indices);
					}
				}
			}

			result.vertices.reserve(tiled_vertices.size() * Dimension);
			for (const auto &vertex : tiled_vertices)
			{
				result.vertices.insert(result.vertices.end(), vertex.begin(), vertex.end());
			}

			result.elements.reserve(tiled_elements.size() * (Dimension + 1));
			for (const auto &element : tiled_elements)
			{
				result.elements.insert(result.elements.end(), element.begin(), element.end());
			}

			result.boundary_shape_velocities.resize(result.shape_parameters.size() * result.vertices.size());
			for (std::size_t parameter = 0; parameter < result.shape_parameters.size(); ++parameter)
			{
				const std::size_t parameter_offset = parameter * result.vertices.size();
				for (std::size_t vertex = 0; vertex < tiled_vertices.size(); ++vertex)
				{
					const std::size_t source_vertex = original[vertex];
					const double normal_velocity = source.normalShapeVelocities()[parameter][source_vertex];
					const auto &normal = source.vertexNormals()[source_vertex];
					const std::size_t vertex_offset = parameter_offset + vertex * Dimension;

					for (std::size_t axis = 0; axis < Dimension; ++axis)
					{
						result.boundary_shape_velocities[vertex_offset + axis] = normal_velocity * normal[axis];
					}
				}
			}
			return result;
		}
	} // namespace

	Mesh inflate(const std::filesystem::path &wire_path, const Options &options,
				 const double *shape_parameters, std::size_t parameter_count)
	{
		static std::mutex mutex;
		const std::lock_guard<std::mutex> lock(mutex);

		const char *type;
		switch (options.cell_type)
		{
		case CellType::Periodic2D:
			type = "2D_doubly_periodic";
			break;
		case CellType::Orthotropic2D:
			type = "2D_orthotropic";
			break;
		case CellType::Orthotropic3D:
			type = "orthotropic";
			break;
		default:
			throw std::invalid_argument("Unsupported inflator cell type.");
		}

		const auto positive = [](double value) { return std::isfinite(value) && value > 0; };
		if (options.tiles == 0 || options.graph_radius == 0 || options.contour_resolution < 2)
		{
			throw std::invalid_argument("Invalid inflator tiling, graph radius, or contour resolution.");
		}

		if (options.cell_type == CellType::Orthotropic3D)
		{
			if (!positive(options.facet_size) || !positive(options.facet_distance)
				|| !positive(options.cell_size) || !positive(options.edge_size))
			{
				throw std::invalid_argument("Meshing sizes must be finite and positive.");
			}
		}
		else if (!positive(options.max_triangle_area))
		{
			throw std::invalid_argument("Meshing sizes must be finite and positive.");
		}

		if (parameter_count == 0 && !positive(options.default_thickness))
		{
			throw std::invalid_argument("Default thickness must be finite and positive.");
		}
		if (parameter_count != 0 && shape_parameters == nullptr)
		{
			throw std::invalid_argument("Nonempty shape parameters require a valid pointer.");
		}

		IsosurfaceInflator source(type, true, wire_path.string(), options.graph_radius);
		auto &meshing = source.meshingOptions();
		meshing.marchingSquaresGridSize = options.contour_resolution;
		meshing.forceMSGridSize = options.force_contour_resolution;
		meshing.maxArea = options.max_triangle_area;
		meshing.facetSize = options.facet_size;
		meshing.facetDistance = options.facet_distance;
		meshing.cellSize = options.cell_size;
		meshing.edgeSize = options.edge_size;

		std::vector<double> parameters;
		if (parameter_count == 0)
		{
			parameters = source.defaultParameters(options.default_thickness);
		}
		else
		{
			if (parameter_count != source.numParams())
			{
				throw std::invalid_argument("Shape parameter count does not match the graph.");
			}

			parameters.assign(shape_parameters, shape_parameters + parameter_count);
			for (const double parameter : parameters)
			{
				if (!std::isfinite(parameter))
				{
					throw std::invalid_argument("Shape parameters must be finite.");
				}
			}
		}

		source.setGenerateFullPeriodCell(true);
		source.inflate(parameters, 0);
		if (source.elements().empty())
		{
			throw std::runtime_error("Inflator produced an empty mesh.");
		}

		if (source.elements().front().size() == 3)
		{
			return copy_and_tile<2>(source, options);
		}
		return copy_and_tile<3>(source, options);
	}
} // namespace inflator
