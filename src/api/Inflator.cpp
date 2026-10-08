#include <inflator/Inflator.hpp>

#include <IsosurfaceInflator.hh>
#include <MeshFEM/PeriodicBoundaryMatcher.hh>
#include <MeshFEM/SimplicialMesh.hh>

#include <cmath>
#include <mutex>
#include <numeric>
#include <stdexcept>

namespace inflator
{
    namespace
    {
        template <std::size_t Dimension>
        Mesh copy_and_tile(const IsosurfaceInflator &source, const Request &request)
        {
            Mesh result;
            result.dimension = Dimension;
            result.parameters = source.inflatedParams();
            for (std::size_t parameter = 0; parameter < result.parameters.size(); ++parameter)
            {
                result.parameter_types.push_back(source.isPositionParam(parameter) ? ParameterType::Position
                    : source.isThicknessParam(parameter) ? ParameterType::Thickness
                    : source.isBlendingParam(parameter) ? ParameterType::Blending : ParameterType::Other);
            }

            const auto &vertices = source.vertices();
            const BBox<PointND<Dimension>> bounds(vertices);
            const auto dimensions = bounds.dimensions();
            SimplicialMesh<Dimension> topology(source.elements(), vertices.size());
            std::vector<PointND<Dimension>> boundary_points;
            for (const auto &vertex : topology.boundaryVertices())
            {
                boundary_points.push_back(truncateFrom3D<PointND<Dimension>>(vertices[vertex.volumeVertex().index()].point));
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

            for (const auto &vertex : vertices)
            {
                result.vertices.push_back({vertex.point[0], vertex.point[1], vertex.point[2]});
            }
            for (const auto &element : source.elements())
            {
                std::array<std::size_t, 4> indices{};
                for (std::size_t corner = 0; corner <= Dimension; ++corner)
                {
                    indices[corner] = element[corner];
                }
                result.elements.push_back(indices);
            }
            std::vector<std::size_t> original(vertices.size());
            std::iota(original.begin(), original.end(), 0);
            for (std::size_t axis = 0; axis < Dimension; ++axis)
            {
                const std::size_t vertex_count = result.vertices.size();
                const std::size_t element_count = result.elements.size();
                for (int tile = 1; tile < request.tiles; ++tile)
                {
                    std::vector<std::size_t> mapping(vertex_count);
                    for (std::size_t vertex = 0; vertex < vertex_count; ++vertex)
                    {
                        auto point = result.vertices[vertex];
                        point[axis] += tile * dimensions[axis];
                        std::size_t destination = result.vertices.size();
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
                                    const double difference = point[component] - result.vertices[candidate][component];
                                    squared_distance += difference * difference;
                                }
                                if (squared_distance < 1e-16)
                                {
                                    destination = candidate;
                                    break;
                                }
                            }
                        }
                        if (destination == result.vertices.size())
                        {
                            result.vertices.push_back(point);
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
                        auto indices = result.elements[element];
                        for (std::size_t corner = 0; corner <= Dimension; ++corner)
                        {
                            indices[corner] = mapping[indices[corner]];
                        }
                        result.elements.push_back(indices);
                    }
                }
            }

            result.shape_velocities.resize(result.parameters.size());
            for (std::size_t parameter = 0; parameter < result.parameters.size(); ++parameter)
            {
                auto &velocities = result.shape_velocities[parameter];
                velocities.resize(result.vertices.size());
                for (std::size_t vertex = 0; vertex < result.vertices.size(); ++vertex)
                {
                    for (std::size_t axis = 0; axis < 3; ++axis)
                    {
                        velocities[vertex][axis] = source.normalShapeVelocities()[parameter][original[vertex]]
                            * source.vertexNormals()[original[vertex]][axis];
                    }
                }
            }
            return result;
        }
    }

    Mesh inflate(const Request &request)
    {
        static std::mutex mutex;
        const std::lock_guard<std::mutex> lock(mutex);
        if (request.tiles < 1)
        {
            throw std::invalid_argument("Inflator tile count must be positive.");
        }
        IsosurfaceInflator source(request.type, true, request.wire_path, request.graph_radius);
        source.meshingOptions().load(nlohmann::json::parse(request.meshing_options));
        source.setGenerateFullPeriodCell(true);
        source.inflate(request.parameters.empty() ? source.defaultParameters(request.default_thickness) : request.parameters, 0);
        if (source.elements().empty())
        {
            throw std::runtime_error("Inflator produced an empty mesh.");
        }
        if (source.elements().front().size() == 3)
        {
            return copy_and_tile<2>(source, request);
        }
        return copy_and_tile<3>(source, request);
    }
}
