# Library interface

`inflator/Inflator.hpp` is the only installed public header. It uses C++17
standard-library types; Eigen, CGAL, MeshFEM, Triangle, and JSON types are private.
The shared library hides implementation symbols. Use a compatible C++ standard
library/ABI in the application and inflator.

```cpp
#include <inflator/Inflator.hpp>

inflator::Request request;
request.type = "orthotropic";
request.wire_path = "cell.obj";
request.meshing_options = R"({"cellSize":0.15,"facetSize":0.025})";
const auto mesh = inflator::inflate(request);
```

`inflate` is synchronous and throws a standard exception on failure. Calls are
serialized inside the library because the existing implementation uses global
benchmarking state. CGAL meshing starts from a fixed random seed for reproducible
evaluations. Mesh-quality optimization uses convergence rather than CGAL's
default CPU-time cutoff, which otherwise changes output between evaluations.
No subprocess or intermediate mesh/velocity file is required.

## Inputs

- `type`: existing inflator symmetry/mesher name. Use `2D_doubly_periodic` for
  unrestricted planar periodic graphs; use `orthotropic` or `cubic` for reflected
  3D cells with matching periodic face triangulations. Arbitrary
  `triply_periodic` CGAL meshes do not guarantee matching interfaces.
- `wire_path`: OBJ line graph. Supply the whole graph, including negative
  coordinates; the existing implementation normalizes its bounding box to
  `[-1,1]` before extracting the symmetry base graph.
- `parameters`: graph parameters in the original inflator ordering. Empty
  selects defaults, with thickness set by `default_thickness`.
- `meshing_options`: a JSON string using the existing meshing option names.
  2D primarily uses `maxArea`; 3D uses `facetSize`, `facetDistance`, `cellSize`,
  `edgeSize`, and `cellRadiusEdgeRatio`.
- `tiles`: repetitions along every active axis. Opposite cell boundaries must
  match; tiling stitches shared vertices rather than concatenating meshes.
- `graph_radius`: graph-edge neighborhood used for implicit geometry evaluation.

## Outputs

- `dimension`: 2 for triangles, 3 for tetrahedra.
- `vertices`: three coordinates per vertex; planar meshes have zero z.
- `elements`: four zero-based indices per simplex; triangles use only the first
  three entries.
- `parameters` and `parameter_types`: actual parameter values and classification.
- `shape_velocities[parameter][vertex][axis]`: normal boundary velocity,
  computed from the implicit function derivative and propagated through
  reflection/tiling. Interior velocities are zero. This is a geometric shape
  derivative, not a derivative of the discrete remeshing algorithm.

## Build and consume

```sh
cmake --preset library
cmake --build build/library --target inflator_api_test
ctest --test-dir build/library --output-on-failure
cmake --install build/library --prefix /path/to/install
```

The `sanitize` preset enables AddressSanitizer and UndefinedBehaviorSanitizer.
CGAL 5.6.2 replaces 4.12, whose tetrahedral vertex-removal path invalidates
unordered-map iterators with current Boost and fails the 3D API test under ASan.

Consumers use `find_package(inflator CONFIG REQUIRED)` and link
`inflator::inflator`. Only the public include directory, C++17 requirement, and
shared library propagate. The install does not export private dependency targets.
Build the inflator separately (or with CMake ExternalProject) to isolate its
legacy dependency configuration from the application.

GMP, MPFR, Boost headers, a C++17 compiler, CMake, and Ninja must be available.
The existing dependency recipes download CGAL, Eigen, and other source
dependencies automatically. CLI programs remain available with
`MICRO_BUILD_BINARIES=ON`; library consumers do not require them.

The GitHub Actions `Library` workflow checks GCC and Clang on Ubuntu 24.04,
including a separately configured installed-package consumer and the exported
symbol/header boundary. Run it locally with `act -W .github/workflows/library.yml
-P ubuntu-24.04=catthehacker/ubuntu:act-latest --concurrent-jobs 1`.
These Linux container tests do not establish macOS or Windows compatibility.
