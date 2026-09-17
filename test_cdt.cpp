// Regression tests for the CDT3 pipeline: tetrahedrize one OFF file and check that the result is a
// sound constrained Delaunay tetrahedrization of it, rather than merely that the process survived.
//
//   test_cdt <file.off> [expected_volume [tolerance]]
//
// Checks, in order of how early they catch damage:
//   1. every tet tagged DT_IN has four finite corners -- no INFINITE_VERTEX, i.e. no ghost was
//      spliced into the solid region,
//   2. every corner's neighbour link points at a tet that exists -- no dangling adjacency,
//   3. the DT_IN region's boundary is closed and manifold: every triangular face is shared by two
//      DT_IN tets or by exactly one (a boundary face), never by three,
//   4. the boundary is consistently oriented -- every directed boundary edge is used exactly once,
//   5. the DT_IN tets all agree on orientation -- none is inverted -- and their volumes sum to the
//      volume the input PLC encloses, when the caller states it.
//
// Prints CDT3 TEST PASSED and returns 0 only if all of them hold. The success line matters: on a
// failure inside the library, ip_error() ends the process -- and upstream's ip_error exits with
// status 0 -- so a test that judged by exit status alone would read a mid-run abort as a pass. The
// CMake test matches the line, not the status.

#ifdef _MSC_VER // Workaround for known bug on MSVC
#define _HAS_STD_BYTE 0
#endif

#include "cdt.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <utility>
#include <vector>

namespace {

int failures = 0;

void check(bool ok, const char* what) {
	if (!ok) { printf("FAILED: %s\n", what); failures++; }
}

// The three corners of face 'c' of tet 't', as a sorted triple, so the two tets sharing a face
// produce the same key.
std::array<uint32_t, 3> faceKey(const TetMesh& tin, uint64_t t, int c) {
	const uint32_t* n = tin.tet_node.data() + (t << 2);
	std::array<uint32_t, 3> f = { n[(c + 1) & 3], n[(c + 2) & 3], n[(c + 3) & 3] };
	std::sort(f.begin(), f.end());
	return f;
}

// The same face with its winding preserved, seen from outside the tet.
std::array<uint32_t, 3> faceLoop(const TetMesh& tin, uint64_t t, int c) {
	const uint32_t* n = tin.tet_node.data() + (t << 2);
	// Corners 0 and 2 see the opposite face with one orientation, 1 and 3 with the other.
	if (c & 1) return { n[(c + 1) & 3], n[(c + 2) & 3], n[(c + 3) & 3] };
	return { n[(c + 3) & 3], n[(c + 2) & 3], n[(c + 1) & 3] };
}

double tetVolume(const TetMesh& tin, uint64_t t) {
	const uint32_t* n = tin.tet_node.data() + (t << 2);
	double p[4][3];
	for (int i = 0; i < 4; i++) tin.vertices[n[i]]->getApproxXYZCoordinates(p[i][0], p[i][1], p[i][2], true);
	const double a[3] = { p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2] };
	const double b[3] = { p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2] };
	const double c[3] = { p[3][0] - p[0][0], p[3][1] - p[0][1], p[3][2] - p[0][2] };
	return (a[0] * (b[1] * c[2] - b[2] * c[1])
	      - a[1] * (b[0] * c[2] - b[2] * c[0])
	      + a[2] * (b[0] * c[1] - b[1] * c[0])) / 6.0;
}

} // namespace

int main(int argc, char** argv) {
	if (argc < 2 || argc > 4) {
		fprintf(stderr, "usage: %s <file.off> [expected_volume [tolerance]]\n", argv[0]);
		return 2;
	}
	const bool have_expected = (argc >= 3);
	const double expected_volume = have_expected ? atof(argv[2]) : 0.0;
	const double tolerance = (argc >= 4) ? atof(argv[3]) : 1e-9;

	// Read the OFF file straight into the tetrahedrizer, WITHOUT inputPLC::postProcess.
	//
	// postProcess() de-duplicates vertices, drops degenerate triangles, and then spatially reorders
	// the vertices. That reordering changes the insertion order, so it changes the Delaunay
	// triangulation, which faces need recovering, and which of them carry vertices lying in their
	// own plane -- the flat-cavity-face regression below does not reproduce through it at all. It is
	// also how a caller embedding this library drives it when the vertex indices it passes in have
	// meaning to it and must survive: TetMesh::init_vertices + PLCx is the public way to do that.
	// So the tests exercise the input exactly as written, which is also the harder case, since
	// nothing has pre-conditioned the point order in the tetrahedrizer's favour.
	//
	// It puts the burden of a clean PLC on the input file: the ones here are closed, manifold and
	// free of duplicate vertices and degenerate triangles, which the audit below would report.
	double* coordinates = nullptr;
	uint32_t* triangle_vertices = nullptr;
	uint32_t num_vertices = 0, num_triangles = 0;
	read_OFF_file(argv[1], &coordinates, &num_vertices, &triangle_vertices, &num_triangles, false);
	if (num_vertices == 0 || num_triangles == 0) { printf("FAILED: empty input\n"); return 1; }
	printf("input: %u vertices, %u triangles\n", num_vertices, num_triangles);

	TetMesh tin;
	tin.init_vertices(coordinates, num_vertices);
	tin.tetrahedrize();

	PLCx steiner_plc(tin, triangle_vertices, num_triangles);
	steiner_plc.segmentRecovery_HSi(true /*quiet*/);
	steiner_plc.faceRecovery(true /*quiet*/);
	const size_t num_inner = steiner_plc.markInnerTets();

	printf("result: %u tets (%zu inner), %u vertices (%u Steiner)\n",
	       tin.numTets(), num_inner, tin.numVertices(), steiner_plc.numSteinerVertices());

	check(num_inner > 0, "the PLC encloses no tetrahedra");

	// 1+2: no ghost and no dangling adjacency inside the solid region.
	//
	// This is what the flat-cavity-face regression used to break. meshCavity() would reconnect a
	// kept tet to a cavity tet it had just turned into a ghost and was about to delete, so the mesh
	// kept a neighbour link into storage that the following truncation released; a later walk
	// followed it and read INFINITE_VERTEX as a vertex index. Nothing faulted at the moment of the
	// damage, only much later, which is why these two invariants are checked over the whole mesh
	// rather than trusted.
	size_t ghost_in_solid = 0, dangling = 0;
	for (uint64_t t = 0; t < tin.numTets(); t++) {
		if (tin.mark_tetrahedra[t] == DT_IN)
			for (int i = 0; i < 4; i++)
				if (tin.tet_node[(t << 2) + i] == INFINITE_VERTEX) ghost_in_solid++;
		for (int i = 0; i < 4; i++)
			if ((tin.tet_neigh[(t << 2) + i] >> 2) >= tin.numTets()) dangling++;
	}
	check(ghost_in_solid == 0, "a tet marked DT_IN has an infinite vertex");
	check(dangling == 0, "a corner's neighbour link points past the end of the mesh");
	if (ghost_in_solid || dangling)
		printf("  (%zu ghost corner(s) in solid tets, %zu dangling neighbour link(s))\n",
		       ghost_in_solid, dangling);

	// 3+4: the solid region is bounded by a closed, manifold, consistently oriented surface.
	std::map<std::array<uint32_t, 3>, int> face_use;
	std::vector<std::array<uint32_t, 3>> boundary;
	for (uint64_t t = 0; t < tin.numTets(); t++) {
		if (tin.mark_tetrahedra[t] != DT_IN) continue;
		for (int i = 0; i < 4; i++) {
			face_use[faceKey(tin, t, i)]++;
			const uint64_t n = tin.tet_neigh[(t << 2) + i] >> 2;
			if (n >= tin.numTets() || tin.mark_tetrahedra[n] != DT_IN) boundary.push_back(faceLoop(tin, t, i));
		}
	}
	size_t overused = 0;
	for (const auto& fu : face_use) if (fu.second > 2) overused++;
	check(overused == 0, "a triangular face is shared by more than two solid tets (non-manifold)");

	std::map<std::pair<uint32_t, uint32_t>, int> directed;
	for (const auto& f : boundary)
		for (int i = 0; i < 3; i++) directed[{ f[i], f[(i + 1) % 3] }]++;
	size_t repeated = 0, unpaired = 0;
	for (const auto& de : directed) {
		if (de.second != 1) repeated++;
		if (directed.find({ de.first.second, de.first.first }) == directed.end()) unpaired++;
	}
	check(repeated == 0, "a directed boundary edge is used more than once (inconsistent winding)");
	check(unpaired == 0, "a boundary edge has no opposite (the solid boundary is not closed)");
	printf("boundary: %zu triangles, %zu directed edges\n", boundary.size(), directed.size());

	// 5: consistent orientation, and the volume the PLC encloses.
	//
	// Signs are counted rather than required to be positive: which way round CDT3 winds a tet's
	// corners is its own convention, and pinning it here would make the test fail on a convention
	// change instead of on a defect. What must hold is that every solid tet agrees -- one inverted
	// tet among them is a folded mesh -- so the enclosed volume is the magnitude of the sum.
	double volume = 0.0;
	size_t positive = 0, negative = 0, degenerate = 0;
	for (uint64_t t = 0; t < tin.numTets(); t++) {
		if (tin.mark_tetrahedra[t] != DT_IN) continue;
		const double v = tetVolume(tin, t);
		if (v > 0.0) positive++; else if (v < 0.0) negative++; else degenerate++;
		volume += v;
	}
	check(degenerate == 0, "a solid tet has zero volume");
	check(positive == 0 || negative == 0, "solid tets disagree on orientation (the mesh folds over)");
	if (positive && negative) printf("  (%zu positive, %zu negative)\n", positive, negative);
	volume = std::fabs(volume);
	printf("volume: %.17g\n", volume);
	if (have_expected) {
		const bool ok = std::fabs(volume - expected_volume) <= tolerance;
		if (!ok) printf("  expected %.17g (tolerance %g), differs by %g\n",
		                expected_volume, tolerance, std::fabs(volume - expected_volume));
		check(ok, "the tetrahedrized volume differs from the volume the PLC encloses");
	}

	// PLCx borrows the triangle array for its whole life, so both survive until here.
	free(coordinates);
	free(triangle_vertices);

	if (failures) { printf("CDT3 TEST FAILED (%d check(s))\n", failures); return 1; }
	printf("CDT3 TEST PASSED\n");
	return 0;
}
