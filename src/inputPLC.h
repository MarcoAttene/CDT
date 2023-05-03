#include <iostream>
#include <fstream>

using namespace std;

class reVertex {
public:
    double c[3];
    uint32_t index;

    reVertex(double* _c, uint32_t _i) : c{ _c[0], _c[1], _c[2] }, index(_i) {}

    static bool lessThanOnX(const reVertex& a, const reVertex& b) { return (a.c[0] < b.c[0]); }
    static bool lessThanOnY(const reVertex& a, const reVertex& b) { return (a.c[1] < b.c[1]); }
    static bool lessThanOnZ(const reVertex& a, const reVertex& b) { return (a.c[2] < b.c[2]); }
};

class vBlock {
public:
    uint32_t begin, end;
    uint32_t dir_split;

    vBlock(uint32_t b, uint32_t e, uint32_t d) : begin(b), end(e), dir_split(d) {}
};

void reorderVertices(double* coordinates, uint32_t numVertices, uint32_t* triVertices, uint32_t numTriangles) {
    std::vector<reVertex> revertices; 
    revertices.reserve(numVertices);
    for (uint32_t i = 0; i < numVertices; i++) revertices.push_back(reVertex(coordinates + (i * 3), i));

    std::vector<vBlock> blocks;
    blocks.push_back(vBlock(0, numVertices, 0));

    for (uint32_t b = 0; b < blocks.size(); b++) {
        const vBlock block = blocks[b];
        const uint32_t dir_split = block.dir_split;
        if (dir_split == 0) std::sort(revertices.begin() + block.begin, revertices.begin() + block.end, reVertex::lessThanOnX);
        else if (dir_split == 1) std::sort(revertices.begin() + block.begin, revertices.begin() + block.end, reVertex::lessThanOnY);
        else std::sort(revertices.begin() + block.begin, revertices.begin() + block.end, reVertex::lessThanOnZ);
        const uint32_t block_mid = (block.begin + block.end) >> 1;
        if (block_mid != block.begin && block_mid != block.end) {
            blocks.push_back(vBlock(block.begin, block_mid, (dir_split + 1) % 3));
            blocks.push_back(vBlock(block_mid, block.end, (dir_split + 1) % 3));
        }
    }

    std::vector<uint32_t> new_index(numVertices);
    for (uint32_t i = 0; i < numVertices; i++) new_index[revertices[i].index] = i;

    for (uint32_t i = 0; i < numTriangles * 3; i++) triVertices[i] = new_index[triVertices[i]];
    for (uint32_t i = 0; i < numVertices; i++) {
        const double* pt = revertices[i].c;
        coordinates[i * 3] = pt[0];
        coordinates[i * 3 + 1] = pt[1];
        coordinates[i * 3 + 2] = pt[2];
    }
}

struct input_vertex_t {
public:
    double coord[3];          // Coordinates
    uint32_t original_index;  // Index to support reordering
};

bool misAlignment(const double* p, const double* q, const double* r)
{
    return orient2d(p[0], p[1], q[0], q[1], r[0], r[1]) ||
        orient2d(p[1], p[2], q[1], q[2], r[1], r[2]) ||
        orient2d(p[0], p[2], q[0], q[2], r[0], r[2]);
}

void read_OFF_file(const char* filename,
    double** vertices_p, uint32_t* npts,
    uint32_t** tri_vertices_p, uint32_t* ntri, bool verbose) {

    FILE* file = fopen(filename, "r");
    if (file == NULL)
        ip_error("read_OFF_file: FATAL ERROR "
            "cannot open input file.\n");

    // Check OFF mark (1st line).
    char file_ext_read[3];
    char file_ext_target[] = { 'O','F','F' };
    if (fscanf(file, "%3c", file_ext_read) == 0)
        ip_error("read_OFF_file: FATAL ERROR "
            "cannot read 1st line of input file\n");

    for (uint32_t i = 0; i < 3; i++)
        if (file_ext_read[i] != file_ext_target[i])
            ip_error("read_OFF_file: FATAL ERROR "
                "1st line of input file is different from OFF\n");

    // Reading number of points and triangles.
    if (fscanf(file, " %d %d %*d ", npts, ntri) == 0)
        ip_error("read_OFF_file: FATAL ERROR 2st line of "
            "input file do not contanins point and triangles numbers.\n");

    if (verbose) std::cout << "file " << filename << " contains " << *npts << " vertices and " << *ntri << " constraints (triangles)\n";

    // Reading points coordinates.
    *vertices_p = (double*)malloc(sizeof(double) * 3 * (*npts));
    *tri_vertices_p = (uint32_t*)malloc(sizeof(uint32_t) * 3 * (*ntri));

    for (uint32_t i = 0; i < (*npts); i++) {
        if (fscanf(file, " %lf %lf %lf ",
            (*vertices_p) + (i * 3), (*vertices_p) + (i * 3 + 1), (*vertices_p) + (i * 3 + 2)) == 0)
            ip_error("error reading input file\n");
    }

    uint32_t nv;
    for (uint32_t i = 0; i < (*ntri); i++) {
        if (fscanf(file, " %u %u %u %u ", &nv,
            (*tri_vertices_p) + (i * 3), (*tri_vertices_p) + (i * 3 + 1), (*tri_vertices_p) + (i * 3 + 2)) == 0)
            ip_error("error reading input file\n");
        if (nv != 3) ip_error("Non-triangular faces not supported\n");
    }
    fclose(file);
}

int vertex_compare(const void* void_v1, const void* void_v2)
{
    const input_vertex_t* v1 = (input_vertex_t*)void_v1;
    const input_vertex_t* v2 = (input_vertex_t*)void_v2;
    const double dx = v1->coord[0] - v2->coord[0];
    const double dy = v1->coord[1] - v2->coord[1];
    const double dz = v1->coord[2] - v2->coord[2];
    return (4 * ((dx > 0) - (dx < 0)) +
        2 * ((dy > 0) - (dy < 0)) +
        ((dz > 0) - (dz < 0)));
}

int triOrder(const void* t1, const void* t2) {
    const uint32_t* a = (uint32_t*)t1;
    const uint32_t* b = (uint32_t*)t2;

    // Here we should pre-order a[3] and b[3] to identify coincident triangles with different vertex ordering!!!
    // To be done!
    if (a[0] < b[0]) return -1;
    if (a[0] > b[0]) return 1;
    if (a[1] < b[1]) return -1;
    if (a[1] > b[1]) return 1;
    if (a[2] < b[2]) return -1;
    return (a[2] > b[2]);
}

inline bool coincident_points(const input_vertex_t* a, const input_vertex_t* b)
{
    return !vertex_compare(a->coord, b->coord);
}

void remove_duplicated_points(input_vertex_t** vertices_p, uint32_t* npts,
    input_vertex_t* vrts_copy,
    uint32_t* map, uint32_t* diff) {

    // Sorting vertices by coordinates lexicographic order (x,y,z)
    qsort(vrts_copy, *npts, sizeof(input_vertex_t), vertex_compare);

    // Count and memory position of duplicated vertices.
    uint32_t vrts_counter = 0;
    for (uint32_t i = 1; i < (*npts); i++) {
        if (coincident_points(vrts_copy + i - 1, vrts_copy + i)) vrts_counter++;
        diff[i] = vrts_counter;
    }

    // Set original_index to follow vertices permutation.
    for (uint32_t i = 0; i < (*npts); i++)  map[vrts_copy[i].original_index] = i;

    // Allocating memory to store uinque mesh vertices (vertices_p).
    *vertices_p = (input_vertex_t*)malloc(sizeof(input_vertex_t) * (*npts - vrts_counter));

    // Fill mesh vertices (vertices_p)
    memcpy(*vertices_p, vrts_copy, sizeof(input_vertex_t));
    for (uint32_t i = vrts_counter = 1; i < (*npts); i++)
        if (!coincident_points(vrts_copy + i - 1, vrts_copy + i))
            memcpy(*vertices_p + (vrts_counter++), vrts_copy + i, sizeof(input_vertex_t));

    // Update uinque vertices number.
    (*npts) = vrts_counter;
}


/// //////////////////////////////////////////////////////////////////////////////////////////

void read_nodes_and_constraints(double* coords_A, uint32_t npts_A, uint32_t* tri_idx_A, uint32_t ntri_A,
    input_vertex_t** vertices_p, uint32_t* npts,
    uint32_t** tri_vertices_p, uint32_t* ntri, bool verbose) {

    // Reading points coordinates.
    *npts = npts_A;
    *ntri = ntri_A;
    input_vertex_t* tmp = (input_vertex_t*)malloc(*npts * sizeof(input_vertex_t));
    uint32_t* diff = (uint32_t*)calloc(*npts, sizeof(uint32_t));
    uint32_t* map = (uint32_t*)malloc(*npts * sizeof(uint32_t));
    *tri_vertices_p = (uint32_t*)malloc(sizeof(uint32_t) * 3 * (*ntri));

    for (uint32_t i = 0; i < (*npts); i++) {
        tmp[i].coord[0] = coords_A[i * 3];
        tmp[i].coord[1] = coords_A[i * 3 + 1];
        tmp[i].coord[2] = coords_A[i * 3 + 2];
        tmp[i].original_index = i;
    }

    if (*npts > 1) remove_duplicated_points(vertices_p, npts, tmp, map, diff);
    if (verbose) std::cout << "Using " << *npts << " unique vertices\n";

    free(tmp);

    for (uint32_t i = 0, j = 0; i < (*ntri); j++) {

        const uint32_t i1 = tri_idx_A[j * 3];
        const uint32_t i2 = tri_idx_A[j * 3 + 1];
        const uint32_t i3 = tri_idx_A[j * 3 + 2];
        (*tri_vertices_p)[3 * i] = map[i1] - diff[map[i1]];
        (*tri_vertices_p)[3 * i + 1] = map[i2] - diff[map[i2]];
        (*tri_vertices_p)[3 * i + 2] = map[i3] - diff[map[i3]];

        const double* v1c = ((*vertices_p) + (*tri_vertices_p)[3 * i])->coord;
        const double* v2c = ((*vertices_p) + (*tri_vertices_p)[3 * i + 1])->coord;
        const double* v3c = ((*vertices_p) + (*tri_vertices_p)[3 * i + 2])->coord;

        if (!misAlignment(v1c, v2c, v3c))
        {
            ip_error("Model has degenerate triangles. Unsupported!\n");
            (*ntri)--;
        }
        else i++;
    }
    free(map);
    free(diff);

    qsort((*tri_vertices_p), *ntri, sizeof(uint32_t) * 3, triOrder);
    uint32_t num_duplicates = 0;
    for (uint32_t i = 0; i < *ntri - 1; i++) if (triOrder((*tri_vertices_p) + i * 3, (*tri_vertices_p) + (i + 1) * 3) == 0) {
        (*tri_vertices_p)[i * 3] = (*tri_vertices_p)[i * 3 + 1] = (*tri_vertices_p)[i * 3 + 2] = UINT32_MAX;
        num_duplicates++;
    }
    qsort((*tri_vertices_p), *ntri, sizeof(uint32_t) * 3, triOrder);
    *ntri -= num_duplicates;

    if (verbose) std::cout << "Using " << *ntri << " non-degenerate constraints\n";
}

class inputPLC {
public:
    std::vector<double> coordinates; // x1,y1,z1,x2,y2,z2, ..., xn,yn,zn
    std::vector<uint32_t> triangle_vertices; // t1_v1, t1_v2, t1_v3, t2_v1, t2_v2, ...
    const char* input_file_name;

    uint32_t numVertices() const { return (uint32_t)coordinates.size() / 3; }
    uint32_t numTriangles() const { return (uint32_t)triangle_vertices.size() / 3; }

    inputPLC() {}

    inputPLC(const char* filename) { initFromFile(filename, true); }

    bool initFromFile(const char* filename, bool verbose) {
        input_file_name = filename;

        // Read OFF file (only triangles supported)
        uint32_t npts, ntri;
        double* vertex_p;
        uint32_t* tri_vertices_p;
        read_OFF_file(filename, &vertex_p, &npts, &tri_vertices_p, &ntri, verbose);
        if (verbose) printf("File read\n");
        if (npts == 0) ip_error("Input file has no vertices\n");
        if (ntri == 0) ip_error("Input file has no triangles\n");

        postProcess(vertex_p, npts, tri_vertices_p, ntri, verbose);

        free(vertex_p);
        free(tri_vertices_p);

        return true;
    }

    bool initFromVectors(double* vertex_p, uint32_t npts, uint32_t* tri_vertices_p, uint32_t ntri, bool verbose) {
        input_file_name = "";

        postProcess(vertex_p, npts, tri_vertices_p, ntri, verbose);

        return true;
    }

    void postProcess(double* vertices_p, uint32_t npts, uint32_t* tri_vertices_p, uint32_t ntri, bool verbose) {
        // Convert OFF to valid set of vertices (no duplications) and constraints (no degeneracies)
        uint32_t* valid_tri_vertices_p;
        uint32_t num_valid_tris;
        input_vertex_t* tmp_vertices; // These have floating point coordinates
        uint32_t num_vertices;
        read_nodes_and_constraints(vertices_p, npts, tri_vertices_p, ntri, &tmp_vertices, &num_vertices, &valid_tri_vertices_p, &num_valid_tris, verbose);
        if (verbose) printf("Valid input built\n");
        if (num_vertices == 0) ip_error("Input file has no valid vertices\n");
        if (num_valid_tris == 0) ip_error("Input file has no valid triangles\n");

        coordinates.resize(3 * num_vertices);
        for (uint32_t i = 0; i < num_vertices; i++) {
            coordinates[i * 3] = tmp_vertices[i].coord[0];
            coordinates[i * 3 + 1] = tmp_vertices[i].coord[1];
            coordinates[i * 3 + 2] = tmp_vertices[i].coord[2];
        }

        triangle_vertices.assign(valid_tri_vertices_p, valid_tri_vertices_p + (num_valid_tris * 3));

        free(tmp_vertices);
        free(valid_tri_vertices_p);

        reorderVertices(coordinates.data(), (uint32_t)coordinates.size() / 3, triangle_vertices.data(), (uint32_t)triangle_vertices.size() / 3);
    }

    // Add eight vertices to enclose the input in a box
    void addBoundingBoxVertices() {
        double bbmin[3] = { DBL_MAX, DBL_MAX, DBL_MAX };
        double bbmax[3] = { -DBL_MAX, -DBL_MAX, -DBL_MAX };
        for (uint32_t i = 0; i < numVertices(); i++) {
            const double *v = coordinates.data() + i * 3;
            for (int j = 0; j < 3; j++) {
                if (v[j] < bbmin[j]) bbmin[j] = v[j];
                if (v[j] > bbmax[j]) bbmax[j] = v[j];
            }
        }
        const double bbox[3] = { bbmax[0] - bbmin[0], bbmax[1] - bbmin[1], bbmax[2] - bbmin[2] };
        for (int j = 0; j < 3; j++) {
            bbmin[j] -= bbox[j] * 0.05;
            bbmax[j] += bbox[j] * 0.05;
        }

        const int idx[] = { 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1 };

        for (int j = 0; j < 24; j++)
            coordinates.push_back(idx[j] ? (bbmax[j%3]) : (bbmin[j%3]));
    }
};
