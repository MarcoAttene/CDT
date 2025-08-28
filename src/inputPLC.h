#ifndef INPUT_PLC_H
#define INPUT_PLC_H

#include <vector>
#include "implicit_point.h"

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

struct input_vertex_t {
public:
    double coord[3];          // Coordinates
    uint32_t original_index;  // Index to support reordering
};

void reorderVertices(double* coordinates, uint32_t numVertices, uint32_t* triVertices, uint32_t numTriangles);
bool misAlignment(const double* p, const double* q, const double* r);
void read_OFF_file(const char* filename, double** vertices_p, uint32_t* npts, uint32_t** tri_vertices_p, uint32_t* ntri, bool verbose);
int vertex_compare(const void* void_v1, const void* void_v2);
int triOrder(const void* t1, const void* t2);
bool coincident_points(const input_vertex_t* a, const input_vertex_t* b);
void remove_duplicated_points(input_vertex_t** vertices_p, uint32_t* npts,
    input_vertex_t* vrts_copy,
    uint32_t* map, uint32_t* diff);
void read_nodes_and_constraints(double* coords_A, uint32_t npts_A, uint32_t* tri_idx_A, uint32_t ntri_A,
    input_vertex_t** vertices_p, uint32_t* npts,
    uint32_t** tri_vertices_p, uint32_t* ntri, bool verbose);

class inputPLC {
public:
    std::vector<double> coordinates; // x1,y1,z1,x2,y2,z2, ..., xn,yn,zn
    std::vector<uint32_t> triangle_vertices; // t1_v1, t1_v2, t1_v3, t2_v1, t2_v2, ...
    const char* input_file_name;

    uint32_t numVertices() const;
    uint32_t numTriangles() const;

    inputPLC();
    inputPLC(const char* filename);

    bool initFromFile(const char* filename, bool verbose);
    bool initFromVectors(double* vertex_p, uint32_t npts, uint32_t* tri_vertices_p, uint32_t ntri, bool verbose);
    void postProcess(double* vertices_p, uint32_t npts, uint32_t* tri_vertices_p, uint32_t ntri, bool verbose);
    void addBoundingBoxVertices();
};

#endif // INPUT_PLC_H