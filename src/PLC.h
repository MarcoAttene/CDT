#ifndef _PLC_
#define _PLC_

#include "delaunay.h"
#include <cstring>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include <fstream>
#include <set>
#include <algorithm>
#include <iomanip>

#pragma intrinsic(fabs)

// NOTES: 1) "both_acute_ep" edges will be immediatelly split by inserting the middle point (each subedge becomes a "one_acute_ep")
//        2) sub-edges of "no_acute_ep" and "one_acute_ep" inherit type
//        3) "flat" edges will be ignored by segment recovery algorithm and will not be further classyfied
typedef enum{
	undet,
	no_acute_ep, one_acute_ep, both_acute_ep, flat
} PLCedge_type;

// An edge or sub-edge of the input triangle mesh (PLC). It is a segment.
// ep[2] are the two endpoints of the edge.
// In case the edge is a sub-edge (at least one endpint is a Steiner point)
// oep[2] are the endpoints of the parent-edge. Otherwise oep[2] = ep[2].
// type=2 edeges are those that have an acute endpoint (i.e.) exists
// another incident edge which forms an acute angle.
// For type=2 edges:
// - oep[0] is the acute vertex,
// - ep[0] = oep[0] OR ope[0]<ep[0]<ep[1] along the straight line for the edge,
// - ep[1] = oep[2] OR ep[0]<ep[1]<oep[1] along the straight line for the edge.
// The above rules do not hold for type=3 edges: they are all split first and
// their sub-edges (type=2) have oep[2] = ep[2].
class PLCedge{
public:
    uint32_t ep[2];                 // Endpoints (vertices-inds wrt TetMesh vertices vector)
    uint32_t oep[2];                // Parent-edge endpoints (vertices-inds wrt TetMesh vertices vector)
    std::vector<uint32_t> inc_tri;  // Incident triangles (triangle-inds wrt input triangles vector)
    PLCedge_type type;              // See notes after definition of PLCedge_type.

    inline PLCedge() {}

    inline PLCedge(const uint32_t e0, const uint32_t e1, const uint32_t oe0, const uint32_t oe1,
        const std::vector<uint32_t>& itri, const PLCedge_type t) : ep{ e0, e1 }, oep{oe0, oe1}, inc_tri(itri), type(t) {}

    bool isFlat() const { return type == flat; }

    inline void swap() { std::swap(ep[0], ep[1]); std::swap(oep[0], oep[1]); }

    inline void fill_preEdge(const uint32_t v0, const uint32_t v1, const uint32_t fi) {
        if (v0 > v1) { oep[0] = ep[0] = v1; oep[1] = ep[1] = v0; }
        else { oep[0] = ep[0] = v0; oep[1] = ep[1] = v1; }
        inc_tri.push_back(fi);
    }

    void replaceIncidentFace(uint32_t old_f, uint32_t new_f) {
        std::replace(inc_tri.begin(), inc_tri.end(), old_f, new_f);
    }

    uint32_t commonVertex(const PLCedge& e) const {
        if (ep[0] == e.ep[0] || ep[0] == e.ep[1]) return ep[0];
        else if (ep[1] == e.ep[0] || ep[1] == e.ep[1]) return ep[1];
        else return UINT32_MAX;
    }

    bool coincident(const PLCedge& e) const {
        return (ep[0] == e.ep[0] && ep[1] == e.ep[1]) || (ep[0] == e.ep[1] && ep[1] == e.ep[0]);
    }

    uint32_t oppositeVertex(uint32_t v) const {
        if (ep[0] == v) return ep[1];
        assert(ep[1] == v);
        return ep[0];
    }

    uint32_t commonOriginalVertex(const PLCedge& e) const {
        if (oep[0] == e.oep[0] || oep[0] == e.oep[1]) return oep[0];
        assert(oep[1] == e.oep[0] || oep[1] == e.oep[1]);
        return oep[1];
    }

    uint32_t oppositeOriginalVertex(uint32_t v) const {
        if (oep[0] == v) return oep[1];
        assert(oep[1] == v);
        return oep[0];
    }

    inline bool isIsolated() const { return (inc_tri.empty()); }

    inline bool hasOriginalVertex(uint32_t v) const {
        return (oep[0] == v && oep[1] == v);
    }

    inline bool hasOriginalVertices(uint32_t v1, uint32_t v2) const {
        return (oep[0] == v1 && oep[1] == v2) || (oep[0] == v2 && oep[1] == v1);
    }

    // Static functions to be used as predicates in std algorithms
    static inline bool isIsolatedPtr(const PLCedge& e) { return e.isIsolated(); }

    static inline bool vertexSortFunc(const PLCedge& e1, const PLCedge& e2) {
        if (e1.ep[0] == e2.ep[0]) return (e1.ep[1] < e2.ep[1]);
        else return (e1.ep[0] < e2.ep[0]);
    }

    static inline bool edgeSortFuncPtr(const PLCedge *e1, const PLCedge *e2) {
        return e1 < e2;
    }
};


/// <summary>
/// PLCface
/// This is a maximal flat face of the input PLC
/// It is built out of edge-adjacent and coplanar input triangles
/// It might be bounded by one or more loops of PLCedges
/// </summary>

class PLCface {
public:
    std::vector<uint32_t> triangles; // Original triangles composing the face
    std::vector<PLCedge *> bounding_edges; // Set of bounding edges
    std::vector<std::pair<uint32_t, uint32_t>> orig_flat_edges; // Original flat edges
    std::vector<uint32_t> vertices; // Ordered mesh vertices bounding the face (see savePLC)
    std::vector<uint32_t> flat_vertices; // Face internal vertices (having a flat neghborhood)
    int max_comp_normal;    // Max component of face normal (0=x, 1=y, 2=z)
    bool is_convex;
    bool is_simply_connected;

    PLCface() {}

    void zip();

    void initConvexity(const class PLCx& plc);

    void replaceEdge(PLCedge* old_e, PLCedge* new_e) {
        std::replace(bounding_edges.begin(), bounding_edges.end(), old_e, new_e);
    }

    void makeVertices();

    void absorb(PLCface& f, PLCedge* e);

    void replaceIncidentEdgeFaces(uint32_t old_f, uint32_t new_f) {
        for (PLCedge* e : bounding_edges) e->replaceIncidentFace(old_f, new_f);
    }

    static inline bool isEmpty(const PLCface& e) { return e.bounding_edges.empty(); }
};


#define UNDET_ORIENTATION   -2

class PLCx{
public:
  const size_t input_nv; // number of input vertices
  const uint32_t input_nt; // number of input triangles
  const uint32_t* input_tv; // input triangles (linearized vertex IDs)

  TetMesh& delmesh; // Delaunay tetrahedrization
  std::vector<PLCedge> edges; // edges of the PLC

  std::vector<int> v_orient; // Pre-computed orientation of vertices wrt one plane
  std::vector<std::vector<std::vector<uint32_t>>> vt_maps; // Set of input triangles incident upon each vertex

  std::vector<PLCface> faces; // Faces of the PLC
  std::vector<uint32_t> v_reindex; // Maps global vertex IDs to local indexes into face vertex vectors
  std::vector<std::pair<uint32_t, uint32_t>> singular_v; // Set of pairs <global_vertex_ID, local_face_vertex_ID>, one per singular vertex in face

  bool is_polyhedron; // TRUE if all the PLC edges have an even number of incident faces


  PLCx(TetMesh& m, const uint32_t* _input_tv, const uint32_t _input_nt) :
      input_nv(m.vertices.size()), input_nt(_input_nt), input_tv(_input_tv), delmesh(m), is_polyhedron(false)
  { initialize(); };

  void initialize();
  void mergePreEdges(); // Removes duplicated pre-edges

  void pushVertex(pointType* p, uint32_t acute_v_id) {
      delmesh.pushVertex(p);
      delmesh.marked_vertex.push_back(0);
  }

  // For each face, for each of its vertices, set of incident face triangles
  void makeVertexTriangleMaps(std::vector<std::vector<std::vector<uint32_t>>>& vt_maps);

  void makeVertexTriangleMap(PLCface& f, std::vector<std::vector<uint32_t>>& vt_map,
      std::vector<bool>& orig_tri_mark);


  // Access functions

  bool isSteinerVertex(uint32_t v) const { return v >= input_nv; }

  uint32_t numSteinerVertices() const { return (uint32_t)(delmesh.vertices.size() - input_nv); }

  bool isAcute(const uint32_t vi, const std::vector<std::vector<uint32_t>>& vv) const;

  bool isFlat(const PLCedge& e) const {
      return e.inc_tri.size() == 2 && delmesh.vOrient3D(e.ep[0], e.ep[1], opposite_vrt(e, e.inc_tri[0]), opposite_vrt(e, e.inc_tri[1])) == 0;
  }

  uint32_t opposite_vrt(const PLCedge& e, const uint32_t ti) const;

  bool faceHasTriangle(const PLCface& f, const uint32_t tv[3]) const;

  // Functions to calculate Steiner point positions
  double getT1(uint32_t oe0i, uint32_t e0i) const;
  double getT2(uint32_t oe1i, uint32_t e1i) const;
  inline implicitPoint_LNC* getProjectionOrMidPoint(uint32_t oe0i, uint32_t oe1i, uint32_t e0i, uint32_t e1i, uint32_t ri, uint32_t& acute_v) const;
  inline implicitPoint_LNC* getProjectionOrMidPoint_noac(uint32_t oe0i, uint32_t oe1i, uint32_t e0i, uint32_t e1i, uint32_t ri, uint32_t& acute_v) const;
  inline implicitPoint_LNC* getProjectionOrMidPoint_noac_rev(uint32_t oe0i, uint32_t oe1i, uint32_t e0i, uint32_t e1i, uint32_t ri, uint32_t& acute_v) const;
  inline implicitPoint_LNC* getMidPoint(uint32_t oe0i, uint32_t oe1i, uint32_t e0i, uint32_t e1i) const;

  bool is_missing_PLCedge(const uint32_t ei) const;
  void find_missing_PLCedges(std::vector<uint32_t>& me) const;
  bool splitMissingEdge(uint32_t me);
  uint32_t findEncroachingPoint(const PLCedge& e, uint64_t& tet) const;

  void edgeSplit(const uint32_t ei, pointType* Pt_c, uint32_t acute_v_id);
  void middleEdgeSplit(const uint32_t ei);
  void splitStrategy1(const uint32_t ei, const uint32_t ref);
  void splitStrategy2(const uint32_t ei, const uint32_t ref);
  void segmentRecovery_HSi(bool quiet =false); // Segment recovery main function


  void makePLCfaces();
  void initFaceFlatEdges(PLCface& f);
  bool faceRecovery(bool quiet =false); // Face recovery main function

  // Exact predicates
  bool segmentCrossesFlatEdge(uint32_t ev[2], const std::vector<std::pair<uint32_t, uint32_t>>& flat_edges, int max_comp_normal);
  bool edgeIntersectsFacePlane(uint32_t v1, uint32_t v2, const PLCface& f);
  bool edgeIntersectsFace(uint32_t v1, uint32_t v2, const PLCface& f);
  bool lineIntersectsFace(uint32_t v1, uint32_t v2, const PLCface& f);
  bool innerEdgeIntersectsFace(uint32_t v1, uint32_t v2, const PLCface& f);
  bool triangleIntersectsFace(uint64_t t, const PLCface& f);
  bool tetIntersectsFace(uint64_t t, const PLCface& f);
  bool isTriangleOnFace(const uint32_t cv[3], uint32_t fi, const std::vector<std::pair<uint32_t, uint32_t>>& orig_flat_edges);
  bool tetIntersectsInnerTriangle(uint64_t t, uint32_t v1, uint32_t v2, uint32_t v3);

  // Collect tetrahedra whose interior intersects a PLC face.
  // If cornerMask is non-null, each tet face that overlaps with the PLC face is marked
  void getTetsIntersectingFace(uint32_t fi, std::vector<uint64_t> *i_tets, std::vector<bool> *cornerMask =NULL);

  // TRUE if v1 and v2 are consecutive in one of the boundary loops of f
  bool adjacentFaceVertices(uint32_t v1, uint32_t v2, const PLCface& f);

  bool isUpperCavityTet(const uint64_t t) const;
  bool isLowerCavityTet(const uint64_t t) const;

  int localOrient3d(uint32_t v1, uint32_t v2, uint32_t v3, uint32_t v4, std::vector<uint32_t>& to_unorient);
  int cachedOrient3D(uint32_t v, uint32_t v1, uint32_t v2, uint32_t v3);

  bool recoverFaceHSi(std::vector<uint64_t>& i_tets, const PLCface& f, bool& sisMethodWorks);


  void giftWrap(std::vector<uint64_t>& bnd, const std::vector<uint32_t>& vertices, std::vector<uint32_t>& newtets);
  uint64_t missingFaceInCavity(const std::vector<uint64_t>& bnd, const std::vector<uint32_t>& vertices);
  uint64_t meshCavity(const std::vector<uint64_t>& bnd, const std::vector<uint32_t>& vertices, std::vector<uint64_t>& base);
  uint64_t expandCavity(std::vector<uint64_t>& bnd, std::vector<uint32_t>& vertices, uint64_t t, const PLCface& f);

  size_t markInnerTets();

  //void getTetsIntersectingFaceSlow(uint32_t fi, std::vector<uint64_t>* i_tets) {
  //    const PLCface& f = faces[fi];

  //    //
  //    //// SLOW VERSION - USE TO CHECK
  //    for (uint32_t v : f.vertices) delmesh.marked_vertex[v] = 1;
  //    for (size_t i = 0; i < v_orient.size(); i++) v_orient[i] = UNDET_ORIENTATION;

  //    for (size_t i = 0; i < delmesh.numTets(); i++)
  //        if (tetIntersectsFace(i, f))
  //            i_tets->push_back(i);

  //    for (size_t i = 0; i < v_orient.size(); i++) v_orient[i] = UNDET_ORIENTATION;
  //    for (uint32_t v : f.vertices) delmesh.marked_vertex[v] = 0;
  //    return;
  //}

  //void saveFaces() const {
  //    FILE* fp = fopen("faces.off", "w");
  //    fprintf(fp, "OFF\n%lu %lu 0\n", delmesh.vertices.size(), faces.size());
  //    for (auto* v : delmesh.vertices) {
  //        double x, y, z;
  //        v->getApproxXYZCoordinates(x, y, z);
  //        fprintf(fp, "%f %f %f\n", x, y, z);
  //    }
  //    for (const PLCface& f : faces) {
  //        fprintf(fp, "%lu ", f.vertices.size());
  //        for (auto vi : f.vertices) fprintf(fp, "%u ", vi);
  //        fprintf(fp, "\n");
  //    }
  //    fclose(fp);
  //}
};

#endif // _PLC_
