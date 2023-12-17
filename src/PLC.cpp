#include "PLC.h"
#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <set>

#ifdef USE_INDIRECT_PREDS
double PLCx::getT1(uint32_t oe0i, uint32_t e0i) const {
    const std::vector<pointType*>& vs = delmesh.vertices;
    if (vs[e0i] == vs[oe0i]) return 0.0;
    else if (&vs[e0i]->toLNC().P() == vs[oe0i]) return vs[e0i]->toLNC().T();
    else return 1.0 - vs[e0i]->toLNC().T();
}

double PLCx::getT2(uint32_t oe1i, uint32_t e1i) const {
    const std::vector<pointType*>& vs = delmesh.vertices;
    if (vs[e1i] == vs[oe1i]) return 1.0;
    else if (&vs[e1i]->toLNC().Q() == vs[oe1i]) return vs[e1i]->toLNC().T();
    else return 1.0 - vs[e1i]->toLNC().T();
}

inline implicitPoint_LNC* PLCx::getProjectionOrMidPoint(uint32_t oe0i, uint32_t oe1i, uint32_t e0i, uint32_t e1i, uint32_t ri, uint32_t& acute_v) const
{
    const std::vector<pointType*>& vs = delmesh.vertices;
    const vector3d e_ref(vs[ri]);
    const vector3d e_oe0(vs[oe0i]);
    const vector3d e_oe1(vs[oe1i]);

    double discr = sqrt(e_ref.dist_sq(e_oe0) / ((e_oe1.dist_sq(e_oe0))));
    const double t1 = getT1(oe0i, e0i);
    const double t2 = getT2(oe1i, e1i);
    const double dv = (t2 - t1) * 0.2;
    if (discr <= (t1 + dv) || discr >= (t2 - dv)) { discr = (t1 + t2) / 2; acute_v = UINT32_MAX; }
    else acute_v = oe0i;

    return new implicitPoint_LNC(vs[oe0i]->toExplicit3D(), vs[oe1i]->toExplicit3D(), discr);
}

inline implicitPoint_LNC* PLCx::getProjectionOrMidPoint_noac(uint32_t oe0i, uint32_t oe1i, uint32_t e0i, uint32_t e1i, uint32_t ri, uint32_t& acute_v) const
{
    const std::vector<pointType*>& vs = delmesh.vertices;
    const vector3d e_ref(vs[ri]);
    const vector3d e_e0(vs[e0i]);
    const vector3d e_oe0(vs[oe0i]);
    const vector3d e_oe1(vs[oe1i]);

    double discr = sqrt(e_ref.dist_sq(e_e0) / ((e_oe1.dist_sq(e_oe0))));
    const double t1 = getT1(oe0i, e0i);
    const double t2 = getT2(oe1i, e1i);
    discr += t1;
    if (discr >= t2) { discr = (t1 + t2) / 2; acute_v = UINT32_MAX; }
    else acute_v = oe0i;

    return new implicitPoint_LNC(vs[oe0i]->toExplicit3D(), vs[oe1i]->toExplicit3D(), discr);
}

inline implicitPoint_LNC* PLCx::getProjectionOrMidPoint_noac_rev(uint32_t oe0i, uint32_t oe1i, uint32_t e0i, uint32_t e1i, uint32_t ri, uint32_t& acute_v) const
{
    const std::vector<pointType*>& vs = delmesh.vertices;
    const vector3d e_ref(vs[ri]);
    const vector3d e_e1(vs[e1i]);
    const vector3d e_oe0(vs[oe0i]);
    const vector3d e_oe1(vs[oe1i]);

    double discr = sqrt(e_ref.dist_sq(e_e1) / ((e_oe1.dist_sq(e_oe0))));
    const double t1 = getT1(oe0i, e0i);
    const double t2 = getT2(oe1i, e1i);
    discr = t2 - discr;
    if (discr <= t1) { discr = (t1 + t2) / 2; acute_v = UINT32_MAX; }
    else acute_v = oe1i;

    return new implicitPoint_LNC(vs[oe0i]->toExplicit3D(), vs[oe1i]->toExplicit3D(), discr);
}

inline implicitPoint_LNC* PLCx::getMidPoint(uint32_t oe0i, uint32_t oe1i, uint32_t e0i, uint32_t e1i) const
{
    const std::vector<pointType*>& vs = delmesh.vertices;
    const double t1 = getT1(oe0i, e0i);
    const double t2 = getT2(oe1i, e1i);
    return new implicitPoint_LNC(vs[oe0i]->toExplicit3D(), vs[oe1i]->toExplicit3D(), (t1 + t2) / 2);
}
#else
implicitPoint_LNC* PLCx::getProjectionOrMidPoint(uint32_t oe0i, uint32_t oe1i, uint32_t e0i, uint32_t e1i, uint32_t ri, uint32_t& acute_v) const
{
    const std::vector<pointType*>& vs = delmesh.vertices;
    const vector3d e_ref(vs[ri]);
    const vector3d e_e0(vs[e0i]);
    const vector3d e_e1(vs[e1i]);
    const vector3d e_oe0(vs[oe0i]);
    const vector3d e_oe1(vs[oe1i]);
    const coord_t elen = e_oe1.dist_sq(e_oe0);
    
    coord_t discr = GET_SQRT(coord_t(e_ref.dist_sq(e_oe0) / elen));
    const coord_t t1 = GET_SQRT(coord_t(e_e0.dist_sq(e_oe0) / elen));
    const coord_t t2 = GET_SQRT(coord_t(e_e1.dist_sq(e_oe0) / elen));
    const coord_t dv = (t2 - t1) * 0.2;
    if (discr <= (t1 + dv) || discr >= (t2 - dv)) { discr = (t1 + t2) / 2; acute_v = UINT32_MAX; }
    else acute_v = oe0i;

    return new pointType(*vs[oe0i], *vs[oe1i], discr);
}

implicitPoint_LNC* PLCx::getProjectionOrMidPoint_noac(uint32_t oe0i, uint32_t oe1i, uint32_t e0i, uint32_t e1i, uint32_t ri, uint32_t& acute_v) const
{
    const std::vector<pointType*>& vs = delmesh.vertices;
    const vector3d e_ref(vs[ri]);
    const vector3d e_e0(vs[e0i]);
    const vector3d e_e1(vs[e1i]);
    const vector3d e_oe0(vs[oe0i]);
    const vector3d e_oe1(vs[oe1i]);
    const coord_t elen = e_oe1.dist_sq(e_oe0);

    coord_t discr = GET_SQRT(coord_t(e_ref.dist_sq(e_e0) / elen));
    const coord_t t1 = GET_SQRT(coord_t(e_e0.dist_sq(e_oe0) / elen));
    const coord_t t2 = GET_SQRT(coord_t(e_e1.dist_sq(e_oe0) / elen));
    discr += t1;
    if (discr >= t2) { discr = (t1 + t2) / 2; acute_v = UINT32_MAX; }
    else acute_v = oe0i;

    return new pointType(*vs[oe0i], *vs[oe1i], discr);
}

implicitPoint_LNC* PLCx::getProjectionOrMidPoint_noac_rev(uint32_t oe0i, uint32_t oe1i, uint32_t e0i, uint32_t e1i, uint32_t ri, uint32_t& acute_v) const
{
    const std::vector<pointType*>& vs = delmesh.vertices;
    const vector3d e_ref(vs[ri]);
    const vector3d e_e0(vs[e0i]);
    const vector3d e_e1(vs[e1i]);
    const vector3d e_oe0(vs[oe0i]);
    const vector3d e_oe1(vs[oe1i]);
    const coord_t elen = e_oe1.dist_sq(e_oe0);

    coord_t discr = GET_SQRT(coord_t(e_ref.dist_sq(e_e1) / elen));
    const coord_t t1 = GET_SQRT(coord_t(e_e0.dist_sq(e_oe0) / elen));
    const coord_t t2 = GET_SQRT(coord_t(e_e1.dist_sq(e_oe0) / elen));
    discr = t2 - discr;
    if (discr <= t1) { discr = (t1 + t2) / 2; acute_v = UINT32_MAX; }
    else acute_v = oe1i;

    return new pointType(*vs[oe0i], *vs[oe1i], discr);
}

implicitPoint_LNC* PLCx::getMidPoint(uint32_t oe0i, uint32_t oe1i, uint32_t e0i, uint32_t e1i) const
{
    const std::vector<pointType*>& vs = delmesh.vertices;
    const vector3d e_e0(vs[e0i]);
    const vector3d e_e1(vs[e1i]);
    const vector3d e_oe0(vs[oe0i]);
    const vector3d e_oe1(vs[oe1i]);
    const coord_t elen = e_oe1.dist_sq(e_oe0);
    const coord_t t1 = GET_SQRT(coord_t(e_e0.dist_sq(e_oe0) / elen));
    const coord_t t2 = GET_SQRT(coord_t(e_e1.dist_sq(e_oe0) / elen));
    const coord_t discr = (t1 + t2) / 2;
    return new pointType(*vs[oe0i], *vs[oe1i], discr);
}
#endif

//void saveTETS(const char* filename, std::vector<uint32_t>& newtets, TetMesh& tin)
//{
//    ofstream f(filename);
//
//    if (!f) ip_error("\nTetMesh::saveTET: FATAL ERROR cannot open the file.\n");
//
//    f << tin.numVertices() << " vertices\n";
//
//    f << newtets.size() / 4 << " tets\n";
//    for (uint32_t i = 0; i < tin.numVertices(); i++)
//        f << *tin.vertices[i] << "\n";
//    const uint32_t* tet_node = newtets.data();
//    for (uint64_t i = 0; i < newtets.size() / 4; i++) {
//        f << "4 " << tet_node[i * 4] << " " << tet_node[i * 4 + 1] << " " << tet_node[i * 4 + 2] << " " << tet_node[i * 4 + 3] << "\n";
//    }
//
//    f.close();
//}

// TRUE if pqr is an acute angle at q
bool isAcuteAngle(const pointType* p, const pointType* q, const pointType* r) {
    return pointType::dotProductSign3D(*p, *r, *q) > 0;
}


// PLCface ---

void PLCface::zip() {

    while (!bounding_edges.empty() && bounding_edges.front() == bounding_edges.back()) {
        std::rotate(bounding_edges.begin(), bounding_edges.begin() + 1, bounding_edges.end());
        bounding_edges.pop_back();
        bounding_edges.pop_back();
    }

    bool iterate;
    do {
        iterate = false;
        for (size_t i = 1; i < bounding_edges.size(); i++) {
            const PLCedge* e1 = bounding_edges[i];
            const PLCedge* e2 = bounding_edges[i - 1];
            if (e1 == e2) {
                bounding_edges.erase(bounding_edges.begin() + i - 1, bounding_edges.begin() + i + 1);
                iterate = true;
            }
        }
    } while (iterate);

    auto ne = bounding_edges;
    std::sort(ne.begin(), ne.end());
    for (size_t i = 1; i < ne.size(); i++)
        if (ne[i - 1] == ne[i]) absorb(*this, ne[i]);
}

void PLCface::makeVertices() {
    size_t last_loop = 0;
    for (size_t i = 0; i < bounding_edges.size(); i++) {
        const PLCedge* e1 = bounding_edges[i];
        const PLCedge* e2 = bounding_edges[(i + 1) % bounding_edges.size()];
        const uint32_t cv = e1->commonVertex(*e2);
        if (cv != UINT32_MAX) vertices.push_back(cv);
        else {
            const uint32_t cv2 = e1->commonVertex(*bounding_edges[last_loop]);
            assert(cv2 != UINT32_MAX);
            vertices.push_back(cv2);
            last_loop = i + 1;
        }
    }
}

void PLCface::initConvexity(const PLCx& plc) {
    if (triangles.size() == 1) is_convex = true;    // A single triangle is always convex
    else {
        // Calculate a valid 2D projection plane for the face.
        const uint32_t t0 = triangles[0];
        const uint32_t tv[3] = { plc.input_tv[t0 * 3], plc.input_tv[t0 * 3 + 1], plc.input_tv[t0 * 3 + 2] };
        const TetMesh& dm = plc.delmesh;
        const std::vector<pointType*>& dmv = dm.vertices;
        const explicitPoint& tp1 = dmv[tv[0]]->toExplicit3D();
        const explicitPoint& tp2 = dmv[tv[1]]->toExplicit3D();
        const explicitPoint& tp3 = dmv[tv[2]]->toExplicit3D();

        max_comp_normal = pointType::maxComponentInTriangleNormal(tp1.X(), tp1.Y(), tp1.Z(), tp2.X(), tp2.Y(), tp2.Z(), tp3.X(), tp3.Y(), tp3.Z());

        for (size_t i = 1; i < bounding_edges.size(); i++) {
            const PLCedge* e1 = bounding_edges[i - 1];
            const PLCedge* e2 = bounding_edges[i];
            if (e1->commonVertex(*e2) == UINT32_MAX) {
                is_convex = false;
                is_simply_connected = false;
                return;
            }
        }

        is_simply_connected = true;

        // Check whether the projected 2D orientation of consecutive 
        // triplets of original vertices has always the same sign.
        bool neg = false, pos = false;
        for (size_t i = 0; i < bounding_edges.size(); i++) {
            const PLCedge* e1 = bounding_edges[i];
            size_t i2;
            for (i2 = (i + 1) % bounding_edges.size(); i2 < bounding_edges.size(); i2++)
                if (!bounding_edges[i2 % bounding_edges.size()]->hasOriginalVertices(e1->oep[0], e1->oep[1])) break;
            const PLCedge* e2 = bounding_edges[i2 % bounding_edges.size()];
            const uint32_t v2 = e1->commonOriginalVertex(*e2);
            const uint32_t v1 = e1->oppositeOriginalVertex(v2);
            const uint32_t v3 = e2->oppositeOriginalVertex(v2);
            const pointType* c1 = dm.vertices[v1];
            const pointType* c2 = dm.vertices[v2];
            const pointType* c3 = dm.vertices[v3];
            int co;
            if (max_comp_normal == 0) co = pointType::orient2Dyz(*c1, *c2, *c3);
            else if (max_comp_normal == 1) co = pointType::orient2Dzx(*c1, *c2, *c3);
            else co = pointType::orient2Dxy(*c1, *c2, *c3);
            if (co < 0) neg = true;
            if (co > 0) pos = true;
        }

        is_convex = (!(neg && pos));
    }
}

//--------------
// PLCx
//--------------

bool PLCx::faceHasTriangle(const PLCface& f, const uint32_t tv[3]) const {
    for (const uint32_t t : f.triangles) {
        const uint32_t* ftv = input_tv + t * 3;
        if ((ftv[0] == tv[0] || ftv[0] == tv[1] || ftv[0] == tv[2]) &&
            (ftv[1] == tv[0] || ftv[1] == tv[1] || ftv[1] == tv[2]) &&
            (ftv[2] == tv[0] || ftv[2] == tv[1] || ftv[2] == tv[2])) return true;
    }
    return false;
}

// For each face, for each of its vertices, set of incident face triangles
void PLCx::makeVertexTriangleMaps(std::vector<std::vector<std::vector<uint32_t>>>& vt_maps) {
    vt_maps.clear();
    vt_maps.resize(faces.size());

    v_reindex.resize(delmesh.vertices.size(), UINT32_MAX);
    std::vector<bool> orig_tri_mark(input_nt, false);
    for (size_t i = 0; i < faces.size(); i++)
        makeVertexTriangleMap(faces[i], vt_maps[i], orig_tri_mark);
}

void PLCx::makeVertexTriangleMap(PLCface& f, std::vector<std::vector<uint32_t>>& vt_map,
    std::vector<bool>& orig_tri_mark) {

    // Assume that vertices are unmarked, and use marks to keep track
    // of non-manifold vertices that may appear several times in f.vertices

    // Explicit map (first = vertex id, second = position in f)
    for (uint32_t v : f.vertices) delmesh.marked_vertex[v] = 0;
    for (uint32_t v : f.vertices) delmesh.marked_vertex[v]++;
    for (uint32_t i = 0; i < (uint32_t)f.vertices.size(); i++) {
        const uint32_t v = f.vertices[i];
        if (delmesh.marked_vertex[v] > 1)
            singular_v.push_back(std::pair<uint32_t, uint32_t>(v, i));
        v_reindex[v] = i;
    }

    for (uint32_t i : f.triangles) orig_tri_mark[i] = true;

    // ve-map is only used to retrieve the vt for steiner points.
    // we may omit the non manifold case check here.
    std::vector<std::vector<const PLCedge*>> ve_map;
    ve_map.resize(f.vertices.size());
    for (size_t i = 0; i < f.bounding_edges.size(); i++) {
        const PLCedge* e = f.bounding_edges[i];
        assert(e->type != flat);
        ve_map[v_reindex[e->ep[0]]].push_back(e);
        ve_map[v_reindex[e->ep[1]]].push_back(e);
    }

    vt_map.resize(f.vertices.size());
    for (size_t vi = 0; vi < f.vertices.size(); vi++)
        if (isSteinerVertex(f.vertices[vi])) {
            assert(delmesh.marked_vertex[f.vertices[vi]] == 1);
            assert(ve_map[vi].size() == 2);
            const PLCedge* e = ve_map[vi][0];
            for (uint32_t t : e->inc_tri)
                if (orig_tri_mark[t])
                    vt_map[vi].push_back(t);
        }

    for (uint32_t i : f.triangles)
        for (uint32_t j = 0; j < 3; j++) {
            const uint32_t v = input_tv[i * 3 + j];
            if (v_reindex[v] != UINT32_MAX) {
                if (delmesh.marked_vertex[v] == 1) vt_map[v_reindex[v]].push_back(i);
                else {
                    for (auto p : singular_v) if (p.first == v) vt_map[p.second].push_back(i);
                }
            }

        }

    for (uint32_t i : f.triangles) orig_tri_mark[i] = false;

    for (uint32_t v : f.vertices) {
        v_reindex[v] = UINT32_MAX;
        delmesh.marked_vertex[v] = 0;
    }

    singular_v.clear();
}

// Removes duplicated pre-PLCedges
void PLCx::mergePreEdges(){

  // sort edges by lexicografic non-descending endpoints
  std::sort(edges.begin(), edges.end(), PLCedge::vertexSortFunc);

  for(uint32_t ei=0; ei < edges.size()-1; ){
    PLCedge& e = edges[ei];

    while( ++ei < edges.size() && e.ep[0]==edges[ei].ep[0] && e.ep[1] == edges[ei].ep[1]){
      PLCedge& ne = edges[ei];
      // each edge initially has only one incident triangle, ne never been merged
      e.inc_tri.push_back( ne.inc_tri[0] );
      ne.inc_tri.clear(); // ne is now isolated
    }

  }

  // remove duplicated edges (no inc_faces) from edges vectcor
  edges.erase(std::remove_if(edges.begin(), edges.end(), PLCedge::isIsolatedPtr), edges.end());
}

// Assumes that the PLCedge e is one of the sides of the input-triangle ti.
// Returns the index of the vertex of ti different from e endpoints.
uint32_t PLCx::opposite_vrt(const PLCedge& e, const uint32_t ti) const {
  uint32_t v = input_tv[3*ti];
  if(v != e.ep[0]  &&  v != e.ep[1] ) return v;

  v = input_tv[3*ti + 1];
  if(v != e.ep[0]  &&  v != e.ep[1] ) return v;

  v = input_tv[3*ti + 2];
  if(v != e.ep[0]  &&  v != e.ep[1] ) return v;

  return UINT32_MAX;
}

// Returns true if at least a couple of non-flat edges incident at the vertex
// indexed as vi forms an acute angle, i.e. the scalar produc between these two
// edges is positive.
bool PLCx::isAcute(const uint32_t vi, const std::vector<std::vector<uint32_t>>& vv) const {
    const pointType* vip = delmesh.vertices[vi];
    for (uint32_t i = 0; i < vv[vi].size(); i++)
        for (uint32_t j = 0; j < i; j++)
            if (isAcuteAngle(delmesh.vertices[vv[vi][i]], vip, delmesh.vertices[vv[vi][j]])) return true;

  return false;
}

// Fill the PLCx data structure by using the input triangulation information
// ASSUMPTIONS: all faces are triangles.
void PLCx::initialize(){
  std::vector<bool> is_acute(delmesh.numVertices(), false);

  // -- Fill edges --
  edges.resize(input_nt * 3); // More than necessary, will be merged later

  // Fill only EV and partial ET relations
  for(uint32_t ti=0; ti<input_nt; ti++){
    const uint32_t ei = 3*ti;
    const uint32_t v0 = input_tv[ei];
    const uint32_t v1 = input_tv[ei+1];
    const uint32_t v2 = input_tv[ei+2];
    edges[ei].fill_preEdge(v0,v1,ti);
    edges[ei+1].fill_preEdge(v1,v2,ti);
    edges[ei+2].fill_preEdge(v2,v0,ti);
  }

  // merge duplicated pre-edges
  mergePreEdges(); // Now edges contains proper edges.

  // Mark flat edges and fill vv (only for non-flat edges)
  std::vector<std::vector<uint32_t>> vv; // vv relation relative to PLC
  vv.resize(delmesh.vertices.size());

  is_polyhedron = true;
  for(uint32_t ei=0; ei<edges.size(); ei++){
    PLCedge& e = edges[ei];
    if (isFlat(e)) {
        e.type = flat;
    }
    else{
      e.type = undet;
      vv[e.ep[0]].push_back(e.ep[1]);
      vv[e.ep[1]].push_back(e.ep[0]);
      // Not a polyhedron if an edge has an odd number of incident triangles
      if (e.inc_tri.size() & 1) is_polyhedron = false;
    }
  }

  // Mark acute vertices
  is_acute.resize(delmesh.vertices.size());
  for(uint32_t vi=0; vi<delmesh.vertices.size(); vi++)
      is_acute[vi] = isAcute(vi, vv);

  // Classify non-flat edges by type
  for(uint32_t ei=0; ei<edges.size(); ei++) if(edges[ei].type != flat){
    PLCedge& e = edges[ei];
    const bool e0_acute = is_acute[ e.ep[0] ];
    const bool e1_acute = is_acute[ e.ep[1] ];
    if(e0_acute && e1_acute) e.type = both_acute_ep;
    else if(e0_acute || e1_acute){
      e.type = one_acute_ep;
      if(e1_acute) e.swap();
    }
    else e.type = no_acute_ep;
  }

  //size_t last = edges.size() - 1;
  //while (last && edges[last].type != flat) last--;
  //for (size_t i = 0; i < last; i++) if (edges[i].type != flat) {
  //    std::swap(edges[i], edges[last--]);
  //    while (last && edges[last].type != flat) last--;
  //}

  delmesh.marked_vertex.resize(delmesh.vertices.size());
}


// -- Segment recovery algorithm --

// Returns the index of the vertex between those that encroach upon
// PLCedge e = <ep0,ep1> that forms with ep0 and ep1 the maximum circle.
// If there are no encroahing vrts UINT32_MAX is returned.
// Note that if e is missing at least one encroaching point must exist.
// DEF. of encroaching point: a vertex that lies inside or on the diameter
//                            sphere of segment <ep0,ep1>.
// In the case that exist more than one encroaching point, we choose the one
// such that "the circumradius of the smallest circumsphere of the triangle
// <ep0,ep1,enc_pt> is maximum over other encroaching points."
// Since for 3 points exists an unique circle C, the "smallest circumsphere of
// the triangle <ep0,ep1,enc_pt>" is the one having C as equatorial circle.
//

uint32_t PLCx::findEncroachingPoint(const PLCedge& e, uint64_t& tet) const {
    return delmesh.findEncroachingPoint(e.ep[0], e.ep[1], tet);
}

// Return true if edges[ei] is a missing edge (i.e. non-flat PLCedge that is not
// a side of some tetrahedron of delmesh)
bool PLCx::is_missing_PLCedge(const uint32_t ei) const{
    return !delmesh.hasEdge(edges[ei].ep[0], edges[ei].ep[1]);
}

// Find all missing edges (i.e. non-flat PLCedges that are not a side of some
// tetrahedron of delmesh) and add their indices to vector me.
void PLCx::find_missing_PLCedges(std::vector<uint32_t>& me) const {
    for (uint32_t ei = 0; ei < edges.size(); ei++) {
        const PLCedge& e = edges[ei];
        if (e.type != flat && is_missing_PLCedge(ei)) me.push_back(ei);
  }
}


// Splits PLCedge e, inserting the point Pt (which is internal to the edge),
// inserts Pt in the vector vertices of TetMesh DS (do not update tetrahedrization),
// updates PLCx DS.
// This split function works only for edges of types: "no_acute_ep", "one_acute_ep"
void PLCx::edgeSplit(const uint32_t ei, pointType* Pt_c, uint32_t acute_v_id){
  PLCedge& e = edges[ei];

  // 1-Create new vertex
  const uint32_t Pt_i = (uint32_t)delmesh.vertices.size(); // New vertex (mid point of e) index.
  // Add the point Pt to TetMesh vertices vector
  pushVertex(Pt_c, acute_v_id);

  // 2-Updates PLCedges
  const uint32_t e1 = e.ep[1]; // Memory parent edge endpoints
  e.ep[1] = Pt_i; // Update PLCedge e endpoints: <e0,e1> becomes <e0,Pt>
  edges.push_back(PLCedge(Pt_i, e1, e.oep[0], e.oep[1], e.inc_tri, e.type));
}

// Splits PLCedge e, inserting its middle point Mpt,
// inserts Mpt in the vector vertices of TetMesh DS (do not update tetrahedrization),
// updates PLCx DS.
// This split function works only for edges of types: "no_acute_ep", "one_acute_ep", "both_acute_ep"
void PLCx::middleEdgeSplit(const uint32_t ei){
  PLCedge& e = edges[ei];
  const uint32_t e0 = e.ep[0], e1 = e.ep[1]; // Memory partent edge endpoints
  const uint32_t oe0 = e.oep[0], oe1 = e.oep[1]; // Memory partent edge endpoints

  implicitPoint_LNC* np = getMidPoint(oe0, oe1, e0, e1);

  // 1-Create new vertex
  const uint32_t Mpt_i = (uint32_t)delmesh.vertices.size(); // New vertex (mid point of e) index.
  // Add the mid point of e to TetMesh vertices vector
  pushVertex(np, UINT32_MAX);

  // 2-Updates PLCedges
  e.ep[1] = Mpt_i; // Update PLCedge e endpoints: <e0,e1> becomes <e0,Mpt>

  if (e.type == both_acute_ep) {
      e.type = one_acute_ep;
      edges.push_back(PLCedge(e1, Mpt_i, e.oep[1], e.oep[0], e.inc_tri, one_acute_ep));
  }
  else edges.push_back(PLCedge(Mpt_i, e1, e.oep[0], e.oep[1], e.inc_tri, e.type));
}

// e = <e0,e1> is of type "no_acute_ep"
void PLCx::splitStrategy1(const uint32_t ei, const uint32_t ref){
  const PLCedge& e = edges[ei];

  const pointType* e0p = delmesh.vertices[e.ep[0]];
  const pointType* e1p = delmesh.vertices[e.ep[1]];
  const pointType* refp = delmesh.vertices[ref];
  implicitPoint_LNC *vc_p;

  uint32_t acute_v = UINT32_MAX;

  if (vector3d::isAtMostTwiceDistanceThan(e0p, refp, e1p))
      vc_p = getProjectionOrMidPoint_noac(e.oep[0], e.oep[1], e.ep[0], e.ep[1], ref, acute_v);
  else if (vector3d::isAtMostTwiceDistanceThan(e1p, refp, e0p))
      vc_p = getProjectionOrMidPoint_noac_rev(e.oep[0], e.oep[1], e.ep[0], e.ep[1], ref, acute_v);
  else vc_p = getMidPoint(e.oep[0], e.oep[1], e.ep[0], e.ep[1]);

  edgeSplit(ei, vc_p, acute_v);
}

//  e = <e0,e1> is of type "one_acute_ep"
void PLCx::splitStrategy2(const uint32_t ei, const uint32_t ref){
  const PLCedge& e = edges[ei];
  implicitPoint_LNC* vc_p;

  uint32_t acute_v;
  vc_p = getProjectionOrMidPoint(e.oep[0], e.oep[1], e.ep[0], e.ep[1], ref, acute_v);

  if (vector3d::isCloserThan(vc_p, delmesh.vertices[e.ep[1]], delmesh.vertices[ref])) {
      delete vc_p;
      // Here we should switch to splitStrategy3, but it is not really necessary
      // Just using midpoint provides better performances
      vc_p = getMidPoint(e.oep[0], e.oep[1], e.ep[0], e.ep[1]);
  }

  edgeSplit(ei, vc_p, acute_v);
}

//
bool PLCx::splitMissingEdge(uint32_t mei) {
    if (!is_missing_PLCedge(mei)) return false;

    uint64_t ct;
    const PLCedge& e = edges[mei];

    // 1-SPLIT the missing edge with a Steiner point
    if (e.type == both_acute_ep) {
        ct = delmesh.inc_tet[e.ep[0]] << 2;
        middleEdgeSplit(mei);
    }
    else {
        const uint32_t ref = delmesh.findEncroachingPoint(e.ep[0], e.ep[1], ct);
        if (ref == UINT32_MAX) { // This should never happen
            printf("WARNING: Could not find a valid encroaching point! Using midpoint...\n");
            ct = delmesh.inc_tet[e.ep[0]] << 2;
            middleEdgeSplit(mei);
        }
        else if (e.type == no_acute_ep) splitStrategy1(mei, ref);
        else splitStrategy2(mei, ref);
    }

    // 2-UPDATE TetMesh by tetrahedralizing around steiner_point (last elem of vertices)
    delmesh.insertExistingVertex((uint32_t)delmesh.vertices.size() - 1, ct);

    return true;
}

int myrand(void)
{
    static int h = 1;
    return(((h = h * 214013L + 2531011L) >> 16) & 0x7fff);
}

template< class T > void shuffle_vec(T first, T last)
{
    for (auto i = (last - first)-1; i > 0; --i) {
        std::swap(first[i], first[myrand() % (i + 1)]);
    }
}

//// 
void PLCx::segmentRecovery_HSi(bool quiet)
{
    for (uint64_t tet_i = 0; tet_i < delmesh.numTets(); tet_i++)
        delmesh.unmark_Tet_1(tet_i);

    // indices of missing PLCedges (flat edges are not considered)
    std::vector<uint32_t> miss_edges;
    find_missing_PLCedges(miss_edges);

    // If edge was not touched by the insertion, no need to check if it is missing
    std::vector<bool> touched_vertex(delmesh.vertices.size(), false);
    std::vector<uint32_t> loval_vv;
    size_t num_st = 0;

    // Main loop
    while (!miss_edges.empty()) {
        const size_t nme = miss_edges.size();
        while (!miss_edges.empty()) {
            const uint32_t ei = miss_edges.back();
            const uint32_t vm = delmesh.numVertices();
            miss_edges.pop_back();
            const uint32_t v1 = edges[ei].ep[0];
            const uint32_t v2 = edges[ei].ep[1];
            if (splitMissingEdge(ei)) {
                loval_vv.push_back(v1);
                loval_vv.push_back(v2);
                delmesh.VV(vm, loval_vv);
                for (uint32_t w : loval_vv) touched_vertex[w] = true;
                touched_vertex.push_back(true);
                loval_vv.clear();

                num_st++;
                if (!quiet && (num_st & 1023) == 1023) {
                    printf("\rMissing edges: %zu (Tot %zu Steiner points inserted)                 ", nme, num_st); fflush(stdout);
                }
            }
        }

        for (uint32_t ei = 0; ei < edges.size(); ei++) {
            const PLCedge& e = edges[ei];
            if (e.type != flat && (touched_vertex[e.ep[0]] && touched_vertex[e.ep[1]]) && is_missing_PLCedge(ei)) miss_edges.push_back(ei);
        }
        for (uint32_t i = 0; i < delmesh.numVertices(); i++) touched_vertex[i] = 0;
        
        shuffle_vec(miss_edges.begin(), miss_edges.end());

        if (!quiet) {
            printf("\rMissing edges: %zu (Tot %zu Steiner points inserted)                 ", nme, num_st); fflush(stdout);
        }
    }
    if (!quiet) printf("\n");

    delmesh.removeDelTets();
}



class iEdge {
public:
    uint32_t v1, v2; // Endpoints
    std::vector<uint32_t> pEdges; // Ordered chain of PLC (sub-)edges
    std::vector<uint32_t> iTris; // incident triangles

    iEdge(uint32_t _v1, uint32_t _v2, uint32_t t) : v1(_v1), v2(_v2), iTris{ t } { if (v1 > v2) std::swap(v1, v2); }
    bool operator<(const iEdge& e) const { return (2 * ((v1 < e.v1) - (v1 > e.v1)) + ((v2 < e.v2) - (v2 > e.v2))) > 0; }
    bool operator==(const iEdge& e) const { return v1 == e.v1 && v2 == e.v2; }
    void absorb(iEdge& e) { iTris.insert(iTris.end(), e.iTris.begin(), e.iTris.end()); e.iTris.clear(); }

    static bool isEmpty(const iEdge& e) { return e.iTris.empty(); }

    void fillEdges(const std::vector<uint32_t>& te, std::vector<PLCedge>& edges) {
        for (size_t i = 0; i < te.size(); i++) {
            PLCedge& e = edges[te[i]];
            if (e.oep[0] == v2 && e.oep[1] == v1) e.swap();
            if (e.oep[0] == v1 && e.oep[1] == v2) pEdges.push_back(te[i]);
        }

        // This is a dumb O(n^2) sorting, but should be good enough for small vectors
        uint32_t pv = v1;
        for (size_t i = 0; i < pEdges.size(); i++) {
            for (size_t j = i; j < pEdges.size(); j++) {
                PLCedge& e = edges[pEdges[j]];
                if (e.ep[0] == pv) {
                    pv = e.ep[1];
                    std::swap(pEdges[i], pEdges[j]);
                    break;
                }
            }
        }

        assert(pv == v2);
    }

    void copyOrderedEdges(PLCface& f, uint32_t tv1, uint32_t tv2, uint32_t tv3, std::vector<PLCedge>& edges) const {
        std::vector<PLCedge*> tie;
        tie.resize(pEdges.size());
        for (size_t i = 0; i < pEdges.size(); i++) tie[i] = &edges[pEdges[i]];

        if (tv1 == v1 && tv2 == v2) { // Insert at the beginning in same order
            f.bounding_edges.insert(f.bounding_edges.begin(), tie.begin(), tie.end());
        } else if (tv1 == v2 && tv2 == v1) { // Insert at the beginning in reverse order
            std::reverse(tie.begin(), tie.end());
            f.bounding_edges.insert(f.bounding_edges.begin(), tie.begin(), tie.end());
        }
        else if (tv2 == v1 && tv3 == v2) { // Insert in the middle in same order
            size_t pos = 0;
            for (; pos < f.bounding_edges.size(); pos++) {
                const PLCedge *e = f.bounding_edges[pos];
                if (e->ep[0] == tv3 || e->ep[1] == tv3) break;
            }
            f.bounding_edges.insert(f.bounding_edges.begin()+pos, tie.begin(), tie.end());
        }
        else if (tv2 == v2 && tv3 == v1) { // Insert in the middle in reverse order
            size_t pos = 0;
            for (; pos < f.bounding_edges.size(); pos++) {
                const PLCedge* e = f.bounding_edges[pos];
                if (e->ep[0] == tv3 || e->ep[1] == tv3) break;
            }
            std::reverse(tie.begin(), tie.end());
            f.bounding_edges.insert(f.bounding_edges.begin() + pos, tie.begin(), tie.end());
        }
        else if (tv3 == v1 && tv1 == v2) { // Insert at the end in same order
            f.bounding_edges.insert(f.bounding_edges.end(), tie.begin(), tie.end());
        }
        else if (tv3 == v2 && tv1 == v1) { // Insert at the end in reverse order
            std::reverse(tie.begin(), tie.end());
            f.bounding_edges.insert(f.bounding_edges.end(), tie.begin(), tie.end());
        }
        else ip_error("Should never get here!\n");
    }
};

void PLCface::absorb(PLCface& f, PLCedge* e) {
    std::vector<PLCedge*> nb; // new bounding edges
    auto& be = bounding_edges;
    auto& ne = f.bounding_edges;

    size_t pos = 0, pos2 = 0;

    pos = std::find(bounding_edges.begin(), bounding_edges.end(), e) - bounding_edges.begin();
    if (pos == bounding_edges.size()) return;

    if (&f == this) {
        pos2 = pos + 1; // Make sure to get the other instance
        while (be[pos2] != e) pos2++; // e's position in ne
        nb.assign(be.begin() + pos + 1, be.begin() + pos2);
        be.erase(be.begin() + pos, be.begin() + pos2 + 1);
        be.insert(be.end(), nb.begin(), nb.end());
    }
    else {
        while (ne[pos2] != e) pos2++; // e's position in ne

        // If faces have opposite orientations -> reverse
        size_t npos = (pos + 1) % be.size();
        size_t npos2 = (pos2 + 1) % ne.size();
        if (be[npos]->commonVertex(*e) == ne[npos2]->commonVertex(*e)) {
            std::reverse(ne.begin(), ne.end());
            pos2 = (ne.size() - 1) - pos2;
        }

        if (pos) nb.insert(nb.end(), be.begin(), be.begin() + pos);
        if (pos2 != (ne.size() - 1)) nb.insert(nb.end(), ne.begin() + pos2 + 1, ne.end());
        if (pos2) nb.insert(nb.end(), ne.begin(), ne.begin() + pos2);
        if (pos != (be.size() - 1)) nb.insert(nb.end(), be.begin() + pos + 1, be.end());
        be.assign(nb.begin(), nb.end());
        triangles.insert(triangles.end(), f.triangles.begin(), f.triangles.end());

        ne.clear();
        f.triangles.clear();
    }
}

void PLCx::makePLCfaces() {
    if (!faces.empty()) return; // Faces were previously built

    // Create original mesh halfedges
    std::vector<iEdge> iEdges;
    for (uint32_t i = 0; i < input_nt; i++) {
        iEdges.push_back(iEdge(input_tv[i * 3], input_tv[i * 3 + 1], i));
        iEdges.push_back(iEdge(input_tv[i * 3 + 1], input_tv[i * 3 + 2], i));
        iEdges.push_back(iEdge(input_tv[i * 3 + 2], input_tv[i * 3], i));
    }

    // Transform halfedges into full edges
    std::sort(iEdges.begin(), iEdges.end());
    iEdge* prev = &iEdges[0];
    for (size_t i = 1; i < iEdges.size(); i++) 
        if (iEdges[i] == *prev) prev->absorb(iEdges[i]);
        else prev = &iEdges[i];
    iEdges.erase(std::remove_if(iEdges.begin(), iEdges.end(), iEdge::isEmpty), iEdges.end());

    // for each plcEdge e, for each inc_tri, associate e to inc_tri in a map-vector triEdges
    std::vector<std::vector<uint32_t>> triEdges;
    triEdges.resize(input_nt);
    for (size_t i = 0; i < edges.size(); i++)
        for (uint32_t& t : edges[i].inc_tri) triEdges[t].push_back((uint32_t)i);

    // for each iEdge e, consider first iTri and copy submap from triEdges (only edges that have v1-v2 as orig endpoints)
    // sort pEdges consistently with v1-v2
    for (iEdge& e : iEdges) e.fillEdges(triEdges[e.iTris[0]], edges);
    

    // Create a face for each input triangle
    faces.resize(input_nt);
    for (uint32_t i = 0; i < input_nt; i++) faces[i].triangles.push_back(i);

    // For each iEdge e, copy its ordered edge chain to each of its incident tris in faces
    for (iEdge& e : iEdges)
        for (uint32_t t : e.iTris)
            e.copyOrderedEdges(faces[t], input_tv[t * 3], input_tv[t * 3 + 1], input_tv[t * 3 + 2], edges);

    // At this point, each face has its ordered chain of bounding edges

    // Merge faces across flat mesh edges
    std::vector<uint32_t> remap;
    remap.resize(input_nt);
    for (uint32_t i = 0; i < input_nt; i++) remap[i] = i;

    for (iEdge& e : iEdges) if (edges[e.pEdges.front()].type == flat) {
        uint32_t t1i = e.iTris[0];
        uint32_t t2i = e.iTris[1];
        while (remap[t1i] != t1i) t1i = remap[t1i];
        while (remap[t2i] != t2i) t2i = remap[t2i];
        if (t1i == t2i) continue;
        if (t2i < t1i) std::swap(t1i, t2i);
        faces[t1i].absorb(faces[t2i], &edges[e.pEdges[0]]);
        remap[t2i] = t1i;
    }

    for (PLCface& f : faces) f.zip();

    // This latest cycle makes non simply-connected faces
    for (iEdge& e : iEdges) if (edges[e.pEdges.front()].type == flat) {
        uint32_t t1i = e.iTris[0];
        uint32_t t2i = e.iTris[1];
        while (remap[t1i] != t1i) t1i = remap[t1i];
        while (remap[t2i] != t2i) t2i = remap[t2i];
        if (t1i == t2i && !PLCface::isEmpty(faces[t1i])) faces[t1i].absorb(faces[t2i], &edges[e.pEdges[0]]);
    }

    faces.erase(std::remove_if(faces.begin(), faces.end(), PLCface::isEmpty), faces.end());

    // Make face vertices
    for (PLCface& f : faces) f.makeVertices();

    // Init flat vertices
    for (size_t i = 0; i < delmesh.vertices.size(); i++) delmesh.marked_vertex[i] = 0;
    for (PLCface& f : faces) {
        for (PLCedge* e : f.bounding_edges)
            delmesh.marked_vertex[e->oep[0]] = delmesh.marked_vertex[e->oep[1]] = 1;

        for (uint32_t t : f.triangles)
            for (int i = 0; i < 3; i++) {
                const uint32_t v = input_tv[t * 3 + i];
                if (!delmesh.marked_vertex[v]) {
                    delmesh.marked_vertex[v] = 1;
                    f.flat_vertices.push_back(v);
                }
            }
        for (uint32_t t : f.triangles)
            for (int i = 0; i < 3; i++)
                delmesh.marked_vertex[input_tv[t * 3 + i]] = 0;
    }

    for (PLCface& f : faces) f.initConvexity(*this);
    for (PLCface& f : faces) initFaceFlatEdges(f);

    // For each face, for each of its vertices, set of incident face triangles
    makeVertexTriangleMaps(vt_maps);
    v_reindex.resize(delmesh.vertices.size());
}



bool PLCx::edgeIntersectsFacePlane(uint32_t v1, uint32_t v2, const PLCface& f) {
    const uint32_t* tv = input_tv + f.triangles[0] * 3;
    const int o1 = cachedOrient3D(v1, tv[0], tv[1], tv[2]);
    const int o2 = cachedOrient3D(v2, tv[0], tv[1], tv[2]);
    return (o1 >= 0 && o2 <= 0) || (o1 <= 0 && o2 >= 0);
}

bool PLCx::edgeIntersectsFace(uint32_t v1, uint32_t v2, const PLCface& f) {
    const uint32_t* fv = input_tv + f.triangles[0]*3;
    const int o1 = cachedOrient3D(v1, fv[0], fv[1], fv[2]);
    const int o2 = cachedOrient3D(v2, fv[0], fv[1], fv[2]);
    if ((o1 > 0 && o2 > 0) || (o1 < 0 && o2 < 0)) return false;
    else return lineIntersectsFace(v1, v2, f);
}

bool PLCx::innerEdgeIntersectsFace(uint32_t v1, uint32_t v2, const PLCface& f) {
    const uint32_t* fv = input_tv + f.triangles[0] * 3;
    const int o1 = cachedOrient3D(v1, fv[0], fv[1], fv[2]);
    const int o2 = cachedOrient3D(v2, fv[0], fv[1], fv[2]);
    if ((o1 >= 0 && o2 >= 0) || (o1 <= 0 && o2 <= 0)) return false;
    else return lineIntersectsFace(v1, v2, f);
}

bool PLCx::lineIntersectsFace(uint32_t v1, uint32_t v2, const PLCface& f) {
    const pointType* vp[2] = { delmesh.vertices[v1], delmesh.vertices[v2] };
    for (uint32_t i : f.triangles) {
        const uint32_t t[3] = { input_tv[i * 3], input_tv[i * 3 + 1], input_tv[i * 3 + 2] };
        const pointType* tv[3] = { delmesh.vertices[t[0]], delmesh.vertices[t[1]], delmesh.vertices[t[2]] };
        if (pointType::lineCrossesTriangle(*vp[0], *vp[1], *tv[0], *tv[1], *tv[2])) return true;
    }

    return false;
}


bool PLCx::triangleIntersectsFace(uint64_t t, const PLCface& f) {
    uint64_t tb = t & (~3);
    uint32_t tv[3];
    for (int j = 0; tb < tb + 4; tb++)
        if (tb != t) tv[j++] = delmesh.tet_node[tb];

    return (
        edgeIntersectsFace(tv[0], tv[1], f) ||
        edgeIntersectsFace(tv[1], tv[2], f) ||
        edgeIntersectsFace(tv[2], tv[0], f)
        );
}

bool PLCx::tetIntersectsFace(uint64_t t, const PLCface& f) {
    const uint32_t* n = delmesh.tet_node.data() + (t<<2);
    if (n[3] == INFINITE_VERTEX) return false;

    for (int i = 0; i < 4; i++)
        for (int j = i + 1; j < 4; j++)
            if (!delmesh.marked_vertex[n[i]] && !delmesh.marked_vertex[n[j]] && innerEdgeIntersectsFace(n[i], n[j], f))
                return true;

    return false;
}

bool PLCx::tetIntersectsInnerTriangle(uint64_t t, uint32_t v1, uint32_t v2, uint32_t v3) {
    const uint32_t* n = delmesh.tet_node.data() + (t << 2);
    if (n[3] == INFINITE_VERTEX) return false;

    for (int i = 0; i < 4; i++)
        for (int j = i + 1; j < 4; j++)
            if (!delmesh.marked_vertex[n[i]] && !delmesh.marked_vertex[n[j]] && 
                (cachedOrient3D(n[i], v1, v2, v3) * cachedOrient3D(n[j], v1, v2, v3)) <= 0 &&
                pointType::lineCrossesTriangle(*delmesh.vertices[n[i]], *delmesh.vertices[n[j]], 
                *delmesh.vertices[v1], *delmesh.vertices[v2], *delmesh.vertices[v3]))
                    return true;

    return false;
}


bool PLCx::adjacentFaceVertices(uint32_t v1, uint32_t v2, const PLCface& f) {
    if (!f.is_simply_connected || delmesh.marked_vertex[v1]>1 || delmesh.marked_vertex[v2]>1) {
        for (PLCedge* e : f.bounding_edges)
            if ((e->ep[0] == v1 && e->ep[1] == v2) || (e->ep[0] == v2 && e->ep[1] == v1))
                return true;
        return false;
    }
    uint32_t idx[2] = { v_reindex[v1],  v_reindex[v2] };
    if (idx[1] < idx[0]) std::swap(idx[0], idx[1]);
    const uint32_t id = idx[1] - idx[0];
    return (id == 1 || id == (f.vertices.size()-1));
}

int PLCx::localOrient3d(uint32_t v1, uint32_t v2, uint32_t v3, uint32_t v4, std::vector<uint32_t>& to_unorient) {
    if (v_orient[v1] != UNDET_ORIENTATION) return v_orient[v1];
    to_unorient.push_back(v1);
    return (v_orient[v1] = delmesh.vOrient3D(v1, v2, v3, v4));
}

int PLCx::cachedOrient3D(uint32_t v, uint32_t v1, uint32_t v2, uint32_t v3) {
    if (v_orient[v] == UNDET_ORIENTATION)
        v_orient[v] = delmesh.vOrient3D(v, v1, v2, v3);
    return v_orient[v];
}

inline void pushAndMark(uint64_t t, TetMesh& m, std::vector<uint64_t>& B) {
    B.push_back(t);
    m.mark_Tet_1(t);
}

void PLCx::getTetsIntersectingFace(uint32_t fi, std::vector<uint64_t> *i_tets, std::vector<bool>* cornerMask) {
    const PLCface& f = faces[fi];

    // Let e=(v1, v2) be a nonflat edge in f
    uint32_t v_t[4];
    const PLCedge* e0 = f.bounding_edges.front();
    v_t[0] = e0->ep[0]; v_t[1] = e0->ep[1];

    std::vector<uint64_t> et;
    delmesh.ET(v_t[0], v_t[1], et);

    const uint32_t *tv = input_tv + f.triangles[0] * 3; // Vertices of one face triangle for orientation

    // If face has only three vertices, just check whether they are in one tet of the ET
    if (f.vertices.size() == 3 && f.flat_vertices.empty()) {
        const uint32_t ov = (tv[0] != v_t[0] && tv[0] != v_t[1]) ? (tv[0]) : ((tv[1] != v_t[0] && tv[1] != v_t[1]) ? (tv[1]) : (tv[2]));
        for (uint64_t t : et) if (delmesh.tetHasVertex(t, ov)) {
            if (cornerMask) {
                const uint64_t c = delmesh.tetOppositeCorner(t, v_t[0], v_t[1], ov);
                (*cornerMask)[c] = (*cornerMask)[delmesh.tet_neigh[c]] = true;
            }
            return;
        }
    }

    // Mark f vertices and init orientations
    for (uint32_t v : f.vertices) {
        delmesh.marked_vertex[v]++;
        v_orient[v] = 0;
    }
    for (uint32_t v : f.flat_vertices) v_orient[v] = 0;

    // init v_reindex with f 

    for (uint32_t i = 0; i < (uint32_t)f.vertices.size(); i++) {
        const uint32_t v = f.vertices[i];
        if (delmesh.marked_vertex[v] > 1)
            singular_v.push_back(std::pair<uint32_t, uint32_t>(v, i));
        v_reindex[v] = i;
    }

    // Flat orig edges in f
    const std::vector<std::pair<uint32_t, uint32_t>>& orig_flat_edges = f.orig_flat_edges;

    std::vector<uint32_t> to_unorient;

    //  Find a tet t0 in ET(e) intersecting the face interior
    uint64_t t0 = UINT64_MAX;
    for (uint64_t t : et) {
        //    - let v3 and v4 be the vertices of t0 opposite wrt e (oppositeTetEdge)
        delmesh.oppositeTetEdge(t<<2, v_t, v_t + 2);

        // If we are using this function to mark faces we need three vertices on the plane
        if (cornerMask) {
            if (v_orient[v_t[2]] == 0 && isTriangleOnFace(v_t, fi, orig_flat_edges)) { t0 = t; break; }
            std::swap(v_t[2], v_t[3]);
            if (v_orient[v_t[2]] == 0 && isTriangleOnFace(v_t, fi, orig_flat_edges)) { t0 = t; break; }
            continue;
        }

        //    - if v_orient[v3]!=0 e v_orient[v3]==v_orient[v4] -> skip (no intersection)
        const int ov3 = localOrient3d(v_t[2], tv[0], tv[1], tv[2], to_unorient);
        const int ov4 = localOrient3d(v_t[3], tv[0], tv[1], tv[2], to_unorient);
        if (ov3 != 0 && ov3 == ov4) continue;

        //const uint64_t tb = t << 2;
        //const uint32_t* n = delmesh.tet_node.data() + tb;

        //    - if v_orient[v3]*v_orient[v4]<0 e lineIntersectsFace(v3,v4,f) -> intersection
        if (ov3 * ov4 < 0 && lineIntersectsFace(v_t[2], v_t[3], f)) { t0 = t; break; }
        //    - if v_orient[v3]==0
        if (ov3 == 0) {
            //      - if v3 not marked -> skip
            if (!delmesh.marked_vertex[v_t[2]]) continue;
            //      - if v1,v2,v3 part of f
            if (isTriangleOnFace(v_t, fi, orig_flat_edges)) { t0 = t; break; }
        }

        //    - if v_orient[v4]==0
        if (ov4 == 0) {
            //      - if v4 not marked -> skip
            if (!delmesh.marked_vertex[v_t[3]]) continue;
            //      - if v1,v2,v4 part of f
            std::swap(v_t[2], v_t[3]);
            if (isTriangleOnFace(v_t, fi, orig_flat_edges)) { t0 = t; break; }
        }
    }

    // B = empty
    std::vector<uint64_t> B;

    // Mark t0 and insert in B
    if (t0!=UINT64_MAX) B.push_back(t0);
    if (f.flat_vertices.size()) {
        for (uint32_t v : f.flat_vertices) delmesh.VT(v, B);
        B.erase(std::unique(B.begin(), B.end()), B.end());
        for (uint64_t t : B) delmesh.mark_Tet_1(t);
    }
    else delmesh.mark_Tet_1(t0);

    // In the remainder, OK means "add n to B, mark it"
    // for each tet t in B
    //int ct = 0;
    for (size_t k = 0; k < B.size(); k++) {
        const uint64_t t = B[k];
        const uint64_t tb = t << 2;
        const uint64_t* nn = delmesh.tet_neigh.data() + tb;
        const uint32_t* tn = delmesh.tet_node.data() + tb;
        // for each of the four neighbors n of t
        for (int i = 0; i < 4; i++) {
            const uint64_t n = nn[i] >> 2;
            //   if n is not a ghost and is not marked
            if (!delmesh.isGhost(n) && !delmesh.is_marked_Tet_1(n)) {
                //      let cv be the three common vertices of t and n
                const uint32_t cv[3] = { tn[(i + 1) & 3], tn[(i + 2) & 3], tn[(i + 3) & 3] };
                const unsigned char mv[3] = { delmesh.marked_vertex[cv[0]], delmesh.marked_vertex[cv[1]], delmesh.marked_vertex[cv[2]] };
                const int o3d[3] = {
                    localOrient3d(cv[0], tv[0], tv[1], tv[2], to_unorient),
                    localOrient3d(cv[1], tv[0], tv[1], tv[2], to_unorient),
                    localOrient3d(cv[2], tv[0], tv[1], tv[2], to_unorient)
                };

                if (mv[0] && mv[1] && mv[2]) pushAndMark(n, delmesh, B);
                else if (mv[0] && mv[1]) {
                    if (!adjacentFaceVertices(cv[0], cv[1], f)) pushAndMark(n, delmesh, B);
                }
                else if (mv[1] && mv[2]) {
                    if (!adjacentFaceVertices(cv[1], cv[2], f)) pushAndMark(n, delmesh, B);
                }
                else if (mv[2] && mv[0]) {
                    if (!adjacentFaceVertices(cv[2], cv[0], f)) pushAndMark(n, delmesh, B);
                }
                else if (mv[0]) {
                    if (o3d[1] * o3d[2] < 0) pushAndMark(n, delmesh, B);
                }
                else if (mv[1]) {
                    if (o3d[2] * o3d[0] < 0) pushAndMark(n, delmesh, B);
                }
                else if (mv[2]) {
                    if (o3d[0] * o3d[1] < 0) pushAndMark(n, delmesh, B);
                }
                else {
                    if (o3d[1] * o3d[2] < 0 ||
                        o3d[2] * o3d[0] < 0 ||
                        o3d[0] * o3d[1] < 0) pushAndMark(n, delmesh, B);
                }
            }
        }
    }

    for (uint64_t t : B) {
        delmesh.unmark_Tet_1(t);
        const uint32_t* n = delmesh.tet_node.data() + (t << 2);
        const int uv[4] = { localOrient3d(n[0], tv[0], tv[1], tv[2], to_unorient),
                            localOrient3d(n[1], tv[0], tv[1], tv[2], to_unorient),
                            localOrient3d(n[2], tv[0], tv[1], tv[2], to_unorient),
                            localOrient3d(n[3], tv[0], tv[1], tv[2], to_unorient) };

        if (cornerMask) {
            int j = 0, nj;
            for (int i = 0; i < 4; i++) if (uv[i] == 0) j++; else nj = i;
            if (j == 3) {
                const uint64_t c = (t << 2) + nj;
                (*cornerMask)[c] = (*cornerMask)[delmesh.tet_neigh[c]] = true;
            }
        }
        else if (!((uv[0] >= 0 && uv[1] >= 0 && uv[2] >= 0 && uv[3] >= 0) || (uv[0] <= 0 && uv[1] <= 0 && uv[2] <= 0 && uv[3] <= 0))) {
            assert(tetIntersectsFace(t, f));
            i_tets->push_back(t);
        }
    }

    for (uint32_t v : to_unorient) v_orient[v] = UNDET_ORIENTATION;

    for (uint32_t v : f.vertices) {
        delmesh.marked_vertex[v] = 0;
        v_orient[v] = UNDET_ORIENTATION;
    }
    for (uint32_t v : f.flat_vertices) v_orient[v] = UNDET_ORIENTATION;
    singular_v.clear();

    //if (i_tets.size()) {
    //    savePLC("face.off", &f);
    //    save_tet_faces("to_remove.off", i_tets);
    //    save_tet_faces("fat.off", B);
    //    std::vector<uint64_t> one(1, t0);
    //    save_tet_faces("first.off", one);
    //}


    //
    //// SLOW VERSION - USE TO CHECK
    ////for (uint32_t v : f.vertices) delmesh.marked_vertex[v] = 1;
    ////for (size_t i = 0; i < v_orient.size(); i++) v_orient[i] = UNDET_ORIENTATION;

    ////for (size_t i = 0; i < delmesh.numTets(); i++)
    ////    if (tetIntersectsFace(i, f))
    ////        i_tets.push_back(i);

    ////for (size_t i = 0; i < v_orient.size(); i++) v_orient[i] = UNDET_ORIENTATION;
    ////for (uint32_t v : f.vertices) delmesh.marked_vertex[v] = 0;
    ////return;
}


bool PLCx::segmentCrossesFlatEdge(uint32_t ev[2], const std::vector<std::pair<uint32_t, uint32_t>>& flat_edges, int max_comp_normal) {
    const pointType& A = *delmesh.vertices[ev[0]];
    const pointType& B = *delmesh.vertices[ev[1]];
    for (auto& e : flat_edges) {
        if ((e.first == ev[0] && e.second == ev[1]) || (e.first == ev[1] && e.second == ev[0])) return false;
        if (e.first == ev[0] || e.second == ev[1] || e.first == ev[1] || e.second == ev[0]) continue;
        const pointType& P = *delmesh.vertices[e.first];
        const pointType& Q = *delmesh.vertices[e.second];

        if (pointType::innerSegmentsCross(A, B, P, Q, max_comp_normal)) return true;
    }
    return false;
}

// A triangle cv = {v1,v2,v3} is part of the face f=faces[fi] if:
//  v1,v2,v3 are all marked (meaning they are on the border of f) AND
//   f is convex OR 
//   v1,v2,v3 is one of its triangles OR
//   v1,v2,v3 is a subtriangle of one of the triangles of f OR
//   one of the sides vi,vi+1 "crosses" one of the flat edges of f
//
// This function assumes that none of the cv is a flat vertex !

bool PLCx::isTriangleOnFace(const uint32_t cv[3], uint32_t fi, const std::vector<std::pair<uint32_t, uint32_t>>& orig_flat_edges) {
    if (delmesh.marked_vertex[cv[0]] && delmesh.marked_vertex[cv[1]] && delmesh.marked_vertex[cv[2]]) {
        const PLCface& f = faces[fi];
        if (f.is_convex || faceHasTriangle(f, cv)) return true;

        for (auto j : f.triangles) {
            bool hascv[3];
            for (int i = 0; i < 3; i++) {
                uint32_t v = v_reindex[cv[i]];
                const auto& vts = vt_maps[fi][v];
                hascv[i] = (std::find(vts.begin(), vts.end(), j) != vts.end());
            }
            if (hascv[0] && hascv[1] && hascv[2]) {
                return true;
            }
        }

        for (int i = 0; i < 3; i++) {
            uint32_t ev[2] = { cv[i], cv[(i + 1) % 3] };
            if (segmentCrossesFlatEdge(ev, orig_flat_edges, f.max_comp_normal)) return true;
        }

        if (f.flat_vertices.size()) {
            for (int i = 0; i < 3; i++)
                if (std::find(f.flat_vertices.begin(), f.flat_vertices.end(), cv[i]) != f.flat_vertices.end())
                    return true;

            for (int i = 0; i < 3; i++) {
                uint32_t ev[2] = { cv[i], cv[(i + 1) % 3] };
                for (uint32_t v : f.flat_vertices)
                    if (pointType::pointInInnerSegment(*delmesh.vertices[v], *delmesh.vertices[ev[0]], *delmesh.vertices[ev[1]])) return true;
            }
        }
    }
    
    return false;
}


void PLCx::initFaceFlatEdges(PLCface& f) {
    std::set<std::pair<uint32_t, uint32_t>> orig_edges;
    for (uint32_t t : f.triangles) {
        const uint32_t* tfv = input_tv + t * 3;
        for (int i = 0; i < 3; i++) {
            uint32_t ev[2] = { tfv[i], tfv[(i + 1) % 3] };
            if (ev[1] < ev[0]) std::swap(ev[0], ev[1]);
            if (!orig_edges.insert(std::pair<uint32_t, uint32_t>(ev[0], ev[1])).second)
                f.orig_flat_edges.push_back(std::pair<uint32_t, uint32_t>(ev[0], ev[1]));
        }
    }
}


// Mark internal tetrahedra
size_t PLCx::markInnerTets() {

    // If the PLC does not define a valid polyhedron, just tag every tet as IN but the ghosts
    if (!is_polyhedron) {
        size_t ng = 0;
        for (size_t i = 0; i < delmesh.numTets(); i++) 
            if (delmesh.isGhost(i)) delmesh.mark_tetrahedra[i] = DT_OUT;
            else {
                delmesh.mark_tetrahedra[i] = DT_IN;
                ng++;
            }
        return ng;
    }

    // 1) Mark constraint corners
    std::vector<bool> cornerMask(delmesh.tet_node.size(), false);



    // Crea relazione VF
    //   per ogni faccia f aggiungi f alla VF di tutti i suoi bounding e internal vertices
    // Per ogni triangolo in delmesh, cerca la(le) faccia comune f in VF(v1), VF(v2) e VF(v3)
    //   se c'è più di una faccia comune marca immediatamente il triangolo e passa oltre
    //   altrimenti
    // scopri se il triangolo sta dentro o fuori dalla faccia comune f (orient2d?)
    //   1) se f è convessa
    //   2) se v1 è interno a f (o v2, o v3)
    //   3) se il baricentro del triangolo è interno a uno dei triangoli di f e al triangolo stesso (check per possibile errore numerico)
    // Se la faccia comune esiste e il triangolo ci sta dentro allora marcala, altrimenti no


    for (size_t fi = 0; fi < faces.size(); fi++)
        getTetsIntersectingFace((uint32_t)fi, NULL, &cornerMask);

    // --- DEBUG
    //int num_constrained_faces = 0;
    //for (uint64_t t = 0; t < delmesh.tet_node.size(); t++) 
    //    if (!delmesh.isGhost(t>>2) && cornerMask[t] && 
    //        (delmesh.isGhost(delmesh.tet_neigh[t]>>2) || t>delmesh.tet_neigh[t])) num_constrained_faces++;

    //ofstream f("constraints.off");
    //if (!f) ip_error("PLCx::savePLC: Cannot open file.\n");

    //f << "OFF\n";
    //f << delmesh.numVertices() << " " << num_constrained_faces << " 0\n";

    //for (uint32_t i = 0; i < delmesh.numVertices(); i++)
    //    f << *delmesh.vertices[i] << "\n";

    //for (uint64_t t = 0; t < delmesh.tet_node.size(); t++) 
    //    if (!delmesh.isGhost(t >> 2) && cornerMask[t] &&
    //        (delmesh.isGhost(delmesh.tet_neigh[t]>>2) || t > delmesh.tet_neigh[t])) {
    //    uint32_t v[3];
    //    delmesh.getFaceVertices(t, v);
    //    f << "3 " << v[0] << " " << v[1] << " " << v[2] << "\n";
    //}

    //f.close();
    // --------

    return delmesh.markInnerTets(cornerMask);
}

bool PLCx::faceRecovery(bool quiet) {
    makePLCfaces();

    v_orient.resize(delmesh.vertices.size(), UNDET_ORIENTATION);

    delmesh.marked_vertex.assign(delmesh.vertices.size(), 0);
    for (uint64_t tet_i = 0; tet_i < delmesh.numTets(); tet_i++) delmesh.unmark_Tet_1(tet_i);

    bool needRecursion, sisMethodWorks = true;
    do {
        needRecursion = false;

        for (size_t i = 0; i < faces.size(); i++) {
            const PLCface& f = faces[i];
            const uint32_t* tv = input_tv + f.triangles[0] * 3; // The vertices of f for orientation
            std::vector<uint64_t> i_tets;
            getTetsIntersectingFace((uint32_t)i, &i_tets);

            if (i_tets.size()) {
                if (!quiet) printf("\rRecovering face: %zu                  ", i); fflush(stdout);
                std::vector<uint32_t> cavity_v;
                for (uint32_t v : f.vertices) v_orient[v] = 0;
                for (uint64_t t : i_tets) {
                    const uint32_t* v = delmesh.getTetNodes(t << 2);
                    for (int j = 0; j < 4; j++)
                        if (v_orient[v[j]] == UNDET_ORIENTATION) {
                            v_orient[v[j]] = delmesh.vOrient3D(v[j], tv[0], tv[1], tv[2]);
                            cavity_v.push_back(v[j]);
                        }
                }

                //delmesh.recoverFaceGiftWrap(i_tets, v_orient);
                //sisMethodWorks = false;

                if (!recoverFaceHSi(i_tets, f, sisMethodWorks)) needRecursion = true;

                i_tets.clear();

                // Reset orientations for the future
                for (uint32_t v : cavity_v) v_orient[v] = UNDET_ORIENTATION;
                for (uint32_t v : f.vertices) v_orient[v] = UNDET_ORIENTATION;
            }
        }
        if (!quiet) printf("\n");
        if (needRecursion && !quiet) printf("RECURSION\n");
    } while (needRecursion);

    return sisMethodWorks;
}

bool PLCx::isUpperCavityTet(const uint64_t t) const {
    uint32_t v[3];
    delmesh.getFaceVertices(t, v);
    return v_orient[v[0]] >= 0 && v_orient[v[1]] >= 0 && v_orient[v[2]] >= 0;
}

bool PLCx::isLowerCavityTet(const uint64_t t) const {
    uint32_t v[3];
    delmesh.getFaceVertices(t, v);
    return v_orient[v[0]] <= 0 && v_orient[v[1]] <= 0 && v_orient[v[2]] <= 0;
}


bool PLCx::recoverFaceHSi(std::vector<uint64_t>& i_tets, const PLCface& f, bool& sisMethodWorks) {

    //std::vector<uint32_t> oldtets;
    //for (uint64_t t : i_tets) {
    //    const uint32_t* tn = delmesh.tet_node.data() + (t * 4);
    //    for (int i = 0; i < 4; i++) oldtets.push_back(tn[i]);
    //}
    //saveTETS("old_tets.tet", oldtets, delmesh);
 
    // Create vector 'top_faces' and 'bottom_faces'
    std::vector<uint64_t> top_faces, bottom_faces;
    for (uint64_t t : i_tets) delmesh.mark_Tet_1(t);

    // Move all tets to remove to tail
    uint32_t last = delmesh.numTets() - 1;
    for (uint64_t& t : i_tets) {
        while (delmesh.is_marked_Tet_1(last)) last--;
        if (t < last) {
            delmesh.swapTets(t, last);
            t = last;
        }
    }
    while (delmesh.is_marked_Tet_1(last)) last--;

    last++; // This is now the numtets after having removed toremove

    for (uint64_t t : i_tets) {
        const uint64_t* neigh = delmesh.getTetNeighs(t << 2);
        for (int i = 0; i < 4; i++) if (!delmesh.is_marked_Tet_1(neigh[i] >> 2)) {
            if (isUpperCavityTet(neigh[i])) top_faces.push_back(neigh[i]);
            else {
                assert(isLowerCavityTet(neigh[i])); // Must be either
                bottom_faces.push_back(neigh[i]);
            }
        }
    }

    for (uint64_t t : i_tets) delmesh.unmark_Tet_1(t);

    std::vector<uint32_t> top_vertices, bottom_vertices;

    for (uint64_t t : i_tets) {
        for (int i = 0; i < 4; i++) {
            const uint32_t v = delmesh.tet_node[t * 4 + i];
            if (!delmesh.marked_vertex[v]) {
                delmesh.marked_vertex[v] = 1;
                if (v_orient[v] >= 0) top_vertices.push_back(v);
                if (v_orient[v] <= 0) bottom_vertices.push_back(v);
            }
        }
    }
    for (uint32_t w : f.flat_vertices) if (!delmesh.marked_vertex[w]) {
        delmesh.marked_vertex[w] = 1; 
        top_vertices.push_back(w);
        bottom_vertices.push_back(w);
    }

    std::sort(top_vertices.begin(), top_vertices.end());
    std::sort(bottom_vertices.begin(), bottom_vertices.end());
    for (uint32_t w : top_vertices) delmesh.marked_vertex[w] = 0;
    for (uint32_t w : bottom_vertices) delmesh.marked_vertex[w] = 0;

    uint64_t t;
    std::vector<uint64_t> toremove;
    bool cavity_ok = true;

    while ((t = missingFaceInCavity(top_faces, top_vertices))!=UINT64_MAX) {
        toremove.push_back(expandCavity(top_faces, top_vertices, t, f));
        if (v_orient[top_vertices.back()] < 0) { // This might happen due to a theoretical bug in Si's paper
            cavity_ok = false; break;
        }
        std::sort(top_vertices.begin(), top_vertices.end());
    }
    while ((t = missingFaceInCavity(bottom_faces, bottom_vertices)) != UINT64_MAX) {
        toremove.push_back(expandCavity(bottom_faces, bottom_vertices, t, f));
        if (v_orient[bottom_vertices.back()] > 0) { // This might happen due to a theoretical bug in Si's paper
            cavity_ok = false; break;
        }
        std::sort(bottom_vertices.begin(), bottom_vertices.end());
    }

//    if (!cavity_ok) exit(10);

    if (cavity_ok) {
        delmesh.tet_node.resize(last * 4);
        delmesh.tet_neigh.resize(last * 4);
        delmesh.mark_tetrahedra.resize(last);

        t = meshCavity(top_faces, top_vertices, bottom_faces);
        assert(t == UINT64_MAX);
        t = meshCavity(bottom_faces, bottom_vertices, top_faces);
        assert(t == UINT64_MAX);
    }

    if (!toremove.empty()) {
        if (cavity_ok) {
            for (uint32_t v : top_vertices) v_orient[v] = UNDET_ORIENTATION;
            for (uint32_t v : bottom_vertices) v_orient[v] = UNDET_ORIENTATION;
            for (uint64_t y : toremove) delmesh.pushAndMarkDeletedTets(y);
            delmesh.removeDelTets();
            return false;
        }
        // If Hang Si's cavity expansion fails revert to slower gift-wrap method
        else {
            if (!f.flat_vertices.empty()) 
                ip_error("Hang Si's cavity expansion fails on a flat face with internal vertices.\nTreatment of this very particular case was not implemented!\n");
            delmesh.recoverFaceGiftWrap(i_tets, v_orient);
            sisMethodWorks = false;
            for (uint32_t v : top_vertices) v_orient[v] = UNDET_ORIENTATION;
            for (uint32_t v : bottom_vertices) v_orient[v] = UNDET_ORIENTATION;
        }
    }

    return true;
}

uint64_t PLCx::expandCavity(std::vector<uint64_t>& bnd, std::vector<uint32_t>& vertices, uint64_t t, const PLCface& f) {
    assert(!delmesh.isGhost(t>>2));

    // Remove 't' corner from bnd
    auto pos = std::find(bnd.begin(), bnd.end(), t);
    assert(pos != bnd.end());
    bnd.erase(pos);

    // For each of the other three corners y
    uint64_t tb = (t & (~3));
    int added = 0;
    for (uint64_t i = 0; i < 4; i++) {
        const uint64_t y = tb + i;
        if (y != t) {
            auto pos = std::find(bnd.begin(), bnd.end(), y);
            //   if y in bnd, remove from bnd
            //   otherwise add tet_neigh[y] to bnd
            if (pos != bnd.end()) bnd.erase(pos);
            else {
                bnd.push_back(delmesh.tet_neigh[y]);
                added++;
            }
        }
    }
    //
    // If three tet_neigh were added to bnd, add their common vertex to 'vertices'
    if (added == 3) {
        const uint32_t v = delmesh.tet_node[t];
        if (std::find(vertices.begin(), vertices.end(), v) == vertices.end()) {
            const uint32_t* tv = input_tv + f.triangles[0] * 3;
            if (v_orient[v] == UNDET_ORIENTATION) {
                v_orient[v] = delmesh.vOrient3D(v, tv[0], tv[1], tv[2]);
            }
            vertices.push_back(v);
        }
    }

    return tb;
}


class bdUpdater {
public:
    uint64_t t1, t2;
    uint64_t bnd;

    bdUpdater(uint64_t _t1, uint64_t _t2, uint64_t _bnd) : t1(_t1), t2(_t2), bnd(_bnd) {}
};

uint64_t PLCx::missingFaceInCavity(const std::vector<uint64_t>& bnd, const std::vector<uint32_t>& vertices) {
    // DT of vertices
    TetMesh dt(true);
    dt.vertices.resize(vertices.size());
    size_t i = 0;
    for (uint32_t v : vertices) dt.vertices[i++] = delmesh.vertices[v];
    dt.inc_tet.resize(vertices.size(), UINT64_MAX);
    dt.tetrahedrize();

    uint64_t nt[2];
    for (size_t i = 0; i < vertices.size(); i++) v_reindex[vertices[i]] = (uint32_t)i;
    for (uint64_t t : bnd) {
        uint32_t v[3];
        delmesh.getFaceVertices(t, v);
        for (int i = 0; i < 3; i++) v[i] = v_reindex[v[i]];

        if (!dt.getTetsFromFaceVertices(v[0], v[1], v[2], nt)) {
            for (size_t i = 0; i < vertices.size(); i++) v_reindex[vertices[i]] = UINT32_MAX;
            return t; // Need cavity expansion
        }
    }
    for (size_t i = 0; i < vertices.size(); i++) v_reindex[vertices[i]] = UINT32_MAX;

    return UINT64_MAX;
}

uint64_t PLCx::meshCavity(const std::vector<uint64_t>& bnd, const std::vector<uint32_t>& vertices, std::vector<uint64_t>& base) {

    // DT of vertices
    TetMesh dt(true);
    dt.vertices.resize(vertices.size());
    size_t i = 0;
    for (uint32_t v : vertices) dt.vertices[i++] = delmesh.vertices[v];
    dt.inc_tet.resize(vertices.size(), UINT64_MAX);
    dt.tetrahedrize();

    // Swap ghosts to tail
    uint32_t last = dt.numTets() - 1;
    while (dt.isGhost(last)) last--;
    for (uint64_t t = 0; t < last; t++)
        if (dt.isGhost(t)) {
            dt.swapTets(t, last--);
            while (dt.isGhost(last)) last--;
        }

    // Identify and mark bnd_faces in DT
    //   for each triplet v1,v2,v3, find in VT(v1) the two tets that share v2 and v3
    //   and mark their opposite corners.
    std::vector<bool> cornerMask(dt.tet_node.size(), false);
    uint64_t nt[2];
    std::vector<bdUpdater> bdpairs;

    for (size_t i = 0; i < vertices.size(); i++) v_reindex[vertices[i]] = (uint32_t)i;
    for (uint64_t t : bnd) {
        uint32_t v[3];
        delmesh.getFaceVertices(t, v);
        for (int i = 0; i < 3; i++) v[i] = v_reindex[v[i]];

        if (!dt.getTetsFromFaceVertices(v[0], v[1], v[2], nt)) {
            for (size_t i = 0; i < vertices.size(); i++) v_reindex[vertices[i]] = UINT32_MAX;
            return t; // Need cavity expansion
        }
        const uint64_t c1 = dt.getCornerFromOppositeTet(nt[0], nt[1]);
        const uint64_t c2 = dt.getCornerFromOppositeTet(nt[1], nt[0]);
        cornerMask[c1] = true;
        cornerMask[c2] = true;
        bdpairs.push_back(bdUpdater(c1, c2, t));
    }
    for (size_t i = 0; i < vertices.size(); i++) v_reindex[vertices[i]] = UINT32_MAX;

    // Mark all tets in DT as UNKNOWN
    for (i = 0; i < dt.numTets(); i++) dt.mark_tetrahedra[i] = DT_UNKNOWN; // Might be redundant... CHECK

    // Make vector C of classified tets
    std::vector<uint64_t> C;

    // Find a ghost tet t in DT that has at least a vertex with v_orient!=0
    for (i = 0; i < dt.tet_node.size(); i += 4)
        if (dt.tet_node[i + 3] == INFINITE_VERTEX &&
            (v_orient[vertices[dt.tet_node[i]]] || v_orient[vertices[dt.tet_node[i+1]]] || v_orient[vertices[dt.tet_node[i+2]]])
            ) break;
    assert(i < dt.tet_node.size());

    // Mark internal tets starting from i
    dt.markInnerTets(cornerMask, i >> 2);

    std::vector<uint64_t> remap;
    remap.resize(dt.numTets());
    for (size_t i = 0; i < remap.size(); i++) remap[i] = i;


    // Swap DT_OUT/ghost tets to tail
    last = dt.numTets() - 1;
    while (dt.mark_tetrahedra[last] == DT_OUT) last--;
    for (uint64_t t = 0; t < last; t++)
        if (dt.mark_tetrahedra[t] == DT_OUT) {
            std::swap(remap[last], remap[t]);
            dt.swapTets(t, last--);
            while (dt.mark_tetrahedra[last] == DT_OUT) last--;
        }

    // Transform DT_OUT tets to ghosts and make sure that all vt* point to DT_IN tets
    uint64_t num_in_tets = 0;
    for (size_t i = 0; i < dt.numTets(); i++) {
        if (dt.mark_tetrahedra[i] == DT_OUT) {
            dt.tet_node[i * 4 + 3] = INFINITE_VERTEX;
        }
        else {
            num_in_tets++;
            assert(dt.mark_tetrahedra[i] == DT_IN);
            dt.inc_tet[dt.tet_node[i * 4 + 0]] = dt.inc_tet[dt.tet_node[i * 4 + 1]] =
                dt.inc_tet[dt.tet_node[i * 4 + 2]] = dt.inc_tet[dt.tet_node[i * 4 + 3]] = i;
        }
    }

    // Use map to update bdpairs
    for (bdUpdater& b : bdpairs) {
        b.t1 = (remap[b.t1 >> 2] << 2) + (b.t1 & 3);
        b.t2 = (remap[b.t2 >> 2] << 2) + (b.t2 & 3);
    }

    // Here DT has its own connectivity and all ghosts are in tail.
    // Old i_tets have already been removed from delmesh.
    // Can now insert DT in delmesh using bnd to rebuild adjacencies.

    // 2) Remap dt.tet_neigh while considering their new positions in delmesh
    for (size_t i = 0; i < dt.tet_neigh.size(); i++) dt.tet_neigh[i] += delmesh.tet_neigh.size();

    // 3) Update inc_tets
    for (size_t i = 0; i < dt.vertices.size(); i++)
        dt.inc_tet[i] += delmesh.numTets();

    // 4) Remap dt.tet_node so that they point to IDs in delmesh 
    //    (if != INF_VERTEX dt.tet_node[i] = vertices[dt.tet_node[i]])
    //
    for (size_t i = 0; i < dt.tet_node.size(); i++)
        if (dt.tet_node[i] != INFINITE_VERTEX) dt.tet_node[i] = vertices[dt.tet_node[i]];

    size_t oNumTets = delmesh.numTets();
    num_in_tets += oNumTets;
    for (size_t i = 0; i < dt.numVertices(); i++) {
        if (dt.inc_tet[i] < num_in_tets)
            delmesh.inc_tet[vertices[i]] = dt.inc_tet[i];
    }

    delmesh.tet_node.insert(delmesh.tet_node.end(), dt.tet_node.begin(), dt.tet_node.end());
    delmesh.tet_neigh.insert(delmesh.tet_neigh.end(), dt.tet_neigh.begin(), dt.tet_neigh.end());
    delmesh.mark_tetrahedra.resize(delmesh.mark_tetrahedra.size() + dt.mark_tetrahedra.size(), 0);

    size_t on4 = oNumTets * 4;
    for (bdUpdater& b : bdpairs) {
        if (dt.mark_tetrahedra[b.t1 >> 2] == DT_IN) {
            delmesh.tet_neigh[on4 + b.t1] = b.bnd;
            delmesh.tet_neigh[b.bnd] = on4 + b.t1;
        }
        else {
            assert(dt.mark_tetrahedra[b.t2 >> 2] == DT_IN);
            delmesh.tet_neigh[on4 + b.t2] = b.bnd;
            delmesh.tet_neigh[b.bnd] = on4 + b.t2;
        }
    }

    // Delete mesh tail
    last = delmesh.numTets();
    while (delmesh.isGhost(--last));
    last++;
    delmesh.tet_node.resize(last*4);
    delmesh.tet_neigh.resize(last*4);
    delmesh.mark_tetrahedra.resize(last);

    for (size_t i = on4; i < delmesh.tet_node.size(); i++) {
        if (delmesh.tet_neigh[i] >= delmesh.tet_node.size()) {
            base.push_back(i);
            uint32_t w[3];
            delmesh.getFaceVertices(i, w);
            assert(v_orient[w[0]] == 0 && v_orient[w[1]] == 0 && v_orient[w[2]] == 0);
        }
    }

    return UINT64_MAX;
}
