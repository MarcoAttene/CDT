#include "PLC.h"
#include <iostream>
#include <fstream>
#include <random>
#include <algorithm>
#include <set>


// ------------------------- //
// Fill cavity gift-wrapping //
// ------------------------- //

inline int invOrient3D(const pointType& t0, const pointType& t1, const pointType& t2, const pointType& t3) {
    return -pointType::orient3D(t0, t1, t2, t3);
}

inline int invInSphere(const pointType& t0, const pointType& t1, const pointType& t2, const pointType& t3, const pointType& t4) {
    return -pointType::inSphere(t0, t1, t2, t3, t4);
}

inline bool chkInnerSegmentsCross(const pointType& u1, const pointType& u2, const pointType& v1, const pointType& v2) {
    return (pointType::orient3D(u1, u2, v1, v2) == 0) && pointType::innerSegmentsCross(u1, u2, v1, v2);
}

inline bool same_triangle(const uint32_t* u, const uint32_t* v) {
    if (u[0] != v[0] && u[0] != v[1] && u[0] != v[2]) return false;
    if (u[1] != v[0] && u[1] != v[1] && u[1] != v[2]) return false;
    if (u[2] != v[0] && u[2] != v[1] && u[2] != v[2]) return false;
    return true;
}


// T = <t0,t1,t2,t3> is a well defined tetrahedron
// returns true if v is in the interior of T, false otherwise.
inline bool is_vrt_inside_valid_tet(const pointType& v, const pointType& t0, const pointType& t1, const pointType& t2, const pointType& t3) {
    const int f = invOrient3D(t0, t1, t2, v); 
    return (f != 0 && invOrient3D(t0, t1, v, t3) == f && invOrient3D(t0, v, t2, t3) == f && invOrient3D(v, t1, t2, t3) == f);
}


//
void TetMesh::fill_memo_o3d_v_origbndt(const uint32_t v, const std::vector<uint64_t>& original_bnd_tri) {
    memo_o3d_v_origbndt[v].resize(original_bnd_tri.size());
    for (size_t i = 0; i < original_bnd_tri.size(); i++) {
        uint32_t u[3]; getFaceVertices(original_bnd_tri[i], u);
        memo_o3d_v_origbndt[v][i] = (v == u[0] || v == u[1] || v == u[2]) ? 0 : vOrient3D(u[0], u[1], u[2], v);
    }
}

//
bool TetMesh::FAST_innerSegmentCrossesInnerTriangle(const uint32_t* s_ep, const uint64_t obndt_j, const std::vector<uint64_t>& original_bnd_tri) {

    uint32_t v[3]; getFaceVertices(original_bnd_tri[obndt_j], v);
    const pointType* cs0 = vertices[s_ep[0]];
    const pointType* cs1 = vertices[s_ep[1]];
    const pointType* cv0 = vertices[v[0]];
    const pointType* cv1 = vertices[v[1]];
    const pointType* cv2 = vertices[v[2]];

    // "out of the Box" check.
    //if (three_vrts_out_of_the_box(cs0, cs1, cv0, cv1, cv2)) return false;

    if (memo_o3d_v_origbndt[s_ep[0]].empty()) fill_memo_o3d_v_origbndt(s_ep[0], original_bnd_tri);
    if (memo_o3d_v_origbndt[s_ep[1]].empty()) fill_memo_o3d_v_origbndt(s_ep[1], original_bnd_tri);
    // memo_o3d_v_origbndt[s_ep[0]][j] = invOrient3D(cv0,cv1,cv2, cs0);
    const int orient_tri_s0 = memo_o3d_v_origbndt[s_ep[0]][obndt_j];
    const int orient_tri_s1 = memo_o3d_v_origbndt[s_ep[1]][obndt_j];

    // Check if triangle vertices and at least one of the segment endpoints are coplanar:
    // in this case there is no proper intersection.
    if (orient_tri_s0 == 0 || orient_tri_s1 == 0)             return false;

    // Endpoints of one segment cannot stay both in one of the same half-space defined by the triangle.
    if (orient_tri_s0 == orient_tri_s1)                   return false;

    // Since now, endpoints are one abouve and one below the triangle-plane:
    // check if line for s0 and s1 intersect the triangle.
    return pointType::lineCrossesInnerTriangle(*cs0, *cs1, *cv0, *cv1, *cv2);
}

//
bool TetMesh::FAST_innerSegmentCrossesInnerTriangle(const pointType& cs0, const pointType& cs1, const pointType& cv0, const pointType& cv1, const pointType& cv2, int& o3d_tri_s0, int& o3d_tri_s1) const {
    // "out of the Box" check.
    //if (three_vrts_out_of_the_box(cs0, cs1, cv0, cv1, cv2)) return false;

    if (o3d_tri_s0 == 2) o3d_tri_s0 = invOrient3D(cv0, cv1, cv2, cs0);
    if (o3d_tri_s1 == 2) o3d_tri_s1 = invOrient3D(cv0, cv1, cv2, cs1);

    // Check if triangle vertices and at least one of the segment endpoints are coplanar:
    // in this case there is no proper intersection.
    if (o3d_tri_s0 == 0 || o3d_tri_s1 == 0)             return false;

    // Endpoints of one segment cannot stay both in one of the same half-space defined by the triangle.
    if (o3d_tri_s0 == o3d_tri_s1)                   return false;

    // Since now, endpoints are one abouve and one below the triangle-plane:
    // check if line for s0 and s1 intersect the triangle.
    return pointType::lineCrossesInnerTriangle(cs0, cs1, cv0, cv1, cv2);
}

//
bool TetMesh::aInnerTriASide_Crosses_InnerTriB(const pointType& vA0, const pointType& vA1, const pointType& vA2,
    const pointType& vB0, const pointType& vB1, const pointType& vB2) {
    int o3d_vA0 = 2, o3d_vA1 = 2, o3d_vA2 = 2;
    return (FAST_innerSegmentCrossesInnerTriangle(vA0, vA1, vB0, vB1, vB2, o3d_vA0, o3d_vA1) ||
        FAST_innerSegmentCrossesInnerTriangle(vA1, vA2, vB0, vB1, vB2, o3d_vA1, o3d_vA2) ||
        FAST_innerSegmentCrossesInnerTriangle(vA2, vA0, vB0, vB1, vB2, o3d_vA2, o3d_vA0));
}

//
bool TetMesh::intersectionTEST_3(const pointType& u0, const pointType& u1, const pointType& u2,
    const pointType& v0, const pointType& v1, const pointType& v2,
    const pointType& y, const int face_ori) {
    // Let <u0,u1,u2> an original cavity-boundary triangle.
    // Let <v0,v1,v2> a face of a tetrahderon X we want to build in cavity.
    // Let y a vertex of the cavity we are testing visibilty from <v0,v1,v2>.
    // If u0,u1 or u2 are inside tetrahderon of vertices v0,v1,v2,y then y is not visible.
    // Let T the tetrahderon of vertices v0,v1,v2,y.
    // Let face_ori the singe of invOrient3D(v0,v1,v2,other_vrt_of_X).
    int o3d_tf_y = invOrient3D(v0, v1, v2, y);
    if (face_ori == o3d_tf_y) {
        if (is_vrt_inside_valid_tet(u0, v0, v1, v2, y)) return true;
        if (is_vrt_inside_valid_tet(u1, v0, v1, v2, y)) return true;
        if (is_vrt_inside_valid_tet(u2, v0, v1, v2, y)) return true;
    }
    else { // if(face_ori==(-1*o3d_tf_y))
        if (is_vrt_inside_valid_tet(u0, v1, v0, v2, y)) return true;
        if (is_vrt_inside_valid_tet(u1, v1, v0, v2, y)) return true;
        if (is_vrt_inside_valid_tet(u2, v1, v0, v2, y)) return true;
    }
    // else invOrient3D(v0,v1,v2,y)==0  =>  <v0,v1,v2,y> is flat IMPOSSIBLE: EXCLUDED BEFORE
    return false;
}

// 
bool TetMesh::isTetLocallyDelaunay(const uint32_t* tet_vrts, const std::vector<uint32_t>& C_vrts, const std::vector<uint64_t>& original_bnd_tri) {
    // <v0,v1,v2> is the bnd_tri, w=v3 is the vertex we want to attach
    const pointType& cv0 = *vertices[tet_vrts[0]];
    const pointType& cv1 = *vertices[tet_vrts[1]];
    const pointType& cv2 = *vertices[tet_vrts[2]];
    const pointType& cw = *vertices[tet_vrts[3]];
    //const int o3dw = 1; // before calling this function be shure that invOrient3D(cv0,cv1,cv2,cw) = +1;

    for (const uint32_t c : C_vrts)if (!marked_vertex[c]) {

        // If c is outside tet's circumsphere there are no problems
        // If c is exactly on tet's circumsphere we use symbolic perturbation
        if (vertexInTetSphere(tet_vrts, c) < 0) continue;
        // c is INSIDE tet's circumsphere.  
        // If c is not visible from w it can stay inside tet's circumsphere.
        const pointType& cc = *vertices[c];

        // NOTE. visibilty is occluded only by cavity "original" boundary triangles.

        // TEST 1. 
        // Let pi = plane for <v0,v1,v2>.
        // If c is on the opposite half-space defined by pi and w, c it is not visible. // TRUE if <v0,v1,v2> BELONGS TO THE INITIAL BORDER. BUT WHAT IF NOT?
        if (memo_o3d[c] == 2) memo_o3d[c] = invOrient3D(cv0, cv1, cv2, cc);
        if (memo_o3d[c] < 0) continue;  // invOrient3D(cv0,cv1,cv2,cw) = +1

        // TEST 2. 
        // If any one of the four segments <v0,c>, <v1,c>, <v2,c>, <w,c> intersects
        // the interior of one original bnd_tri then c is not visible.
        // Also if one original bnd_tri edge intersects a tet_face then c is not visioble.
        bool is_visible = true;
        for (uint64_t j = 0; j < original_bnd_tri.size(); j++) {

            if (memo_o3d_v_origbndt[c].empty()) fill_memo_o3d_v_origbndt(c, original_bnd_tri);
            const int bnd_o3dc = memo_o3d_v_origbndt[c][j];
            if (bnd_o3dc == 0) continue; // tri and c are coplanar: c visibility is not occluded.

            // HALF-SPACES CHECKS
            if (memo_o3d_v_origbndt[tet_vrts[0]].empty()) fill_memo_o3d_v_origbndt(tet_vrts[0], original_bnd_tri);
            if (memo_o3d_v_origbndt[tet_vrts[1]].empty()) fill_memo_o3d_v_origbndt(tet_vrts[1], original_bnd_tri);
            if (memo_o3d_v_origbndt[tet_vrts[2]].empty()) fill_memo_o3d_v_origbndt(tet_vrts[2], original_bnd_tri);
            if (memo_o3d_v_origbndt[tet_vrts[3]].empty()) fill_memo_o3d_v_origbndt(tet_vrts[3], original_bnd_tri);

            // If all tet_vrts and c belong to the same half-space defined by orig_tri_j visibility is not occluded.
            if (memo_o3d_v_origbndt[tet_vrts[0]][j] == memo_o3d_v_origbndt[tet_vrts[1]][j] &&
                memo_o3d_v_origbndt[tet_vrts[0]][j] == memo_o3d_v_origbndt[tet_vrts[2]][j] &&
                memo_o3d_v_origbndt[tet_vrts[0]][j] == memo_o3d_v_origbndt[tet_vrts[3]][j] &&
                memo_o3d_v_origbndt[tet_vrts[0]][j] == memo_o3d_v_origbndt[c][j]) continue;


            // INTERSECTION CHECKS

            const uint32_t c_v0[2] = { c, tet_vrts[0] };
            const uint32_t c_v1[2] = { c, tet_vrts[1] };
            const uint32_t c_v2[2] = { c, tet_vrts[2] };
            const uint32_t c_w[2] = { c, tet_vrts[3] };
            if (FAST_innerSegmentCrossesInnerTriangle(c_v0, j, original_bnd_tri)) { is_visible = false; break; }
            if (FAST_innerSegmentCrossesInnerTriangle(c_v1, j, original_bnd_tri)) { is_visible = false; break; }
            if (FAST_innerSegmentCrossesInnerTriangle(c_v2, j, original_bnd_tri)) { is_visible = false; break; }
            if (FAST_innerSegmentCrossesInnerTriangle(c_w, j, original_bnd_tri)) { is_visible = false; break; }

            const uint64_t tri = original_bnd_tri[j];
            uint32_t u[3]; getFaceVertices(tri, u);
            const pointType& cu0 = *vertices[u[0]];
            const pointType& cu1 = *vertices[u[1]];
            const pointType& cu2 = *vertices[u[2]];

            //int o3d_u0 = 2, o3d_u1 = 2, o3d_u2 = 2;
            if (aInnerTriASide_Crosses_InnerTriB(cu0, cu1, cu2, cv0, cv1, cc)) { is_visible = false; break; }
            if (aInnerTriASide_Crosses_InnerTriB(cu0, cu1, cu2, cv1, cv2, cc)) { is_visible = false; break; }
            if (aInnerTriASide_Crosses_InnerTriB(cu0, cu1, cu2, cv2, cv0, cc)) { is_visible = false; break; }
            if (aInnerTriASide_Crosses_InnerTriB(cu0, cu1, cu2, cv0, cw, cc)) { is_visible = false; break; }
            if (aInnerTriASide_Crosses_InnerTriB(cu0, cu1, cu2, cv1, cw, cc)) { is_visible = false; break; }
            if (aInnerTriASide_Crosses_InnerTriB(cu0, cu1, cu2, cv2, cw, cc)) { is_visible = false; break; }

            // invOrient3D(cv0,cv1,cw,cv2) = -invOrient3D(cv0,cv1,cv2,cw) = -1 [1-perm]
            if (intersectionTEST_3(cu0, cu1, cu2, cv0, cv1, cw, cc, -1)) { is_visible = false; break; }
            // invOrient3D(cv1,cv2,cw,cv0) = -invOrient3D(cv0,cv1,cv2,cw) = -1 [3-perm]
            if (intersectionTEST_3(cu0, cu1, cu2, cv1, cv2, cw, cc, -1)) { is_visible = false; break; }
            // invOrient3D(cv2,cv0,cw,cv1) = -invOrient3D(cv0,cv1,cv2,cw) = -1 [3-perm]
            if (intersectionTEST_3(cu0, cu1, cu2, cv2, cv0, cw, cc, -1)) { is_visible = false; break; }
        }
        if (!is_visible) continue;

        return false; // c is visible and in tet's circumsphere
    }

    return true;
}

//// Returns TRUE if the tet obtained by replacing t[pos] with nv is flipped
//bool modTetFlips(const uint32_t* t, int pos, const uint32_t nv, const std::vector<pointType *>& vertices) {
//    const pointType* ct[] = { vertices[t[0]], vertices[t[1]], vertices[t[2]], vertices[t[3]] };
//    ct[pos] = vertices[nv];
//    return invOrient3D(*ct[0], *ct[1], *ct[2], *ct[3]) <= 0;
//}
//
//bool TetMesh::isTetIntersecting(const uint32_t* t, const std::vector<uint64_t>& C_bnd_tri) {
//    const pointType* ct[] = {vertices[t[0]], vertices[t[1]], vertices[t[2]], vertices[t[3]]};
//
//    for (const uint64_t tri : C_bnd_tri) {
//        uint32_t v[3]; getFaceVertices(tri, v);
//        const pointType* vc[] = { vertices[v[0]], vertices[v[1]], vertices[v[2]] };
//        const bool cvs[3] = { v[0] == t[0] || v[0] == t[1] || v[0] == t[2] || v[0] == t[3] ,
//            v[1] == t[0] || v[1] == t[1] || v[1] == t[2] || v[1] == t[3],
//            v[2] == t[0] || v[2] == t[1] || v[2] == t[2] || v[2] == t[3]
//        };
//
//        if (cvs[0] && cvs[1] && cvs[2]) continue;
//
//        for (int j1 = 0; j1 < 3; j1++) {
//            const int j2 = (j1 + 1) % 3, j3 = (j1 + 2) % 3;
//            if (cvs[j1] && cvs[j2] && !cvs[j3]) { // If two vertices are shared
//                int nf = 0;
//                int i = 0; while (t[i] == v[j1] || t[i] == v[j2]) i++;
//                if (modTetFlips(t, i, v[j3], vertices)) nf++;
//                i++; while (t[i] == v[j1] || t[i] == v[j2]) i++;
//                if (modTetFlips(t, i, v[j3], vertices)) nf++;
//                if (nf == 0) return true; // tri cuts through t
//                break; // Only one pair can be shared
//            }
//            else if (cvs[j1] && !cvs[j2] && !cvs[j3]) { // If only one vertex is shared
//                if (!modTetFlips(t, 0, v[j2], vertices) &&
//                    !modTetFlips(t, 1, v[j2], vertices) &&
//                    !modTetFlips(t, 2, v[j2], vertices) &&
//                    !modTetFlips(t, 3, v[j2], vertices)) return true; // v[j2] in volume of t
//                if (!modTetFlips(t, 0, v[j3], vertices) &&
//                    !modTetFlips(t, 1, v[j3], vertices) &&
//                    !modTetFlips(t, 2, v[j3], vertices) &&
//                    !modTetFlips(t, 3, v[j3], vertices)) return true; // v[j3] in volume of t
//
//                int cv = 0; while (t[cv] != v[j1]) cv++;
//                int cv1 = (cv + 1) & 3;
//                int cv2 = (cv + 2) & 3;
//                int cv3 = (cv + 3) & 3;
//                if (pointType::innerSegmentCrossesTriangle(*vc[j1], *vc[j2], *ct[cv1], *ct[cv2], *ct[cv3])) return true;
//                if (pointType::innerSegmentCrossesTriangle(*vc[j1], *vc[j3], *ct[cv1], *ct[cv2], *ct[cv3])) return true;
//                if (pointType::innerSegmentCrossesTriangle(*vc[j2], *vc[j3], *ct[cv1], *ct[cv2], *ct[cv3])) return true;
//                if (pointType::innerSegmentCrossesTriangle(*vc[j2], *vc[j3], *ct[cv], *ct[cv2], *ct[cv3])) return true;
//                if (pointType::innerSegmentCrossesTriangle(*vc[j2], *vc[j3], *ct[cv1], *ct[cv], *ct[cv3])) return true;
//                if (pointType::innerSegmentCrossesTriangle(*vc[j2], *vc[j3], *ct[cv1], *ct[cv2], *ct[cv])) return true;
//                if (pointType::innerSegmentCrossesTriangle(*ct[cv1], *ct[cv2], *vc[0], *vc[1], *vc[2])) return true;
//                if (pointType::innerSegmentCrossesTriangle(*ct[cv2], *ct[cv3], *vc[0], *vc[1], *vc[2])) return true;
//                if (pointType::innerSegmentCrossesTriangle(*ct[cv3], *ct[cv1], *vc[0], *vc[1], *vc[2])) return true;
//                break; // Only one pair can be unshared
//            }
//            else if (!cvs[j1] && !cvs[j2] && !cvs[j3]) { // No vertex is shared
//                for (int i=0; i<3; i++)
//                    if (!modTetFlips(t, 0, v[i], vertices) &&
//                        !modTetFlips(t, 1, v[i], vertices) &&
//                        !modTetFlips(t, 2, v[i], vertices) &&
//                        !modTetFlips(t, 3, v[i], vertices)) return true; // v[i] in volume of t
//
//                for (int k=0; k<4; k++) for (int i = 0; i < 3; i++) // If an edge of tri intersects a face of t -> true
//                    if (pointType::innerSegmentCrossesTriangle(*vc[i], *vc[(i+1)%3], *ct[(k + 1) & 3], *ct[(k + 2) & 3], *ct[(k + 3) & 3])) return true;
//
//                for (int k = 0; k < 4; k++) for (int l = k+1; l < 4; l++) // If an edge of t intersects tri
//                    if (pointType::innerSegmentCrossesTriangle(*ct[k], *ct[l], *vc[0], *vc[1], *vc[2])) return true;
//
//                break;
//            }
//        }
//    }
//
//    return false;
//}

// NOTE: in the case of the upper cavity there are also comm_tris to consider,
//        by the way it is not possible that a tetrahedron in cavity crosses such a triangle.
bool TetMesh::isTetIntersecting(const uint32_t* tet_vrts, const std::vector<uint64_t>& C_bnd_tri) {
    // For each face f in C_bnd_tetface
    // if (f is {t1,t2,t3} or {t1,t2,t4} or ...) there is no intersection, check another bnd_tri.
    // if (triangle_intersects_inner_tet(f, t)) return true; else check another bnd_tri.
    //const uint32_t tet_fv0[] = { tet_vrts[1],tet_vrts[2],tet_vrts[3] };
    //const uint32_t tet_fv1[] = { tet_vrts[0],tet_vrts[2],tet_vrts[3] };
    //const uint32_t tet_fv2[] = { tet_vrts[0],tet_vrts[1],tet_vrts[3] };
    //const uint32_t tet_fv3[] = { tet_vrts[0],tet_vrts[1],tet_vrts[2] };
    const pointType& ct0 = *vertices[tet_vrts[0]];
    const pointType& ct1 = *vertices[tet_vrts[1]];
    const pointType& ct2 = *vertices[tet_vrts[2]];
    const pointType& ct3 = *vertices[tet_vrts[3]];


    // NOTE. remember that:
    // 1) tet_fv3 = <t0,t1,t2> is the bnd_tri we are using to create a new 
    //    half-cavity tetrahedron with a suitable t3; tet_fv3 has been removed
    //    from C_bnd_tri.
    // 2) memo_o3d contains the already computed invOrient3D of cavity vertices
    //    wrt tet_fv3, i.e. tet_vrts[0],tet_vrts[1],tet_vrts[2]
    for (const uint64_t tri : C_bnd_tri) {
        uint32_t v[3]; getFaceVertices(tri, v);
        bool v0_is_t[] = { v[0] == tet_vrts[0], v[0] == tet_vrts[1], v[0] == tet_vrts[2], v[0] == tet_vrts[3] };
        bool v1_is_t[] = { v[1] == tet_vrts[0], v[1] == tet_vrts[1], v[1] == tet_vrts[2], v[1] == tet_vrts[3] };
        bool v2_is_t[] = { v[2] == tet_vrts[0], v[2] == tet_vrts[1], v[2] == tet_vrts[2], v[2] == tet_vrts[3] };

        // Check if tri a face of tet: in that case no intersection between tet and tri can occour.
        if ((v0_is_t[0] || v0_is_t[1] || v0_is_t[2] || v0_is_t[3]) &&
            (v1_is_t[0] || v1_is_t[1] || v1_is_t[2] || v1_is_t[3]) &&
            (v2_is_t[0] || v2_is_t[1] || v2_is_t[2] || v2_is_t[3])) continue;

        const pointType& cv0 = *vertices[v[0]];
        const pointType& cv1 = *vertices[v[1]];
        const pointType& cv2 = *vertices[v[2]];

        // HALF-SPACES CHECKS

        // If v0,v1 and v2 are all in the half-space defined by t0,t1,t2 and opposite to t3 there is no intersection.
        // Since both <v0,v1,v2> and <t0,t1,t2> are cavity boundary tringles none of vi can stay
        // on <t0,t1,t2> unless they coincide with t0,t1 or t2.
        if (memo_o3d[v[0]] == 2) memo_o3d[v[0]] = invOrient3D(ct0, ct1, ct2, cv0);
        if (memo_o3d[v[1]] == 2) memo_o3d[v[1]] = invOrient3D(ct0, ct1, ct2, cv1);
        if (memo_o3d[v[2]] == 2) memo_o3d[v[2]] = invOrient3D(ct0, ct1, ct2, cv2);
        // REMEMBER: memo_o3d[ tet_vrts[3] ] = invOrient3D(ct0,ct1,ct2,ct3) = +1;
        if (memo_o3d[v[0]] <= 0 && memo_o3d[v[1]] <= 0 && memo_o3d[v[2]] <= 0) continue;


        // invOrient3D(ct0,ct1,ct2,ct3) = invOrient3D(ct2,ct1,ct3,ct0) = +1
        int o3d_tf0_v0, o3d_tf0_v1, o3d_tf0_v2;

        if (v0_is_t[2] || v0_is_t[1] || v0_is_t[3]) o3d_tf0_v0 = 0;
        else if (v0_is_t[0]) o3d_tf0_v0 = 1;
        else o3d_tf0_v0 = invOrient3D(ct2, ct1, ct3, cv0);

        if (v1_is_t[2] || v1_is_t[1] || v1_is_t[3]) o3d_tf0_v1 = 0;
        else if (v1_is_t[0]) o3d_tf0_v1 = 1;
        else o3d_tf0_v1 = invOrient3D(ct2, ct1, ct3, cv1);

        if (v2_is_t[2] || v2_is_t[1] || v2_is_t[3]) o3d_tf0_v2 = 0;
        else if (v2_is_t[0]) o3d_tf0_v2 = 1;
        else o3d_tf0_v2 = invOrient3D(ct2, ct1, ct3, cv2);

        if (o3d_tf0_v0 < 0 && o3d_tf0_v1 < 0 && o3d_tf0_v2 < 0) continue;


        // invOrient3D(ct0,ct1,ct2,ct3) = invOrient3D(ct0,ct2,ct3,ct1) = +1
        int o3d_tf1_v0, o3d_tf1_v1, o3d_tf1_v2;

        if (v0_is_t[0] || v0_is_t[2] || v0_is_t[3]) o3d_tf1_v0 = 0;
        else if (v0_is_t[1]) o3d_tf1_v0 = 1;
        else o3d_tf1_v0 = invOrient3D(ct0, ct2, ct3, cv0);

        if (v1_is_t[0] || v1_is_t[2] || v1_is_t[3]) o3d_tf1_v1 = 0;
        else if (v1_is_t[1]) o3d_tf1_v1 = 1;
        else o3d_tf1_v1 = invOrient3D(ct0, ct2, ct3, cv1);

        if (v2_is_t[0] || v2_is_t[2] || v2_is_t[3]) o3d_tf1_v2 = 0;
        else if (v2_is_t[1]) o3d_tf1_v2 = 1;
        else o3d_tf1_v2 = invOrient3D(ct0, ct2, ct3, cv2);

        if (o3d_tf1_v0 < 0 && o3d_tf1_v1 < 0 && o3d_tf1_v2 < 0) continue;


        // invOrient3D(ct0,ct1,ct2,ct3) = invOrient3D(ct1,ct0,ct3,ct2) = +1
        int o3d_tf2_v0, o3d_tf2_v1, o3d_tf2_v2;

        if (v0_is_t[1] || v0_is_t[0] || v0_is_t[3]) o3d_tf2_v0 = 0;
        else if (v0_is_t[2]) o3d_tf2_v0 = 1;
        else o3d_tf2_v0 = invOrient3D(ct1, ct0, ct3, cv0);

        if (v1_is_t[1] || v1_is_t[0] || v1_is_t[3]) o3d_tf2_v1 = 0;
        else if (v1_is_t[2]) o3d_tf2_v1 = 1;
        else o3d_tf2_v1 = invOrient3D(ct1, ct0, ct3, cv1);

        if (v2_is_t[1] || v2_is_t[0] || v2_is_t[3]) o3d_tf2_v2 = 0;
        else if (v2_is_t[2]) o3d_tf2_v2 = 1;
        else o3d_tf2_v2 = invOrient3D(ct1, ct0, ct3, cv2);

        if (o3d_tf2_v0 < 0 && o3d_tf2_v1 < 0 && o3d_tf2_v2 < 0) continue;


        int o3d_v_t0, o3d_v_t1, o3d_v_t2, o3d_v_t3;
        if (v0_is_t[0] || v1_is_t[0] || v2_is_t[0]) o3d_v_t0 = 0;
        else o3d_v_t0 = invOrient3D(cv0, cv1, cv2, ct0);
        if (v0_is_t[1] || v1_is_t[1] || v2_is_t[1]) o3d_v_t1 = 0;
        else o3d_v_t1 = invOrient3D(cv0, cv1, cv2, ct1);
        if (v0_is_t[2] || v1_is_t[2] || v2_is_t[2]) o3d_v_t2 = 0;
        else o3d_v_t2 = invOrient3D(cv0, cv1, cv2, ct2);
        if (v0_is_t[3] || v1_is_t[3] || v2_is_t[3]) o3d_v_t3 = 0;
        else o3d_v_t3 = invOrient3D(cv0, cv1, cv2, ct3);

        if ((o3d_v_t0 <= 0 && o3d_v_t1 <= 0 && o3d_v_t2 <= 0 && o3d_v_t3 <= 0) ||
            (o3d_v_t0 >= 0 && o3d_v_t1 >= 0 && o3d_v_t2 >= 0 && o3d_v_t3 >= 0)) continue;


        // INTERSECTIONS CHECKS

        // NOTE. Check if tri intersects tet:
        // tet_fv3 cannot intesect tri i.e. we have to check intersections only with:
        // <t0,t3> and tri,
        // <t1,t3> and tri,
        // <t2,t3> and tri,
        // interior(tet_fv0) and tri,
        // interior(tet_fv1) and tri,
        // interior(tet_fv2) and tri,

        // Check if <t0,t3> intersects a side of tri <v0,v1,v2> 
        if (chkInnerSegmentsCross(cv0, cv1, ct0, ct3)) return true;
        if (chkInnerSegmentsCross(cv1, cv2, ct0, ct3)) return true;
        if (chkInnerSegmentsCross(cv2, cv0, ct0, ct3)) return true;
        // Check if <t1,t3> intersects a side of tri <v0,v1,v2> 
        if (chkInnerSegmentsCross(cv0, cv1, ct1, ct3)) return true;
        if (chkInnerSegmentsCross(cv1, cv2, ct1, ct3)) return true;
        if (chkInnerSegmentsCross(cv2, cv0, ct1, ct3)) return true;
        // Check if <t2,t3> intersects a side of tri <v0,v1,v2> 
        if (chkInnerSegmentsCross(cv0, cv1, ct2, ct3)) return true;
        if (chkInnerSegmentsCross(cv1, cv2, ct2, ct3)) return true;
        if (chkInnerSegmentsCross(cv2, cv0, ct2, ct3)) return true;

        // Check if a tet edge (not of the face <t0,t1,t2> inner intersects tri
        if (pointType::innerSegmentCrossesInnerTriangle(ct0, ct3, cv0, cv1, cv2)) return true;
        if (pointType::innerSegmentCrossesInnerTriangle(ct1, ct3, cv0, cv1, cv2)) return true;
        if (pointType::innerSegmentCrossesInnerTriangle(ct2, ct3, cv0, cv1, cv2)) return true;

        // A tri edge inner intersects a tet face (but not <t0,t1,t2>)
        if (pointType::innerSegmentCrossesInnerTriangle(cv0, cv1, ct0, ct1, ct3) ||
            pointType::innerSegmentCrossesInnerTriangle(cv0, cv1, ct1, ct2, ct3) ||
            pointType::innerSegmentCrossesInnerTriangle(cv0, cv1, ct2, ct0, ct3)) return true;

        if (pointType::innerSegmentCrossesInnerTriangle(cv1, cv2, ct0, ct1, ct3) ||
            pointType::innerSegmentCrossesInnerTriangle(cv1, cv2, ct1, ct2, ct3) ||
            pointType::innerSegmentCrossesInnerTriangle(cv1, cv2, ct2, ct0, ct3)) return true;

        if (pointType::innerSegmentCrossesInnerTriangle(cv2, cv0, ct0, ct1, ct3) ||
            pointType::innerSegmentCrossesInnerTriangle(cv2, cv0, ct1, ct2, ct3) ||
            pointType::innerSegmentCrossesInnerTriangle(cv2, cv0, ct2, ct0, ct3)) return true;

        // A tri vertices are on a tet faces (on AND inside - boundary have been checked above) 
        // [RARELY HAPPENS]
        if (o3d_tf0_v0 == 0 && pointType::pointInInnerTriangle(cv0, ct2, ct1, ct3)) return true;
        if (o3d_tf1_v0 == 0 && pointType::pointInInnerTriangle(cv0, ct0, ct2, ct3)) return true;
        if (o3d_tf2_v0 == 0 && pointType::pointInInnerTriangle(cv0, ct1, ct0, ct3)) return true;
        if (o3d_tf0_v1 == 0 && pointType::pointInInnerTriangle(cv1, ct2, ct1, ct3)) return true;
        if (o3d_tf1_v1 == 0 && pointType::pointInInnerTriangle(cv1, ct0, ct2, ct3)) return true;
        if (o3d_tf2_v1 == 0 && pointType::pointInInnerTriangle(cv1, ct1, ct0, ct3)) return true;
        if (o3d_tf0_v2 == 0 && pointType::pointInInnerTriangle(cv2, ct2, ct1, ct3)) return true;
        if (o3d_tf1_v2 == 0 && pointType::pointInInnerTriangle(cv2, ct0, ct2, ct3)) return true;
        if (o3d_tf2_v2 == 0 && pointType::pointInInnerTriangle(cv2, ct1, ct0, ct3)) return true;

        // NOTE. it is impossible that t0,t1,t2 or t3 belong to any 
        // of <v0,v1>, <v1,v2>, <v2,v0> sice all those segments are 
        // edges of some delmesh tetrahedron.
        // Conversly, since <t0,t1,t2,t3> is only a "potential" tetrahedron
        // it may happens that v0,v1 or v2 belong to one between <t0,t3>, <t1,t3>, <t2,t3>
        // and this would mean that <t0,t1,t2,t3> cannot be build. [RARELY HAPPENS]
        //if (!v0_is_t[0] && !v0_is_t[1] && !v0_is_t[2] && !v0_is_t[3] &&
        //    (pointType::pointInInnerSegment(cv0, ct0, ct3) || pointType::pointInInnerSegment(cv0, ct1, ct3) || pointType::pointInInnerSegment(cv0, ct2, ct3))) return true;
        //if (!v1_is_t[0] && !v1_is_t[1] && !v1_is_t[2] && !v1_is_t[3] &&
        //    (pointType::pointInInnerSegment(cv1, ct0, ct3) || pointType::pointInInnerSegment(cv1, ct1, ct3) || pointType::pointInInnerSegment(cv1, ct2, ct3))) return true;
        //if (!v2_is_t[0] && !v2_is_t[1] && !v2_is_t[2] && !v2_is_t[3] &&
        //    (pointType::pointInInnerSegment(cv2, ct0, ct3) || pointType::pointInInnerSegment(cv2, ct1, ct3) || pointType::pointInInnerSegment(cv2, ct2, ct3))) return true;
    }

    return false;
}

// Extract and orient correctly bnd_tri vertices, in order to build a new tetrahedron:
// bnd_tri vertices will be placed in the first three places of tet_node so they
// have to define a face "overlooking" the interior of the cavity.
void TetMesh::orient_bnd_tri(const uint64_t bnd_tri, uint32_t* v) const {

    const uint32_t out_k = bnd_tri & 3; // bnd_tri "face index" wrt out_tet
    const uint64_t out_tet4 = (bnd_tri >> 2) * 4; // 4 times index of tet facing bnd_tri outside the cavity

    // let ti = delmesh.tet_node[out_tet4+i]; 
    // out_tet is T = <t0,t1,t2,t3> which means that
    // for face out_k=0 -> or3D( <t1,t2,t3> , <t0> ) = -1 since <t1,t2,t3,t0> is an odd permutation of T
    // for face out_k=0 -> or3D( <t2,t1,t3> , <t0> ) = +1 since <t2,t1,t3,t0> is an even permutation of T
    // for face out_k=1 -> or3D( <t2,t0,t3> , <t1> ) = -1 since <t2,t0,t3,t1> is an odd permutation of T
    // for face out_k=1 -> or3D( <t0,t2,t3> , <t1> ) = +1 since <t0,t2,t3,t1> is an even permutation of T
    // for face out_k=2 -> or3D( <t0,t1,t3> , <t2> ) = -1 since <t0,t1,t3,t2> is an odd permutation of T
    // for face out_k=2 -> or3D( <t1,t0,t3> , <t2> ) = +1 since <t1,t0,t3,t2> is an even permutation of T
    // for face out_k=3 -> or3D( <t1,t0,t2> , <t3> ) = -1 since <t1,t0,t2,t3> is an odd permutation of T
    // for face out_k=3 -> or3D( <t0,t1,t2> , <t3> ) = +1 since <t0,t1,t2,t3> is an even permutation of T
    // NOTE. that faces with +1 or3D are overlook on the interior of out_tet so outside the cavity
    //            faces with -1 or3D are overlook on the exterior of out_tet so inside the cavity
    // NOTE. we want to make a base for a new tetrahedron inside the cavity,
    //       so we must take an odd permutation. 

    if (out_k == 0) {  // <t1,t2,t3>|<t0> face opposite to t0 overlooking cavity
        v[0] = tet_node[out_tet4 + 1];
        v[1] = tet_node[out_tet4 + 2];
        v[2] = tet_node[out_tet4 + 3];
    }
    else if (out_k == 1) { // <t2,t0,t3>|<t1> face opposite to t1 overlooking cavity
        v[0] = tet_node[out_tet4 + 2];
        v[1] = tet_node[out_tet4];
        v[2] = tet_node[out_tet4 + 3];
    }
    else if (out_k == 2) {  // <t0,t1,t3>|<t2> face opposite to t2 overlooking cavity
        v[0] = tet_node[out_tet4];
        v[1] = tet_node[out_tet4 + 1];
        v[2] = tet_node[out_tet4 + 3];
    }
    else /* k=3 */ { // <t1,t0,t2>|<t3> face opposite to t3 overlooking cavity
        v[0] = tet_node[out_tet4 + 1];
        v[1] = tet_node[out_tet4];
        v[2] = tet_node[out_tet4 + 2];
    }
}

// Returns true if the tetrahedron <bnd_tri[1,2,3],w> is locally Delaunay and do
// not intersects any triangles of the cavity boundary.
bool TetMesh::is_the_connecting_vrt(const uint32_t* bnd_tri_v, const uint32_t w,
    const std::vector<uint64_t>& C_bnd_tetfaces,
    const std::vector<uint32_t>& C_vrts,
    const std::vector<uint64_t>& original_C_bnd) {
    marked_vertex[w] = 1;
    const uint32_t tet_vrts[] = { bnd_tri_v[0], bnd_tri_v[1], bnd_tri_v[2], w };
    const pointType& cv0 = *vertices[bnd_tri_v[0]];
    const pointType& cv1 = *vertices[bnd_tri_v[1]];
    const pointType& cv2 = *vertices[bnd_tri_v[2]];
    const pointType& cw = *vertices[w];

    if (memo_o3d[w] == 2)
        memo_o3d[w] = invOrient3D(cv0, cv1, cv2, cw);

    // If invOrient3D(v0,v1,v2,w) <= 0 then w is behind bnd_tri, i.e. not visible from local cavity interior. 
    if (memo_o3d[w] <= 0) { marked_vertex[w] = 0; return false; }
    // If <v0,v1,v2,w> intersects other boundary faces cannot be created.
    if (isTetIntersecting(tet_vrts, C_bnd_tetfaces)) { marked_vertex[w] = 0; return false; }

    // If <v0,v1,v2,w> is not locally Delaunay cannot be created.
    if (!isTetLocallyDelaunay(tet_vrts, C_vrts, original_C_bnd)) { marked_vertex[w] = 0; return false; }

    marked_vertex[w] = 0;
    return true;
}

// Search a cavity vertex w to connect to the cavity boundary triangle bnd_tri
// in order to create a vew tetrahedron filling the cavity.
// This tetrahedron must be found and delmesh and cavity boundary 
// are updated accordling.
void TetMesh::connect_bnd_tri(const uint64_t bnd_tri, std::vector<uint64_t>& C_bnd_tetfaces, std::vector<uint32_t>& C_vrts,
    const std::vector<uint64_t>& original_C_bnd) {

    // NOTE. untill this function ends the vector memo_o3d will be used to save 
    // the invOrient3D signe of cavity vertices wrt bnd_tri:
    // if v is the index of a delmesh vertex belonging to cavity, memo_o3d[v] <- invOrient3D(bnd_tri,v)
    // Initially ecah element of the vector is equal to 2, then each time an orient is
    // computed is saved in the right place.

    uint32_t v[3];
    orient_bnd_tri(bnd_tri, v);

    // Mark bnd_tri's vertices
    marked_vertex[v[0]] = 1;
    marked_vertex[v[1]] = 1;
    marked_vertex[v[2]] = 1;

    // Search for one cavity vertex v3 that can be added to bnd_tri to create a new tetrahedron

    // Optimization: try first vertices that belong to cavity boundary triangles which 
    // share with bnd_tri an edge.
    uint32_t pos = 0;
    for (const uint64_t tri : C_bnd_tetfaces) {
        uint32_t u[3]; getFaceVertices(tri, u);

        for (uint32_t j = 0; j < 3; j++) {
            if (marked_vertex[u[j % 3]] && marked_vertex[u[(j + 1) % 3]]) {
                // <u[j%3],u[(j+1)%3]> is an edge of bnd_tri:
                // swap C_vrts[pos] with C_vrts[i]=u[(j+2)%3];
                uint32_t i = 0; for (; i < C_vrts.size(); i++) if (C_vrts[i] == u[(j + 2) % 3]) break;
                C_vrts[i] = C_vrts[pos];
                C_vrts[pos++] = u[(j + 2) % 3];
            }
        }

        if (pos == 3) break; // bnd_tri has three "adjacent" boundary triangles.
    }

    uint32_t v3 = UINT32_MAX;
    uint32_t C_vrts_i = 0;
    for (; C_vrts_i < C_vrts.size(); C_vrts_i++)if (!marked_vertex[C_vrts[C_vrts_i]]) {
        const uint32_t w = C_vrts[C_vrts_i];
        if (is_the_connecting_vrt(v, w, C_bnd_tetfaces, C_vrts, original_C_bnd)) { v3 = w; break; }
    }

    if (v3 == UINT32_MAX) {
        // if no valid vertex could be found: 
        // ERROR cavity does not have a CDT but this is theoretically impossible.
        ip_error("GIFT-WRAPPING FAILS\n");
    }

    if (++C_vrts_i < C_vrts.size()) {
        const pointType& cv0 = *vertices[v[0]];
        const pointType& cv1 = *vertices[v[1]];
        const pointType& cv2 = *vertices[v[2]];
        // Search for other cavity vertex that belong to new-tetrahedron insphere:
        // each time such a vertex w is founded: 
        // check if can be added to bnd_tri to create a new tetrahedron:
        // if true w is the new v3.
        for (; C_vrts_i < C_vrts.size(); C_vrts_i++)if (!marked_vertex[C_vrts[C_vrts_i]]) {
            const uint32_t w = C_vrts[C_vrts_i];
            const pointType& cw = *vertices[w];
            const pointType& cv3 = *vertices[v3];
            if (invInSphere(cv0, cv1, cv2, cv3, cw) <= 0) continue;

            if (is_the_connecting_vrt(v, w, C_bnd_tetfaces, C_vrts, original_C_bnd))
            {
                v3 = w;
            }
        }
    }

    // Reset vertices markers
    marked_vertex[v[0]] = 0;
    marked_vertex[v[1]] = 0;
    marked_vertex[v[2]] = 0;

    // Tetrahedron {v[0], v[1], v[2], v3} is constrained-Delaunay AND 
    // does not intersects any cavity-boundary-triangle: 
    // it can be added to the mesh.

    const uint64_t tet_ind = numTets();
    mark_tetrahedra.push_back(0);
    tet_node.insert(tet_node.end(), { v[0],v[1],v[2],v3 });
    tet_neigh.insert(tet_neigh.end(), { UINT64_MAX, UINT64_MAX, UINT64_MAX, UINT64_MAX });

    // update delmesh vrts -> incident tetrahedra RELATION
    inc_tet[v[0]] = tet_ind;
    inc_tet[v[1]] = tet_ind;
    inc_tet[v[2]] = tet_ind;
    inc_tet[v3] = tet_ind;

    // Update neighbours relations for <v0,v1,v2>
    tet_neigh[bnd_tri] = 4 * tet_ind + 3;
    tet_neigh[4 * tet_ind + 3] = bnd_tri;

    // Ubdate cavity boundaries: 
    // face opposite to v3 has been already removed from boundary.
    // Only <v0,v1,v3>, <v1,v2,v3>, <v2,v0,v3> have to be checked.
    for (uint32_t j = 0; j < 3; j++) {
        // Consider face opposite to vj 
        bool found_on_bnd = false;
        const uint32_t tet_fv[] = { v[(j + 1) % 3], v[(j + 2) % 3], v3 };

        for (uint64_t k = 0; k < C_bnd_tetfaces.size(); k++) {
            uint32_t bnd_fv[3]; getFaceVertices(C_bnd_tetfaces[k], bnd_fv);
            if (same_triangle(bnd_fv, tet_fv)) {
                found_on_bnd = true;
                // Update neighbours relations for <tet_fv[0],tet_fv[1],v3>
                tet_neigh[C_bnd_tetfaces[k]] = 4 * tet_ind + j;
                tet_neigh[4 * tet_ind + j] = C_bnd_tetfaces[k];
                // Remove C_bnd_tetfaces[k] from cavity boundary
                C_bnd_tetfaces.erase(C_bnd_tetfaces.begin() + k);
                break;
            }
        }

        // If tet_fv is not a boundary face add it to cavity boundary.
        if (!found_on_bnd) { C_bnd_tetfaces.push_back(4 * tet_ind + j); }
    }

    for (uint32_t i = 0; i < vertices.size(); i++) memo_o3d[i] = 2;
}

#define mark_as_common(v) MARKBIT(marked_vertex[v], 32)
#define unmark_as_common(v) UNMARKBIT(marked_vertex[v], 32)
#define belong_to_common(v) (ISMARKEDBIT(marked_vertex[v], 32))

void TetMesh::giftWrapping(const std::vector<uint32_t>& comm_vrts, std::vector<uint32_t>& C1_vrts, std::vector<uint32_t>& C2_vrts,
    const std::vector<uint64_t>& C_bnd_tetface, const uint64_t n_cav_tets, const uint64_t n_C1_bnd_tetface) {

    //std::cout << "GIFT-WRAPPING\n";

    if (memo_o3d.empty()) {
        memo_o3d.resize(vertices.size(), 2); // initialize a vector to avoid recomputing orient3d long as delmesh vertices
        memo_o3d_v_origbndt.resize(vertices.size());
    }

    // DELETE CAVITY TETS (we assume them to be at the end of the vector)
    const uint64_t n_valid_delmesh_tets = numTets() - n_cav_tets;
    tet_node.resize(n_valid_delmesh_tets * 4);
    tet_neigh.resize(n_valid_delmesh_tets * 4);
    mark_tetrahedra.resize(n_valid_delmesh_tets);

    // FILL THE "UPPER" HOLE AND COLLECT TRIANGLES ON THE PLCFACE
    std::vector<uint64_t> C1_bnd_tri(C_bnd_tetface.begin(), C_bnd_tetface.begin() + n_C1_bnd_tetface);
    std::vector<uint64_t> original_C1_bnd_tri(C1_bnd_tri);
    std::vector<uint64_t> comm_tris;
    C1_vrts.insert(C1_vrts.begin(), comm_vrts.begin(), comm_vrts.end());

    while (!C1_bnd_tri.empty()) {
        const uint64_t bnd_tri = C1_bnd_tri.back();
        C1_bnd_tri.pop_back();

        // check wherever bnd_tri is on the plc_face or not
        uint32_t v[3]; getFaceVertices(bnd_tri, v);
        bool is_comm_tri = false;
        for (uint32_t i = 0; i < comm_vrts.size(); i++) mark_as_common(comm_vrts[i]);
        if (belong_to_common(v[0]) && belong_to_common(v[1]) && belong_to_common(v[2])) {
            comm_tris.push_back(bnd_tri);
            is_comm_tri = true;
        }
        for (uint32_t i = 0; i < comm_vrts.size(); i++) unmark_as_common(comm_vrts[i]);

        if (!is_comm_tri) connect_bnd_tri(bnd_tri, C1_bnd_tri, C1_vrts, original_C1_bnd_tri);
    }

    for (const uint32_t v : C1_vrts) memo_o3d_v_origbndt[v].clear(); // Reset memo_o3d_v_origbndt

    // 5 - FILL THE "LOWER" HOLE
    std::vector<uint64_t> C2_bnd_tri(C_bnd_tetface.begin() + n_C1_bnd_tetface, C_bnd_tetface.end());
    std::vector<uint64_t> original_C2_bnd_tri(C2_bnd_tri);
    // Add bnd_tris on the plcface
    C2_bnd_tri.insert(C2_bnd_tri.end(), comm_tris.begin(), comm_tris.end());
    C2_vrts.insert(C2_vrts.begin(), comm_vrts.begin(), comm_vrts.end());

    while (!C2_bnd_tri.empty()) {
        //printf("\r%lu bnd_tri to connect (lower cav)      ", C2_bnd_tri.size()); fflush(stdout);
        const uint64_t bnd_tri = C2_bnd_tri.back();
        C2_bnd_tri.pop_back();

        connect_bnd_tri(bnd_tri, C2_bnd_tri, C2_vrts, original_C2_bnd_tri);
    }
    for (const uint32_t v : C2_vrts) memo_o3d_v_origbndt[v].clear(); // Reset memo_o3d_v_origbndt

    // from now on delmesh connectivity is restored
}

bool TetMesh::isUpperCavityTet(const uint64_t t, std::vector<int>& v_orient) const {
    uint32_t v[3];
    getFaceVertices(t, v);
    return v_orient[v[0]] >= 0 && v_orient[v[1]] >= 0 && v_orient[v[2]] >= 0;
}

bool TetMesh::isLowerCavityTet(const uint64_t t, std::vector<int>& v_orient) const {
    uint32_t v[3];
    getFaceVertices(t, v);
    return v_orient[v[0]] <= 0 && v_orient[v[1]] <= 0 && v_orient[v[2]] <= 0;
}

void TetMesh::recoverFaceGiftWrap(std::vector<uint64_t>& i_tets, std::vector<int>& v_orient) {
    if (marked_vertex.empty()) marked_vertex.resize(numVertices(), 0);

    // Create vector 'top_faces' and 'bottom_faces'
    std::vector<uint64_t> top_faces, bottom_faces;
    for (uint64_t t : i_tets) mark_Tet_1(t);

    // Move all tets to remove to tail
    uint32_t last = numTets() - 1;
    for (uint64_t& t : i_tets) {
        while (is_marked_Tet_1(last)) last--;
        if (t < last) {
            swapTets(t, last);
            t = last;
        }
    }
    while (is_marked_Tet_1(last)) last--;

    last++; // This is now the numtets after having removed toremove

    for (uint64_t t : i_tets) {
        const uint64_t* neigh = getTetNeighs(t << 2);
        for (int i = 0; i < 4; i++) if (!is_marked_Tet_1(neigh[i] >> 2)) {
            if (isUpperCavityTet(neigh[i], v_orient)) top_faces.push_back(neigh[i]);
            else {
                assert(isLowerCavityTet(neigh[i], v_orient)); // Must be either
                bottom_faces.push_back(neigh[i]);
            }
        }
    }

    for (uint64_t t : i_tets) unmark_Tet_1(t);

    std::vector<uint32_t> top_vertices, bottom_vertices, face_vertices;

    for (uint64_t t : i_tets) {
        for (int i = 0; i < 4; i++) {
            const uint32_t v = tet_node[t * 4 + i];
            if (!marked_vertex[v]) {
                marked_vertex[v] = 1;
                if (v_orient[v] > 0) top_vertices.push_back(v);
                else if (v_orient[v] < 0) bottom_vertices.push_back(v);
                else face_vertices.push_back(v);
            }
        }
    }

    for (uint32_t w : top_vertices) marked_vertex[w] = 0;
    for (uint32_t w : bottom_vertices) marked_vertex[w] = 0;
    for (uint32_t w : face_vertices) marked_vertex[w] = 0;

    std::vector<uint64_t> all_bnd_faces(top_faces.begin(), top_faces.end());
    all_bnd_faces.insert(all_bnd_faces.end(), bottom_faces.begin(), bottom_faces.end());

    giftWrapping(face_vertices, top_vertices, bottom_vertices, all_bnd_faces, i_tets.size(), top_faces.size());
}
