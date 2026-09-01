#pragma once

// Uncommenting the following macro definition makes the code use modified parts of hxt_SeqDel (Copyright (C) 2018 Célestin Marot).
// hxt_SeqDel is a sequential Delaunay triangulator hosted at https://git.immc.ucl.ac.be/hextreme/hxt_seqdel as of 2020.
// hxt_SeqDel is GPL licensed, meaning that if you uncomment the following line you accept the terms of the GPL license for
// the whole code which uses this library.
// If you need to use this code under the less restrictive LGPL license, please comment the following line.
// This will make the code slightly slower.
#define USE_MAROTS_METHOD

#include "delaunay.h"

typedef genericPoint pointType;
typedef implicitPoint3D_LNC implicitPoint_LNC;
typedef explicitPoint3D explicitPoint;

typedef TetMesh_t<pointType *> TetMesh_base;


// This extends basicVec3d by adding an approximated interface to convert from pointType

class vector3d : public basicVec3d {
public:
    inline vector3d() : basicVec3d() { }
    inline vector3d(const double x, const double y, const double z) : basicVec3d(x, y, z) {}
    inline vector3d(const pointType* p) : basicVec3d() { p->getApproxXYZCoordinates(coord[0], coord[1], coord[2]); }

    // TRUE if smallest sphere by p,q,r is larger than smallest sphere by p,q,s
    static inline bool hasLargerSphere(const pointType* p, const pointType* q, const pointType* r, const pointType* s) {
        return basicVec3d::hasLargerSphere(vector3d(p), vector3d(q), vector3d(r), vector3d(s));
    }

    // TRUE if p is closer to q than to r
    static bool isCloserThan(const pointType* p, const pointType* q, const pointType* r) {
        return basicVec3d::isCloserThan(vector3d(p), vector3d(q), vector3d(r));
    }

    // TRUE if distance p-q is at most half the distance p-r
    static bool isAtMostTwiceDistanceThan(const pointType* p, const pointType* q, const pointType* r) {
        return basicVec3d::isAtMostTwiceDistanceThan(vector3d(p), vector3d(q), vector3d(r));
    }
};


// This extends the basic TetMesh data structure provided by Delaunay3D

class TetMesh : public TetMesh_base {
public:
    mutable std::vector<bool> cornerMask; // Marks on corners to represent constraints

    // Gift-wrapping fields
    std::vector<int> memo_o3d;
    std::vector<std::vector<int>> memo_o3d_v_origbndt; // i-th vector is {orient3d(original_cav_tri_1,v_i), ..., orient3d(original_cav_tri_n,v_i)}

    // Constructors
    TetMesh() : TetMesh_base() {};
    TetMesh(bool h) : TetMesh_base(h) {};

    /////// Global functions ///////

    // Marks internal tets ad DT_IN and external as DT_OUT and return the number of internal tets.
    // cornerMask must be TRUE for each corner whose opposite face is a constraint.
    size_t markInnerTets(uint64_t single_start = UINT64_MAX);

    // Same as above, but marks as DT_OUT only ghosts and tets that can be reached from ghosts
    // without crossing any constraint. All other tets are marked as DT_IN.
    void markInnerTetsNonManifold();

    // Clear deleted vertices after removal
    void removeDelVertices();

    /////// Local (element-based) functions ///////

    // TRUE if tet contains p (boundary included)
    bool tetContainsPoint(uint64_t tet, const pointType* p) const;

    void getConstrainedCavity(const uint32_t v_id, uint64_t& tet, std::vector<uint64_t>& cavityCorners, uint32_t& cv1, uint32_t& cv2, uint32_t& cv3);
    // fill 'adjacencies' with consecutive pairs of edge-adjacent tets in cavityCorners

    // Collect all the vertices contained in the smallest sphere by ep0 and ep1
    // and return the one generating the largest circumcircle with ep0 and ep1.
    // Init tet with one tet having the encroaching point
    uint32_t findEncroachingPoint_inexact(const uint32_t ep0, const uint32_t ep1, uint64_t& tet) const;
    uint32_t findEncroachingPoint_exact(const uint32_t ep0, const uint32_t ep1, uint64_t& tet) const;
    uint32_t findEncroachingPoint(const uint32_t ep0, const uint32_t ep1, uint64_t& tet) const;

    // Set of functions implementing the face recovery by gift-wrapping
    void fill_memo_o3d_v_origbndt(const uint32_t v, const std::vector<uint64_t>& original_bnd_tri);
    bool FAST_innerSegmentCrossesInnerTriangle(const uint32_t* s_ep, const uint64_t obndt_j, const std::vector<uint64_t>& original_bnd_tri);
    bool FAST_innerSegmentCrossesInnerTriangle(const pointType& cs0, const pointType& cs1, const pointType& cv0, const pointType& cv1, const pointType& cv2, int& o3d_tri_s0, int& o3d_tri_s1) const;
    bool aInnerTriASide_Crosses_InnerTriB(const pointType& vA0, const pointType& vA1, const pointType& vA2, const pointType& vB0, const pointType& vB1, const pointType& vB2);
    bool intersectionTEST_3(const pointType& u0, const pointType& u1, const pointType& u2,
        const pointType& v0, const pointType& v1, const pointType& v2,
        const pointType& y, const int face_ori);
    bool isTetLocallyDelaunay(const uint32_t* tet_vrts, const std::vector<uint32_t>& C_vrts, const std::vector<uint64_t>& original_bnd_tri);
    bool isTetIntersecting(const uint32_t* tet_vrts, const std::vector<uint64_t>& C_bnd_tri);
    void orient_bnd_tri(const uint64_t bnd_tri, uint32_t* v) const;
    bool is_the_connecting_vrt(const uint32_t* bnd_tri_v, const uint32_t w, const std::vector<uint64_t>& C_bnd_tetfaces,
        const std::vector<uint32_t>& C_vrts, const std::vector<uint64_t>& original_C_bnd);
    void connect_bnd_tri(const uint64_t bnd_tri, std::vector<uint64_t>& C_bnd_tetfaces, std::vector<uint32_t>& C_vrts,
        const std::vector<uint64_t>& original_C_bnd);

    void giftWrapping(const std::vector<uint32_t>& comm_vrts, std::vector<uint32_t>& C1_vrts, std::vector<uint32_t>& C2_vrts,
        const std::vector<uint64_t>& C_bnd_tetface, const uint64_t n_cav_tets, const uint64_t n_C1_bnd_tetface);
    bool isUpperCavityTet(const uint64_t t, std::vector<int>& v_orient) const;
    bool isLowerCavityTet(const uint64_t t, std::vector<int>& v_orient) const;
    void recoverFaceGiftWrap(std::vector<uint64_t>& i_tets, std::vector<int>& v_orient);
};

inline size_t TetMesh::markInnerTets(uint64_t single_start) {
    std::vector<uint64_t> C;

    // All ghosts are DT_OUT
    for (size_t i = 0; i < numTets(); i++)
        mark_tetrahedra[i] = (isGhost(i)) ? DT_OUT : DT_UNKNOWN;

    if (single_start != UINT64_MAX) C.push_back(single_start);
    else for (size_t i = 0; i < numTets(); i++)
        if (mark_tetrahedra[i] == DT_OUT) C.push_back(i);

    for (size_t i = 0; i < C.size(); i++) {
        uint64_t t = C[i];
        for (int j = 0; j < 4; j++) {
            const uint64_t n = tet_neigh[t * 4 + j];
            const uint64_t n2 = n >> 2;
            if (mark_tetrahedra[n2] == DT_UNKNOWN) {
                if (!cornerMask[n]) {
                    mark_tetrahedra[n2] = mark_tetrahedra[t];
                }
                else {
                    mark_tetrahedra[n2] = ((mark_tetrahedra[t] == DT_IN) ? (DT_OUT) : (DT_IN));
                }
                C.push_back(n2);
            }
        }
    }

    return std::count(mark_tetrahedra.begin(), mark_tetrahedra.end(), DT_IN);
}

inline uint32_t TetMesh::findEncroachingPoint_inexact(const uint32_t ep0, const uint32_t ep1, uint64_t& tet_e) const {
    static std::vector<uint64_t> enc_queue; // Static to avoid reallocation upon each call

    // Start collecting tetrahedra incident at the endpoints
    VT(ep0, enc_queue);

    for (uint64_t j : enc_queue) mark_Tet_1(j);

    const vector3d p0 = vertices[ep0];
    const vector3d p1 = vertices[ep1];
    const double eslen = (p0 - p1).sq_length();

    vector3d ep;
    uint32_t enc_pt_i = UINT32_MAX;

    marked_vertex[ep0] = marked_vertex[ep1] = 1;

    // Collect all encroaching points while expanding around insphere vertices
    for (uint32_t ti = 0; ti < enc_queue.size(); ti++) {
        const uint64_t tet = enc_queue[ti];
        const uint64_t tb = tet << 2;

        // Check each tet vertex for 'insphereness' and keep track of the one with largest sphere
        const uint32_t* tn = tet_node.data() + tb;
        for (uint32_t i = 0; i < 4; i++) {
            const uint32_t ui = tn[i];
            if (!marked_vertex[ui]) {
                const vector3d& pui = vertices[ui];
                if (((pui - p0).sq_length() + (pui - p1).sq_length()) <= eslen) {
                    marked_vertex[ui] = 1;
                    if (enc_pt_i == UINT32_MAX || basicVec3d::hasLargerSphere(p0, p1, pui, ep)) {
                        ep = pui; enc_pt_i = ui;
                        tet_e = tb;
                    }
                }
                else marked_vertex[ui] = 2;
            }
        }

        const int nvmask[] = { (marked_vertex[tn[0]] == 1), (marked_vertex[tn[1]] == 1), (marked_vertex[tn[2]] == 1), (marked_vertex[tn[3]] == 1) };
        const int totmarkeda = nvmask[0] + nvmask[1] + nvmask[2] + nvmask[3];

        // Expand on adjacent tets if at least one common vertex is insphere
        const uint64_t* tg = tet_neigh.data() + tb;
        for (uint32_t i = 0; i < 4; i++) {
            const uint64_t nc = tg[i];
            const uint64_t n = nc >> 2;
            if (is_marked_Tet_1(n) == 2 || tet_node[nc] == INFINITE_VERTEX) continue;
            const int totmarked = totmarkeda - nvmask[i];
            if (totmarked) {
                mark_Tet_1(n);
                enc_queue.push_back(n);
            }
        }
    }

    // Clear all marks
    marked_vertex[ep0] = marked_vertex[ep1] = 0;
    for (uint64_t j : enc_queue) {
        unmark_Tet_1(j);
        j <<= 2;
        marked_vertex[tet_node[j++]] = 0;
        marked_vertex[tet_node[j++]] = 0;
        marked_vertex[tet_node[j++]] = 0;
        marked_vertex[tet_node[j]] = 0;
    }
    enc_queue.clear();

    return enc_pt_i;
}

inline uint32_t TetMesh::findEncroachingPoint_exact(const uint32_t ep0, const uint32_t ep1, uint64_t& tet_e) const {
    static std::vector<uint64_t> enc_queue; // Static to avoid reallocation upon each call

    // Start collecting tetrahedra incident at the endpoints
    VT(ep0, enc_queue);

    for (uint64_t j : enc_queue) mark_Tet_1(j);

    uint32_t enc_pt_i = UINT32_MAX;

    marked_vertex[ep0] = marked_vertex[ep1] = 1;

    // Collect all encroaching points while expanding around insphere vertices
    for (uint32_t ti = 0; ti < enc_queue.size(); ti++) {
        const uint64_t tet = enc_queue[ti];
        const uint64_t tb = tet << 2;

        // Check each tet vertex for 'insphereness' and keep track of the one with largest sphere
        const uint32_t* tn = tet_node.data() + tb;
        for (uint32_t i = 0; i < 4; i++) {
            const uint32_t ui = tn[i];
            if (!marked_vertex[ui]) {
                if (genericPoint::dotProductSign3D(*vertices[ep0], *vertices[ep1], *vertices[ui]) <= 0) {
                    marked_vertex[ui] = 1;
                    if (enc_pt_i == UINT32_MAX || pointType::inGabrielSphere(*vertices[ui], *vertices[enc_pt_i], *vertices[ep0], *vertices[ep1]) < 0) {
                        enc_pt_i = ui;
                        tet_e = tb;
                    }
                }
                else marked_vertex[ui] = 2;
            }
        }

        const int nvmask[] = { (marked_vertex[tn[0]] == 1), (marked_vertex[tn[1]] == 1), (marked_vertex[tn[2]] == 1), (marked_vertex[tn[3]] == 1) };
        const int totmarkeda = nvmask[0] + nvmask[1] + nvmask[2] + nvmask[3];

        // Expand on adjacent tets if at least one common vertex is insphere
        const uint64_t* tg = tet_neigh.data() + tb;
        for (uint32_t i = 0; i < 4; i++) {
            const uint64_t nc = tg[i];
            const uint64_t n = nc >> 2;
            if (is_marked_Tet_1(n) == 2 || tet_node[nc] == INFINITE_VERTEX) continue;
            const int totmarked = totmarkeda - nvmask[i];
            if (totmarked) {
                mark_Tet_1(n);
                enc_queue.push_back(n);
            }
        }
    }

    // Clear all marks
    marked_vertex[ep0] = marked_vertex[ep1] = 0;
    for (uint64_t j : enc_queue) {
        unmark_Tet_1(j);
        j <<= 2;
        marked_vertex[tet_node[j++]] = 0;
        marked_vertex[tet_node[j++]] = 0;
        marked_vertex[tet_node[j++]] = 0;
        marked_vertex[tet_node[j]] = 0;
    }
    enc_queue.clear();

    return enc_pt_i;
}


// Filtered function to find encroaching points.
// First, try with approximated but fast vector3d functions. If no points are found this way, revert to
// the slower but exact method using implecit points.
inline uint32_t TetMesh::findEncroachingPoint(const uint32_t ep0, const uint32_t ep1, uint64_t& tet_e) const {
    uint32_t enc_pt_i = findEncroachingPoint_inexact(ep0, ep1, tet_e);
    if (enc_pt_i == UINT32_MAX) enc_pt_i = findEncroachingPoint_exact(ep0, ep1, tet_e);
    return enc_pt_i;
}
