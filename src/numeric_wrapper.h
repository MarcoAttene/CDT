#ifndef NUMERIC_WRAPPER
#define NUMERIC_WRAPPER

#define USE_INDIRECT_PREDS
//#define USE_DOUBLE
//#define USE_LAZY_CORE
//#define USE_PLAIN_CORE
//#define USE_LAZY_GMPQ
//#define USE_PLAIN_GMPQ

#include "implicit_point.h"

#ifdef USE_INDIRECT_PREDS
typedef genericPoint pointType;
typedef implicitPoint3D_LNC implicitPoint_LNC;
typedef explicitPoint3D explicitPoint;
#endif

#ifdef USE_DOUBLE
typedef double coord_t;
#define GET_DOUBLE_VAL(a) (a)
#define GET_SIGN(a) (((a)>0) - ((a)<0))
#define GET_FABS(a) (((a)>0) ? (a) : (-a))
#define GET_SQRT(a) (sqrt(a))
#endif

#ifdef USE_LAZY_CORE
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/CORE_Expr.h>
typedef CGAL::Lazy_exact_nt< CORE::Expr> coord_t;
#define GET_DOUBLE_VAL(a) (a).approx().inf()
#define GET_SIGN(a) (((a)>0) - ((a)<0))
#define GET_FABS(a) (((a)>0) ? (a) : (-a))
#define GET_SQRT(a) (sqrt(a))
#endif

#ifdef USE_PLAIN_CORE
#include <CGAL/CORE_Expr.h>
typedef CORE::Expr coord_t;
#define GET_DOUBLE_VAL(a) (a).doubleValue()
#define GET_SIGN(a) ((a).sign())
#define GET_FABS(a) (fabs(a))
#define GET_SQRT(a) (sqrt(a))
#endif

#ifdef USE_LAZY_GMPQ
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/Gmpq.h>
typedef CGAL::Lazy_exact_nt< CGAL::Gmpq > coord_t;
#define GET_DOUBLE_VAL(a) (a).approx().inf()
#define GET_SIGN(a) (((a)>0) - ((a)<0))
#define GET_FABS(a) (((a)>0) ? (a) : (-a))
#define GET_SQRT(a) (::sqrt(GET_DOUBLE_VAL(a)))
#endif


#ifdef USE_PLAIN_GMPQ
#include <CGAL/Gmpxx.h>
typedef mpq_class coord_t;
#define GET_DOUBLE_VAL(a) (a).get_d()
#define GET_SIGN(a) (((a)>0) - ((a)<0))
#define GET_FABS(a) (((a)>0) ? (a) : (-a))
#define GET_SQRT(a) (::sqrt(GET_DOUBLE_VAL(a)))
#endif

#ifndef USE_INDIRECT_PREDS

#include <cstring>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

extern void ip_error(const char* msg);

#pragma intrinsic(fabs)

#define INFINITE_VERTEX UINT32_MAX
#define DT_UNKNOWN  0
#define DT_OUT  1
#define DT_IN  2
#define EXPECTED_VT_SIZE 128

class pointType {
public:
    coord_t coords[3];

    pointType() {}

    pointType(double a, double b, double c) : coords{ a, b, c } {}

    pointType(const pointType& a, const pointType& b, const coord_t& t) {
        for (int i = 0; i < 3; i++)
            coords[i] = (a.coords[i] * (1 - t) + (b.coords[i] * t));
    }

    const coord_t& X() const { return coords[0]; }
    const coord_t& Y() const { return coords[1]; }
    const coord_t& Z() const { return coords[2]; }

    bool isExplicit3D() const { return GET_DOUBLE_VAL(coords[0]) == coords[0] && GET_DOUBLE_VAL(coords[1]) == coords[1] && GET_DOUBLE_VAL(coords[2]) == coords[2]; }
    pointType toExplicit3D() const { return *this; }
    void apapExplicit(pointType& p) const { 
        p.coords[0] = GET_DOUBLE_VAL(coords[0]); 
        p.coords[1] = GET_DOUBLE_VAL(coords[1]);
        p.coords[2] = GET_DOUBLE_VAL(coords[2]);
    }

    coord_t& operator[](int i) { return coords[i]; }
    const coord_t& operator[](int i) const { return coords[i]; }

    bool operator==(const pointType& p) const {
        return coords[0] == p.coords[0] && coords[1] == p.coords[1] && coords[2] == p.coords[2]; 
    }

    void getApproxXYZCoordinates(double& x, double& y, double& z, bool apap=true) const {
        x = GET_DOUBLE_VAL(coords[0]);
        y = GET_DOUBLE_VAL(coords[1]);
        z = GET_DOUBLE_VAL(coords[2]);
    }

    static int dotProductSign3D(const pointType& a, const pointType& b, const pointType& c) {
        coord_t ac[3] = { a.coords[0] - c.coords[0], a.coords[1] - c.coords[1], a.coords[2] - c.coords[2] };
        coord_t bc[3] = { b.coords[0] - c.coords[0], b.coords[1] - c.coords[1], b.coords[2] - c.coords[2] };
        return GET_SIGN(ac[0] * bc[0] + ac[1] * bc[1] + ac[2] * bc[2]);
    }

    static int maxComponentInTriangleNormal(
        const coord_t& ov10, const coord_t& ov11, const coord_t& ov12,
        const coord_t& ov20, const coord_t& ov21, const coord_t& ov22,
        const coord_t& ov30, const coord_t& ov31, const coord_t& ov32)
    {
        const coord_t v3x = ov30 - ov20;
        const coord_t v3y = ov31 - ov21;
        const coord_t v3z = ov32 - ov22;
        const coord_t v2x = ov20 - ov10;
        const coord_t v2y = ov21 - ov11;
        const coord_t v2z = ov22 - ov12;
        const coord_t nvx1 = v2y * v3z;
        const coord_t nvx2 = v2z * v3y;
        const coord_t nvx = nvx1 - nvx2;
        const coord_t nvy1 = v3x * v2z;
        const coord_t nvy2 = v3z * v2x;
        const coord_t nvy = nvy1 - nvy2;
        const coord_t nvz1 = v2x * v3y;
        const coord_t nvz2 = v2y * v3x;
        const coord_t nvz = nvz1 - nvz2;
        const coord_t nvxc = GET_FABS(nvx);
        const coord_t nvyc = GET_FABS(nvy);
        const coord_t nvzc = GET_FABS(nvz);
        coord_t nv = nvxc;
        if (nvyc > nv) nv = nvyc;
        if (nvzc > nv) nv = nvzc;

        if (nv == nvxc) return 0;
        if (nv == nvyc) return 1;
        assert(nv == nvzc);
        return 2;
    }

    static int orient2D(const coord_t& p0, const coord_t& p1, const coord_t& q0, const coord_t& q1, const coord_t& r0, const coord_t& r1)
    {
#ifdef USE_DOUBLE
        return orient2d(p0, p1, q0, q1, r0, r1);
#else
        const coord_t dl = (q0 - p0) * (r1 - p1);
        const coord_t dr = (q1 - p1) * (r0 - p0);
        const coord_t det = dl - dr;

        return (det < 0) - (det > 0);
#endif
    }

    static int orient2Dxy(const pointType& p, const pointType& q, const pointType& r)
    {
        return orient2D(p[0], p[1], q[0], q[1], r[0], r[1]);
    }

    static int orient2Dyz(const pointType& p, const pointType& q, const pointType& r)
    {
        return orient2D(p[1], p[2], q[1], q[2], r[1], r[2]);
    }

    static int orient2Dzx(const pointType& p, const pointType& q, const pointType& r)
    {
        return orient2D(p[2], p[0], q[2], q[0], r[2], r[0]);
    }

    static bool misaligned(const pointType& A, const pointType& B, const pointType& C) {
        return (orient2Dxy(A, B, C) || orient2Dyz(A, B, C) || orient2Dzx(A, B, C));
    }

    static int orient3D(const pointType& p, const pointType& q, const pointType& r, const pointType& s) {
#ifdef USE_DOUBLE
        return orient3d(p.X(), p.Y(), p.Z(), q.X(), q.Y(), q.Z(), r.X(), r.Y(), r.Z(), s.X(), s.Y(), s.Z());
#else
        const coord_t fadx = q[0] - p[0], fbdx = r[0] - p[0], fcdx = s[0] - p[0];
        const coord_t fady = q[1] - p[1], fbdy = r[1] - p[1], fcdy = s[1] - p[1];
        const coord_t fadz = q[2] - p[2], fbdz = r[2] - p[2], fcdz = s[2] - p[2];

        const coord_t fbdxcdy = fbdx * fcdy * fadz; const coord_t fcdxbdy = fcdx * fbdy * fadz;
        const coord_t fcdxady = fcdx * fady * fbdz; const coord_t fadxcdy = fadx * fcdy * fbdz;
        const coord_t fadxbdy = fadx * fbdy * fcdz; const coord_t fbdxady = fbdx * fady * fcdz;

        const coord_t det = (fbdxcdy - fcdxbdy) + (fcdxady - fadxcdy) + (fadxbdy - fbdxady);

        return (det < 0) - (det > 0);
#endif
    }

    static int inSphere(const pointType& pa, const pointType& pb, const pointType& pc, const pointType& pd, const pointType& pe)
    {
#ifdef USE_DOUBLE
        return ::inSphere(pa.X(), pa.Y(), pa.Z(), pb.X(), pb.Y(), pb.Z(), pc.X(), pc.Y(), pc.Z(), pd.X(), pd.Y(), pd.Z(), pe.X(), pe.Y(), pe.Z());
#else
        const coord_t aex = pa[0] - pe[0], bex = pb[0] - pe[0], cex = pc[0] - pe[0], dex = pd[0] - pe[0];
        const coord_t aey = pa[1] - pe[1], bey = pb[1] - pe[1], cey = pc[1] - pe[1], dey = pd[1] - pe[1];
        const coord_t aez = pa[2] - pe[2], bez = pb[2] - pe[2], cez = pc[2] - pe[2], dez = pd[2] - pe[2];

        const coord_t aexbey = aex * bey;
        const coord_t bexaey = bex * aey;
        const coord_t ab = aexbey - bexaey;
        const coord_t bexcey = bex * cey;
        const coord_t cexbey = cex * bey;
        const coord_t bc = bexcey - cexbey;
        const coord_t cexdey = cex * dey;
        const coord_t dexcey = dex * cey;
        const coord_t cd = cexdey - dexcey;
        const coord_t dexaey = dex * aey;
        const coord_t aexdey = aex * dey;
        const coord_t da = dexaey - aexdey;
        const coord_t aexcey = aex * cey;
        const coord_t cexaey = cex * aey;
        const coord_t ac = aexcey - cexaey;
        const coord_t bexdey = bex * dey;
        const coord_t dexbey = dex * bey;
        const coord_t bd = bexdey - dexbey;
        const coord_t abc1 = aez * bc;
        const coord_t abc2 = bez * ac;
        const coord_t abc3 = cez * ab;
        const coord_t abc4 = abc1 + abc3;
        const coord_t abc = abc4 - abc2;
        const coord_t bcd1 = bez * cd;
        const coord_t bcd2 = cez * bd;
        const coord_t bcd3 = dez * bc;
        const coord_t bcd4 = bcd1 + bcd3;
        const coord_t bcd = bcd4 - bcd2;
        const coord_t cda1 = cez * da;
        const coord_t cda2 = dez * ac;
        const coord_t cda3 = aez * cd;
        const coord_t cda4 = cda1 + cda3;
        const coord_t cda = cda4 + cda2;
        const coord_t dab1 = dez * ab;
        const coord_t dab2 = aez * bd;
        const coord_t dab3 = bez * da;
        const coord_t dab4 = dab1 + dab3;
        const coord_t dab = dab4 + dab2;
        const coord_t al1 = aex * aex;
        const coord_t al2 = aey * aey;
        const coord_t al3 = aez * aez;
        const coord_t al4 = al1 + al2;
        const coord_t alift = al4 + al3;
        const coord_t bl1 = bex * bex;
        const coord_t bl2 = bey * bey;
        const coord_t bl3 = bez * bez;
        const coord_t bl4 = bl1 + bl2;
        const coord_t blift = bl4 + bl3;
        const coord_t cl1 = cex * cex;
        const coord_t cl2 = cey * cey;
        const coord_t cl3 = cez * cez;
        const coord_t cl4 = cl1 + cl2;
        const coord_t clift = cl4 + cl3;
        const coord_t dl1 = dex * dex;
        const coord_t dl2 = dey * dey;
        const coord_t dl3 = dez * dez;
        const coord_t dl4 = dl1 + dl2;
        const coord_t dlift = dl4 + dl3;
        const coord_t ds1 = dlift * abc;
        const coord_t ds2 = clift * dab;
        const coord_t dl = ds2 - ds1;
        const coord_t dr1 = blift * cda;
        const coord_t dr2 = alift * bcd;
        const coord_t dr = dr2 - dr1;
        const coord_t det = dl + dr;

        return ((det < 0) - (det > 0));
#endif
    }

    static bool lineCrossesTriangle(const pointType& s1, const pointType& s2, const pointType& v1, const pointType& v2, const pointType& v3)
    {
        const int o1 = orient3D(s1, s2, v1, v2);
        const int o2 = orient3D(s1, s2, v2, v3);
        if ((o1 > 0 && o2 < 0) || (o1 < 0 && o2 > 0)) return false;
        const int o3 = orient3D(s1, s2, v3, v1);
        if ((o1 > 0 && o3 < 0) || (o1 < 0 && o3 > 0)) return false;
        if ((o2 > 0 && o3 < 0) || (o2 < 0 && o3 > 0)) return false;
        return true;
    }

    static bool innerSegmentsCross(const pointType& A, const pointType& B, const pointType& P, const pointType& Q, int xyz)
    {
        int o11, o12, o21, o22;

        if (xyz == 2)
        {
            o11 = orient2Dxy(P, A, B);
            o12 = orient2Dxy(Q, B, A);
            o21 = orient2Dxy(A, P, Q);
            o22 = orient2Dxy(B, Q, P);
        }
        else if (xyz == 0)
        {
            o11 = orient2Dyz(P, A, B);
            o12 = orient2Dyz(Q, B, A);
            o21 = orient2Dyz(A, P, Q);
            o22 = orient2Dyz(B, Q, P);
        }
        else
        {
            o11 = orient2Dzx(P, A, B);
            o12 = orient2Dzx(Q, B, A);
            o21 = orient2Dzx(A, P, Q);
            o22 = orient2Dzx(B, Q, P);
        }

        return (o11 && o21 && o11 == o12 && o21 == o22);
    }

    static bool pointInInnerSegment(const pointType& p, const pointType& v1, const pointType& v2)
    {
        if (misaligned(p, v1, v2)) return false;

        int lt2, lt3;
        lt2 = v1[0] < p[0];
        lt3 = p[0] < v2[0];
        if (lt2) return (lt2 == lt3);
        lt2 = v1[1] < p[1];
        lt3 = p[1] < v2[1];
        if (lt2) return (lt2 == lt3);
        lt2 = v1[2] < p[2];
        lt3 = p[2] < v2[2];
        if (lt2) return (lt2 == lt3);
        return false;
    }

    static inline bool innerSegmentsCross(const pointType& A, const pointType& B, const pointType& P, const pointType& Q)
    {
        int o11, o12, o21, o22;

        o11 = orient2Dxy(P, A, B);
        o12 = orient2Dxy(Q, B, A);
        o21 = orient2Dxy(A, P, Q);
        o22 = orient2Dxy(B, Q, P);
        if (o11 || o21 || o12 || o22) return (o11 == o12 && o21 == o22);

        o11 = orient2Dyz(P, A, B);
        o12 = orient2Dyz(Q, B, A);
        o21 = orient2Dyz(A, P, Q);
        o22 = orient2Dyz(B, Q, P);
        if (o11 || o21 || o12 || o22) return (o11 == o12 && o21 == o22);

        o11 = orient2Dzx(P, A, B);
        o12 = orient2Dzx(Q, B, A);
        o21 = orient2Dzx(A, P, Q);
        o22 = orient2Dzx(B, Q, P);
        if (o11 || o21 || o12 || o22) return (o11 == o12 && o21 == o22);

        return false;
    }

    static inline bool lineCrossesInnerTriangle(const pointType& s1, const pointType& s2, const pointType& v1, const pointType& v2, const pointType& v3)
    {
        const int o1 = pointType::orient3D(s1, s2, v1, v2);
        const int o2 = pointType::orient3D(s1, s2, v2, v3);
        if ((o1 >= 0 && o2 <= 0) || (o1 <= 0 && o2 >= 0)) return false;
        const int o3 = pointType::orient3D(s1, s2, v3, v1);
        if ((o1 >= 0 && o3 <= 0) || (o1 <= 0 && o3 >= 0)) return false;
        if ((o2 >= 0 && o3 <= 0) || (o2 <= 0 && o3 >= 0)) return false;
        return true;
    }

    static inline bool innerSegmentCrossesInnerTriangle(const pointType& s1, const pointType& s2, const pointType& v1, const pointType& v2, const pointType& v3)
    {
        int o1 = orient3D(s1, v1, v2, v3); if (o1 == 0) return false;
        int o2 = orient3D(s2, v1, v2, v3); if (o2 == 0) return false;

        if ((o1 > 0 && o2 > 0) || (o1 < 0 && o2 < 0)) return false;
        o1 = orient3D(s1, s2, v1, v2);
        o2 = orient3D(s1, s2, v2, v3);
        if ((o1 >= 0 && o2 <= 0) || (o1 <= 0 && o2 >= 0)) return false;
        int o3 = orient3D(s1, s2, v3, v1);
        if ((o1 >= 0 && o3 <= 0) || (o1 <= 0 && o3 >= 0)) return false;
        if ((o2 >= 0 && o3 <= 0) || (o2 <= 0 && o3 >= 0)) return false;
        return true;
    }

    static inline bool pointInInnerTriangle(const pointType& P, const pointType& A, const pointType& B, const pointType& C)
    {
        int o1, o2, o3;
        o1 = orient2Dxy(P, A, B);
        o2 = orient2Dxy(P, B, C);
        o3 = orient2Dxy(P, C, A);
        if (o1 || o2 || o3) return ((o1 > 0 && o2 > 0 && o3 > 0) || (o1 < 0 && o2 < 0 && o3 < 0));
        o1 = orient2Dyz(P, A, B);
        o2 = orient2Dyz(P, B, C);
        o3 = orient2Dyz(P, C, A);
        if (o1 || o2 || o3) return ((o1 > 0 && o2 > 0 && o3 > 0) || (o1 < 0 && o2 < 0 && o3 < 0));
        o1 = orient2Dzx(P, A, B);
        o2 = orient2Dzx(P, B, C);
        o3 = orient2Dzx(P, C, A);
        return ((o1 > 0 && o2 > 0 && o3 > 0) || (o1 < 0 && o2 < 0 && o3 < 0));
    }

};

typedef pointType explicitPoint;
typedef pointType implicitPoint_LNC;

inline std::ostream& operator<<(std::ostream& os, const pointType& p)
{
    return os << GET_DOUBLE_VAL(p[0]) << " " << GET_DOUBLE_VAL(p[1]) << " " << GET_DOUBLE_VAL(p[2]);
}

#endif // USE_INDIRECT_PREDS

#endif
