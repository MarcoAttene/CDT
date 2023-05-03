#include "implicit_point.h"

int orient3d_LEEE(double p1x, double p1y, double p1z, double p2x, double p2y, double p2z, double pt, double qx, double qy, double qz, double rx, double ry, double rz, double sx, double sy, double sz);
int orient3d_LLEE(double p1x, double p1y, double p1z, double p2x, double p2y, double p2z, double pt, double q1x, double q1y, double q1z, double q2x, double q2y, double q2z, double qt, double rx, double ry, double rz, double sx, double sy, double sz);
int orient3d_LLLE(double p1x, double p1y, double p1z, double p2x, double p2y, double p2z, double pt, double q1x, double q1y, double q1z, double q2x, double q2y, double q2z, double qt, double r1x, double r1y, double r1z, double r2x, double r2y, double r2z, double rt, double sx, double sy, double sz);
int orient3d_LLLL(double p1x, double p1y, double p1z, double p2x, double p2y, double p2z, double pt, double q1x, double q1y, double q1z, double q2x, double q2y, double q2z, double qt, double r1x, double r1y, double r1z, double r2x, double r2y, double r2z, double rt, double s1x, double s1y, double s1z, double s2x, double s2y, double s2z, double st);
int inSphere_LEEEE(double a1x, double a1y, double a1z, double a2x, double a2y, double a2z, double at, double bx, double by, double bz, double cx, double cy, double cz, double dx, double dy, double dz, double ex, double ey, double ez);
int inSphere_LLEEE(double a1x, double a1y, double a1z, double a2x, double a2y, double a2z, double at, double b1x, double b1y, double b1z, double b2x, double b2y, double b2z, double bt, double cx, double cy, double cz, double dx, double dy, double dz, double ex, double ey, double ez);
int inSphere_LLLEE(double a1x, double a1y, double a1z, double a2x, double a2y, double a2z, double at, double b1x, double b1y, double b1z, double b2x, double b2y, double b2z, double bt, double c1x, double c1y, double c1z, double c2x, double c2y, double c2z, double ct, double dx, double dy, double dz, double ex, double ey, double ez);
int inSphere_LLLLE(double a1x, double a1y, double a1z, double a2x, double a2y, double a2z, double at, double b1x, double b1y, double b1z, double b2x, double b2y, double b2z, double bt, double c1x, double c1y, double c1z, double c2x, double c2y, double c2z, double ct, double d1x, double d1y, double d1z, double d2x, double d2y, double d2z, double dt, double ex, double ey, double ez);
int inSphere_LLLLL(double a1x, double a1y, double a1z, double a2x, double a2y, double a2z, double at, double b1x, double b1y, double b1z, double b2x, double b2y, double b2z, double bt, double c1x, double c1y, double c1z, double c2x, double c2y, double c2z, double ct, double d1x, double d1y, double d1z, double d2x, double d2y, double d2z, double dt, double e1x, double e1y, double e1z, double e2x, double e2y, double e2z, double et);

inline int lnc_orient3d_EEEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	const explicitPoint3D &ae = a.toExplicit3D(), &be = b.toExplicit3D(), &ce = c.toExplicit3D(), &de = d.toExplicit3D();
	return orient3d(ae.X(), ae.Y(), ae.Z(), be.X(), be.Y(), be.Z(), ce.X(), ce.Y(), ce.Z(), de.X(), de.Y(), de.Z());
}

inline int lnc_orient3d_IEEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	const implicitPoint3D_LNC& al = a.toLNC();
	const explicitPoint3D& be = b.toExplicit3D(), & ce = c.toExplicit3D(), & de = d.toExplicit3D();
	return orient3d_LEEE(al.P().X(), al.P().Y(), al.P().Z(), al.Q().X(), al.Q().Y(), al.Q().Z(), al.T(), 
		be.X(), be.Y(), be.Z(), ce.X(), ce.Y(), ce.Z(), de.X(), de.Y(), de.Z());
}

inline int lnc_orient3d_IIEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	const implicitPoint3D_LNC& al = a.toLNC(), &bl = b.toLNC();
	const explicitPoint3D& ce = c.toExplicit3D(), & de = d.toExplicit3D();
	return orient3d_LLEE(al.P().X(), al.P().Y(), al.P().Z(), al.Q().X(), al.Q().Y(), al.Q().Z(), al.T(), 
		bl.P().X(), bl.P().Y(), bl.P().Z(), bl.Q().X(), bl.Q().Y(), bl.Q().Z(), bl.T(), 
		ce.X(), ce.Y(), ce.Z(), de.X(), de.Y(), de.Z());
}

inline int lnc_orient3d_IIIE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	const implicitPoint3D_LNC& al = a.toLNC(), & bl = b.toLNC(), & cl = c.toLNC();
	const explicitPoint3D& de = d.toExplicit3D();
	return orient3d_LLLE(al.P().X(), al.P().Y(), al.P().Z(), al.Q().X(), al.Q().Y(), al.Q().Z(), al.T(),
		bl.P().X(), bl.P().Y(), bl.P().Z(), bl.Q().X(), bl.Q().Y(), bl.Q().Z(), bl.T(),
		cl.P().X(), cl.P().Y(), cl.P().Z(), cl.Q().X(), cl.Q().Y(), cl.Q().Z(), cl.T(),
		de.X(), de.Y(), de.Z());
}

inline int lnc_orient3d_IIII(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	const implicitPoint3D_LNC& al = a.toLNC(), & bl = b.toLNC(), & cl = c.toLNC(), & dl = d.toLNC();
	return orient3d_LLLL(al.P().X(), al.P().Y(), al.P().Z(), al.Q().X(), al.Q().Y(), al.Q().Z(), al.T(),
		bl.P().X(), bl.P().Y(), bl.P().Z(), bl.Q().X(), bl.Q().Y(), bl.Q().Z(), bl.T(),
		cl.P().X(), cl.P().Y(), cl.P().Z(), cl.Q().X(), cl.Q().Y(), cl.Q().Z(), cl.T(),
		dl.P().X(), dl.P().Y(), dl.P().Z(), dl.Q().X(), dl.Q().Y(), dl.Q().Z(), dl.T());
}

// This version assumes that points are either explicit3D or LNC
inline int lnc_orient3D(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	const int i = a.isExplicit3D() + b.isExplicit3D() + c.isExplicit3D() + d.isExplicit3D();

	if (i == 4) return lnc_orient3d_EEEE(a, b, c, d);

	if (i == 3)
	{
		if (!a.isExplicit3D()) return lnc_orient3d_IEEE(a, b, c, d);
		if (!b.isExplicit3D()) return lnc_orient3d_IEEE(b, c, a, d);
		if (!c.isExplicit3D()) return lnc_orient3d_IEEE(c, d, a, b);
		return lnc_orient3d_IEEE(d, a, c, b);
	}

	if (i == 2)
	{
		if (c.isExplicit3D() && d.isExplicit3D()) return lnc_orient3d_IIEE(a, b, c, d);
		if (b.isExplicit3D() && d.isExplicit3D()) return lnc_orient3d_IIEE(a, c, d, b);
		if (a.isExplicit3D() && d.isExplicit3D()) return lnc_orient3d_IIEE(b, c, a, d);
		if (b.isExplicit3D() && c.isExplicit3D()) return lnc_orient3d_IIEE(d, a, c, b);
		if (a.isExplicit3D() && c.isExplicit3D()) return lnc_orient3d_IIEE(d, b, a, c);
		return lnc_orient3d_IIEE(c, d, a, b);
	}

	if (i == 1)
	{
		if (d.isExplicit3D()) return lnc_orient3d_IIIE(a, b, c, d);
		if (c.isExplicit3D()) return lnc_orient3d_IIIE(d, b, a, c);
		if (b.isExplicit3D()) return lnc_orient3d_IIIE(a, c, d, b);
		return lnc_orient3d_IIIE(b, d, c, a);
	}

	return lnc_orient3d_IIII(a, b, c, d);
}



inline int lnc_inSphere_IEEEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d, const genericPoint& e) {
	const implicitPoint3D_LNC& al = a.toLNC();
	const explicitPoint3D& be = b.toExplicit3D(), & ce = c.toExplicit3D(), & de = d.toExplicit3D(), & ee = e.toExplicit3D();
	return inSphere_LEEEE(al.P().X(), al.P().Y(), al.P().Z(), al.Q().X(), al.Q().Y(), al.Q().Z(), al.T(),
		be.X(), be.Y(), be.Z(), ce.X(), ce.Y(), ce.Z(), de.X(), de.Y(), de.Z(), ee.X(), ee.Y(), ee.Z());
}

inline int lnc_inSphere_IIEEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d, const genericPoint& e) {
	const implicitPoint3D_LNC& al = a.toLNC(), & bl = b.toLNC();
	const explicitPoint3D& ce = c.toExplicit3D(), & de = d.toExplicit3D(), & ee = e.toExplicit3D();
	return inSphere_LLEEE(al.P().X(), al.P().Y(), al.P().Z(), al.Q().X(), al.Q().Y(), al.Q().Z(), al.T(),
		bl.P().X(), bl.P().Y(), bl.P().Z(), bl.Q().X(), bl.Q().Y(), bl.Q().Z(), bl.T(),
		ce.X(), ce.Y(), ce.Z(), de.X(), de.Y(), de.Z(), ee.X(), ee.Y(), ee.Z());
}

inline int lnc_inSphere_IIIEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d, const genericPoint& e) {
	const implicitPoint3D_LNC& al = a.toLNC(), & bl = b.toLNC(), & cl = c.toLNC();
	const explicitPoint3D& de = d.toExplicit3D(), & ee = e.toExplicit3D();
	return inSphere_LLLEE(al.P().X(), al.P().Y(), al.P().Z(), al.Q().X(), al.Q().Y(), al.Q().Z(), al.T(),
		bl.P().X(), bl.P().Y(), bl.P().Z(), bl.Q().X(), bl.Q().Y(), bl.Q().Z(), bl.T(),
		cl.P().X(), cl.P().Y(), cl.P().Z(), cl.Q().X(), cl.Q().Y(), cl.Q().Z(), cl.T(),
		de.X(), de.Y(), de.Z(), ee.X(), ee.Y(), ee.Z());
}

inline int lnc_inSphere_IIIIE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d, const genericPoint& e) {
	const implicitPoint3D_LNC& al = a.toLNC(), & bl = b.toLNC(), & cl = c.toLNC(), & dl = d.toLNC();
	const explicitPoint3D& ee = e.toExplicit3D();
	return inSphere_LLLLE(al.P().X(), al.P().Y(), al.P().Z(), al.Q().X(), al.Q().Y(), al.Q().Z(), al.T(),
		bl.P().X(), bl.P().Y(), bl.P().Z(), bl.Q().X(), bl.Q().Y(), bl.Q().Z(), bl.T(),
		cl.P().X(), cl.P().Y(), cl.P().Z(), cl.Q().X(), cl.Q().Y(), cl.Q().Z(), cl.T(),
		dl.P().X(), dl.P().Y(), dl.P().Z(), dl.Q().X(), dl.Q().Y(), dl.Q().Z(), dl.T(),
		ee.X(), ee.Y(), ee.Z());
}

inline int lnc_inSphere_IIIII(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d, const genericPoint& e) {
	const implicitPoint3D_LNC& al = a.toLNC(), & bl = b.toLNC(), & cl = c.toLNC(), & dl = d.toLNC(), & el = e.toLNC();
	return inSphere_LLLLL(al.P().X(), al.P().Y(), al.P().Z(), al.Q().X(), al.Q().Y(), al.Q().Z(), al.T(),
		bl.P().X(), bl.P().Y(), bl.P().Z(), bl.Q().X(), bl.Q().Y(), bl.Q().Z(), bl.T(),
		cl.P().X(), cl.P().Y(), cl.P().Z(), cl.Q().X(), cl.Q().Y(), cl.Q().Z(), cl.T(),
		dl.P().X(), dl.P().Y(), dl.P().Z(), dl.Q().X(), dl.Q().Y(), dl.Q().Z(), dl.T(),
		el.P().X(), el.P().Y(), el.P().Z(), el.Q().X(), el.Q().Y(), el.Q().Z(), el.T());
}

inline int lnc_inSphere(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d, const genericPoint& e)
{
	const int num_explicit = a.isExplicit3D() + b.isExplicit3D() + c.isExplicit3D() + d.isExplicit3D() + e.isExplicit3D();

	if (num_explicit == 5) return ::inSphere(a.toExplicit3D().X(), a.toExplicit3D().Y(), a.toExplicit3D().Z(),
		b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z(),
		c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z(),
		d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z(),
		e.toExplicit3D().X(), e.toExplicit3D().Y(), e.toExplicit3D().Z());

	const genericPoint* A[5] = { &a, &b, &c, &d, &e };

	// Sort points so that I < E
	bool swapped = true;
	int sign_swap = 1;

	while (swapped) {
		swapped = false;
		for (int i = 0; i < 4; i++) {
			if (A[i]->isExplicit3D() && !A[i + 1]->isExplicit3D()) {
				std::swap(A[i], A[i + 1]);
				swapped = true;
				sign_swap *= -1;
			}
		}
	}

	if (num_explicit == 4) return sign_swap * lnc_inSphere_IEEEE(*A[0], *A[1], *A[2], *A[3], *A[4]);
	if (num_explicit == 3) return sign_swap * lnc_inSphere_IIEEE(*A[0], *A[1], *A[2], *A[3], *A[4]);
	if (num_explicit == 2) return sign_swap * lnc_inSphere_IIIEE(*A[0], *A[1], *A[2], *A[3], *A[4]);
	if (num_explicit == 1) return sign_swap * lnc_inSphere_IIIIE(*A[0], *A[1], *A[2], *A[3], *A[4]);
	return sign_swap * lnc_inSphere_IIIII(*A[0], *A[1], *A[2], *A[3], *A[4]);
}

#include "lnc_predicates.hpp"
