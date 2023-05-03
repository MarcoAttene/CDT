/****************************************************************************
* Indirect predicates for geometric constructions					        *
*                                                                           *
* Consiglio Nazionale delle Ricerche                                        *
* Istituto di Matematica Applicata e Tecnologie Informatiche                *
* Sezione di Genova                                                         *
* IMATI-GE / CNR                                                            *
*                                                                           *
* Authors: Marco Attene                                                     *
* Copyright(C) 2019: IMATI-GE / CNR                                         *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU Lesser General Public License as published  *
* by the Free Software Foundation; either version 3 of the License, or (at  *
* your option) any later version.                                           *
*                                                                           *
* This program is distributed in the hope that it will be useful, but       *
* WITHOUT ANY WARRANTY; without even the implied warranty of                *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser  *
* General Public License for more details.                                  *
*                                                                           *
* You should have received a copy of the GNU Lesser General Public License  *
* along with this program.  If not, see http://www.gnu.org/licenses/.       *
*                                                                           *
****************************************************************************/

#include "implicit_point.h"
#include "indirect_predicates.h"

int orient2d(double p1x, double p1y, double p2x, double p2y, double p3x, double p3y);
int orient3d(double px, double py, double pz, double qx, double qy, double qz, double rx, double ry, double rz, double sx, double sy, double sz);

inline int lessThan_EE(double x1, double y1, double z1, double x2, double y2, double z2)
{
	int ret;
	if ((ret = ((x1 > x2) - (x1 < x2)))) return ret;
	if ((ret = ((y1 > y2) - (y1 < y2)))) return ret;
	return ((z1 > z2) - (z1 < z2));
}

inline int lessThan_IE(const genericPoint& p1, double x2, double y2, double z2)
{
	int ret;
	if ((ret = lessThanOnX_IE(p1, x2))) return ret;
	if ((ret = lessThanOnY_IE(p1, y2))) return ret;
	return lessThanOnZ_IE(p1, z2);
}

inline int lessThan_II(const genericPoint& p1, const genericPoint& p2)
{
	int ret;
	if ((ret = lessThanOnX_II(p1, p2))) return ret;
	if ((ret = lessThanOnY_II(p1, p2))) return ret;
	return lessThanOnZ_II(p1, p2);
}

inline int lessThan_EE(const genericPoint& a, const genericPoint& b) { return lessThan_EE(a.toExplicit3D().X(), a.toExplicit3D().Y(), a.toExplicit3D().Z(), b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z()); }
inline int lessThan_IE(const genericPoint& a, const genericPoint& b) { return lessThan_IE(a, b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z()); }

inline int genericPoint::lessThan(const genericPoint& a, const genericPoint& b)
{
		if (a.isExplicit3D() && b.isExplicit3D()) return lessThan_EE(a, b);
		if (!a.isExplicit3D() && b.isExplicit3D()) return lessThan_IE(a, b);
		if (a.isExplicit3D() && !b.isExplicit3D()) return -lessThan_IE(b, a);
		return lessThan_II(a, b);
}

inline int lessThanOnX_IE(const genericPoint& a, const genericPoint& b) { return lessThanOnX_IE(a, b.toExplicit3D().X()); }

inline int genericPoint::lessThanOnX(const genericPoint& a, const genericPoint& b)
{
	if (a.isExplicit3D() && b.isExplicit3D()) return a.toExplicit3D().X() < b.toExplicit3D().X();
	if (!a.isExplicit3D() && b.isExplicit3D()) return lessThanOnX_IE(a, b);
	if (a.isExplicit3D() && !b.isExplicit3D()) return -lessThanOnX_IE(b, a);
	return lessThanOnX_II(a, b);
}

inline int lessThanOnY_IE(const genericPoint& a, const genericPoint& b) { return lessThanOnY_IE(a, b.toExplicit3D().Y()); }

inline int genericPoint::lessThanOnY(const genericPoint& a, const genericPoint& b)
{
	if (a.isExplicit3D() && b.isExplicit3D()) return a.toExplicit3D().Y() < b.toExplicit3D().Y();
	if (!a.isExplicit3D() && b.isExplicit3D()) return lessThanOnY_IE(a, b);
	if (a.isExplicit3D() && !b.isExplicit3D()) return -lessThanOnY_IE(b, a);
	return lessThanOnY_II(a, b);
}

inline int lessThanOnZ_IE(const genericPoint& a, const genericPoint& b) { return lessThanOnZ_IE(a, b.toExplicit3D().Z()); }

inline int genericPoint::lessThanOnZ(const genericPoint& a, const genericPoint& b)
{
	if (a.isExplicit3D() && b.isExplicit3D()) return a.toExplicit3D().Z() < b.toExplicit3D().Z();
	if (!a.isExplicit3D() && b.isExplicit3D()) return lessThanOnZ_IE(a, b);
	if (a.isExplicit3D() && !b.isExplicit3D()) return -lessThanOnZ_IE(b, a);
	return lessThanOnZ_II(a, b);
}

inline int dotproductSign2D_EEE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return dotProductSign2D(a.toExplicit2D().X(), a.toExplicit2D().Y(), b.toExplicit2D().X(), b.toExplicit2D().Y(), c.toExplicit2D().X(), c.toExplicit2D().Y()); }
inline int dotproductSign2D_IEE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return dotProductSign2D_IEE(a, b.toExplicit2D().X(), b.toExplicit2D().Y(), c.toExplicit2D().X(), c.toExplicit2D().Y()); }
inline int dotproductSign2D_IIE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return dotProductSign2D_IIE(a, b, c.toExplicit2D().X(), c.toExplicit2D().Y()); }
inline int dotproductSign2D_III(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return dotProductSign2D_III(a, b, c); }
inline int dotproductSign2D_EEI(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return dotProductSign2D_EEI(c, a.toExplicit2D().X(), a.toExplicit2D().Y(), b.toExplicit2D().X(), b.toExplicit2D().Y()); }
inline int dotproductSign2D_IEI(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return dotProductSign2D_IEI(a, c, b.toExplicit2D().X(), b.toExplicit2D().Y()); }

inline int genericPoint::dotProductSign2D(const genericPoint& a, const genericPoint& b, const genericPoint& c)
{
	if (a.isExplicit2D() && b.isExplicit2D() && c.isExplicit2D()) return dotproductSign2D_EEE(a, b, c);
	if (a.isSSI() && b.isExplicit2D() && c.isExplicit2D()) return dotproductSign2D_IEE(a, b, c);
	if (a.isExplicit2D() && b.isSSI() && c.isExplicit2D()) return dotproductSign2D_IEE(b, a, c);
	if (a.isExplicit2D() && b.isExplicit2D() && c.isSSI()) return dotproductSign2D_EEI(a, b, c);
	if (a.isSSI() && b.isSSI() && c.isExplicit2D()) return dotproductSign2D_IIE(a, b, c);
	if (a.isExplicit2D() && b.isSSI() && c.isSSI())	return dotproductSign2D_IEI(b, a, c);
	if (a.isSSI() && b.isExplicit2D() && c.isSSI())	return dotproductSign2D_IEI(a, b, c);
	return dotproductSign2D_III(a, b, c);
}

inline int dotproductSign3D_EEE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return dotProductSign3D(a.toExplicit3D().X(), a.toExplicit3D().Y(), a.toExplicit3D().Z(), b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z(), c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z()); }
inline int dotproductSign3D_IEE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return dotProductSign3D_IEE(a, b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z(), c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z()); }
inline int dotproductSign3D_IIE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return dotProductSign3D_IIE(a, b, c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z()); }
inline int dotproductSign3D_III(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return dotProductSign3D_III(a, b, c); }
inline int dotproductSign3D_EEI(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return dotProductSign3D_EEI(c, a.toExplicit3D().X(), a.toExplicit3D().Y(), a.toExplicit3D().Z(), b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z()); }
inline int dotproductSign3D_IEI(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return dotProductSign3D_IEI(a, c, b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z()); }

inline int genericPoint::dotProductSign3D(const genericPoint& a, const genericPoint& b, const genericPoint& c)
{
	if (a.isExplicit3D() && b.isExplicit3D() && c.isExplicit3D()) return dotproductSign3D_EEE(a, b, c);
	if (!a.isExplicit3D() && b.isExplicit3D() && c.isExplicit3D()) return dotproductSign3D_IEE(a, b, c);
	if (a.isExplicit3D() && !b.isExplicit3D() && c.isExplicit3D()) return dotproductSign3D_IEE(b, a, c);
	if (a.isExplicit3D() && b.isExplicit3D() && !c.isExplicit3D()) return dotproductSign3D_EEI(a, b, c);
	if (!a.isExplicit3D() && !b.isExplicit3D() && c.isExplicit3D()) return dotproductSign3D_IIE(a, b, c);
	if (a.isExplicit3D() && !b.isExplicit3D() && !c.isExplicit3D())	return dotproductSign3D_IEI(b, a, c);
	if (!a.isExplicit3D() && b.isExplicit3D() && !c.isExplicit3D())	return dotproductSign3D_IEI(a, b, c);
	return dotproductSign3D_III(a, b, c);
}

inline int orient2d_EEE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2d(a.toExplicit2D().X(), a.toExplicit2D().Y(), b.toExplicit2D().X(), b.toExplicit2D().Y(), c.toExplicit2D().X(), c.toExplicit2D().Y()); }
inline int orient2d_IEE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2d_indirect_IEE(a, b.toExplicit2D().X(), b.toExplicit2D().Y(), c.toExplicit2D().X(), c.toExplicit2D().Y()); }
inline int orient2d_IIE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2d_indirect_IIE(a, b, c.toExplicit2D().X(), c.toExplicit2D().Y()); }
inline int orient2d_III(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2d_indirect_III(a, b, c); }

inline int genericPoint::orient2D(const genericPoint& a, const genericPoint& b, const genericPoint& c)
{
	// Here we implicitly assume that points are 2D. Do not check.

	if (a.isExplicit2D() && b.isExplicit2D() && c.isExplicit2D()) return orient2d_EEE(a, b, c);
	if (a.isSSI() && b.isExplicit2D() && c.isExplicit2D()) return orient2d_IEE(a, b, c);
	if (a.isExplicit2D() && b.isSSI() && c.isExplicit2D()) return orient2d_IEE(b, c, a);
	if (a.isExplicit2D() && b.isExplicit2D() && c.isSSI()) return orient2d_IEE(c, a, b);
	if (a.isSSI() && b.isSSI() && c.isExplicit2D()) return orient2d_IIE(a, b, c);
	if (a.isExplicit2D() && b.isSSI() && c.isSSI())	return orient2d_IIE(b, c, a);
	if (a.isSSI() && b.isExplicit2D() && c.isSSI())	return orient2d_IIE(c, a, b);
	return orient2d_III(a, b, c);
}


inline int orient2dxy_EEE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2d(a.toExplicit3D().X(), a.toExplicit3D().Y(), b.toExplicit3D().X(), b.toExplicit3D().Y(), c.toExplicit3D().X(), c.toExplicit3D().Y()); }
inline int orient2dxy_IEE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2dxy_indirect_IEE(a, b.toExplicit3D().X(), b.toExplicit3D().Y(), c.toExplicit3D().X(), c.toExplicit3D().Y()); }
inline int orient2dxy_IIE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2dxy_indirect_IIE(a, b, c.toExplicit3D().X(), c.toExplicit3D().Y()); }
inline int orient2dxy_III(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2dxy_indirect_III(a, b, c); }

inline int genericPoint::orient2Dxy(const genericPoint& a, const genericPoint& b, const genericPoint& c)
{
	if (a.isExplicit3D() && b.isExplicit3D() && c.isExplicit3D()) return orient2dxy_EEE(a, b, c);

	if (!a.isExplicit3D() && b.isExplicit3D() && c.isExplicit3D()) return orient2dxy_IEE(a, b, c);
	if (a.isExplicit3D() && !b.isExplicit3D() && c.isExplicit3D()) return orient2dxy_IEE(b, c, a);
	if (a.isExplicit3D() && b.isExplicit3D() && !c.isExplicit3D()) return orient2dxy_IEE(c, a, b);

	if (!a.isExplicit3D() && !b.isExplicit3D() && c.isExplicit3D()) return orient2dxy_IIE(a, b, c);
	if (!a.isExplicit3D() && b.isExplicit3D() && !c.isExplicit3D()) return orient2dxy_IIE(c, a, b);
	if (a.isExplicit3D() && !b.isExplicit3D() && !c.isExplicit3D()) return orient2dxy_IIE(b, c, a);

	return orient2dxy_III(a, b, c);
}


inline int orient2dyz_EEE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2d(a.toExplicit3D().Y(), a.toExplicit3D().Z(), b.toExplicit3D().Y(), b.toExplicit3D().Z(), c.toExplicit3D().Y(), c.toExplicit3D().Z()); }
inline int orient2dyz_IEE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2dyz_indirect_IEE(a, b.toExplicit3D().Y(), b.toExplicit3D().Z(), c.toExplicit3D().Y(), c.toExplicit3D().Z()); }
inline int orient2dyz_IIE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2dyz_indirect_IIE(a, b, c.toExplicit3D().Y(), c.toExplicit3D().Z()); }
inline int orient2dyz_III(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2dyz_indirect_III(a, b, c); }

inline int genericPoint::orient2Dyz(const genericPoint& a, const genericPoint& b, const genericPoint& c)
{
	if (a.isExplicit3D() && b.isExplicit3D() && c.isExplicit3D()) return orient2dyz_EEE(a, b, c);

	if (!a.isExplicit3D() && b.isExplicit3D() && c.isExplicit3D()) return orient2dyz_IEE(a, b, c);
	if (a.isExplicit3D() && !b.isExplicit3D() && c.isExplicit3D()) return orient2dyz_IEE(b, c, a);
	if (a.isExplicit3D() && b.isExplicit3D() && !c.isExplicit3D()) return orient2dyz_IEE(c, a, b);

	if (!a.isExplicit3D() && !b.isExplicit3D() && c.isExplicit3D()) return orient2dyz_IIE(a, b, c);
	if (!a.isExplicit3D() && b.isExplicit3D() && !c.isExplicit3D()) return orient2dyz_IIE(c, a, b);
	if (a.isExplicit3D() && !b.isExplicit3D() && !c.isExplicit3D()) return orient2dyz_IIE(b, c, a);

	return orient2dyz_III(a, b, c);
}


inline int orient2dzx_EEE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2d(a.toExplicit3D().Z(), a.toExplicit3D().X(), b.toExplicit3D().Z(), b.toExplicit3D().X(), c.toExplicit3D().Z(), c.toExplicit3D().X()); }
inline int orient2dzx_IEE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2dzx_indirect_IEE(a, b.toExplicit3D().Z(), b.toExplicit3D().X(), c.toExplicit3D().Z(), c.toExplicit3D().X()); }
inline int orient2dzx_IIE(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2dzx_indirect_IIE(a, b, c.toExplicit3D().Z(), c.toExplicit3D().X()); }
inline int orient2dzx_III(const genericPoint& a, const genericPoint& b, const genericPoint& c) { return orient2dzx_indirect_III(a, b, c); }

inline int genericPoint::orient2Dzx(const genericPoint& a, const genericPoint& b, const genericPoint& c)
{
	if (a.isExplicit3D() && b.isExplicit3D() && c.isExplicit3D()) return orient2dzx_EEE(a, b, c);

	if (!a.isExplicit3D() && b.isExplicit3D() && c.isExplicit3D()) return orient2dzx_IEE(a, b, c);
	if (a.isExplicit3D() && !b.isExplicit3D() && c.isExplicit3D()) return orient2dzx_IEE(b, c, a);
	if (a.isExplicit3D() && b.isExplicit3D() && !c.isExplicit3D()) return orient2dzx_IEE(c, a, b);

	if (!a.isExplicit3D() && !b.isExplicit3D() && c.isExplicit3D()) return orient2dzx_IIE(a, b, c);
	if (!a.isExplicit3D() && b.isExplicit3D() && !c.isExplicit3D()) return orient2dzx_IIE(c, a, b);
	if (a.isExplicit3D() && !b.isExplicit3D() && !c.isExplicit3D()) return orient2dzx_IIE(b, c, a);

	return orient2dzx_III(a, b, c);
}


inline int orient3d_EEEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d) 
{ 
	return orient3d(a.toExplicit3D().X(), a.toExplicit3D().Y(), a.toExplicit3D().Z(), b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z(),
	c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z(), d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z());
}

inline int orient3d_IEEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return orient3d_indirect_IEEE(a, b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z(),
		c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z(), d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z());
}

inline int orient3d_IIEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return orient3d_indirect_IIEE(a, b,
		c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z(), d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z());
}

inline int orient3d_IIIE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return orient3d_indirect_IIIE(a, b, c, d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z());
}


inline int genericPoint::orient3D(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	// Here we implicitly assume that points are 3D. Do not check.

	const int i = a.isExplicit3D() + b.isExplicit3D() + c.isExplicit3D() + d.isExplicit3D();

	if (i == 4) return orient3d_EEEE(a, b, c, d);
	
	if (i == 3)
	{
		if (!a.isExplicit3D()) return orient3d_IEEE(a, b, c, d);
		if (!b.isExplicit3D()) return orient3d_IEEE(b, c, a, d);
		if (!c.isExplicit3D()) return orient3d_IEEE(c, d, a, b);
		return orient3d_IEEE(d, a, c, b);
	}

	if (i == 2)
	{
		if (c.isExplicit3D() && d.isExplicit3D()) return orient3d_IIEE(a, b, c, d);
		if (b.isExplicit3D() && d.isExplicit3D()) return orient3d_IIEE(a, c, d, b);
		if (a.isExplicit3D() && d.isExplicit3D()) return orient3d_IIEE(b, c, a, d);
		if (b.isExplicit3D() && c.isExplicit3D()) return orient3d_IIEE(d, a, c, b);
		if (a.isExplicit3D() && c.isExplicit3D()) return orient3d_IIEE(d, b, a, c);
		return orient3d_IIEE(c, d, a, b);
	}

	if (i == 1)
	{
		if (d.isExplicit3D()) return orient3d_IIIE(a, b, c, d);
		if (c.isExplicit3D()) return orient3d_IIIE(d, b, a, c);
		if (b.isExplicit3D()) return orient3d_IIIE(a, c, d, b);
		return orient3d_IIIE(b, d, c, a);
	}

	return orient3d_indirect_IIII(a, b, c, d);
}

inline int inSphere_IEEEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d, const genericPoint& e) {
	return inSphere_IEEEE(a,
		b.toExplicit3D().X(), b.toExplicit3D().Y(), b.toExplicit3D().Z(),
		c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z(),
		d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z(),
		e.toExplicit3D().X(), e.toExplicit3D().Y(), e.toExplicit3D().Z());
}

inline int inSphere_IIEEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d, const genericPoint& e) {
	return inSphere_IIEEE(a, b,
		c.toExplicit3D().X(), c.toExplicit3D().Y(), c.toExplicit3D().Z(),
		d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z(),
		e.toExplicit3D().X(), e.toExplicit3D().Y(), e.toExplicit3D().Z());
}

inline int inSphere_IIIEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d, const genericPoint& e) {
	return inSphere_IIIEE(a, b, c,
		d.toExplicit3D().X(), d.toExplicit3D().Y(), d.toExplicit3D().Z(),
		e.toExplicit3D().X(), e.toExplicit3D().Y(), e.toExplicit3D().Z());
}

inline int inSphere_IIIIE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d, const genericPoint& e) {
	return inSphere_IIIIE(a, b, c, d,
		e.toExplicit3D().X(), e.toExplicit3D().Y(), e.toExplicit3D().Z());
}

inline int genericPoint::inSphere(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d, const genericPoint& e)
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

	if (num_explicit == 4) return sign_swap * inSphere_IEEEE(*A[0], *A[1], *A[2], *A[3], *A[4]);
	if (num_explicit == 3) return sign_swap * inSphere_IIEEE(*A[0], *A[1], *A[2], *A[3], *A[4]);
	if (num_explicit == 2) return sign_swap * inSphere_IIIEE(*A[0], *A[1], *A[2], *A[3], *A[4]);
	if (num_explicit == 1) return sign_swap * inSphere_IIIIE(*A[0], *A[1], *A[2], *A[3], *A[4]);
	return sign_swap * inSphere_IIIII(*A[0], *A[1], *A[2], *A[3], *A[4]);
}

inline int incircle2d_EEEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return incircle(a.toExplicit2D().X(), a.toExplicit2D().Y(), b.toExplicit2D().X(), b.toExplicit2D().Y(), c.toExplicit2D().X(), c.toExplicit2D().Y(), d.toExplicit2D().X(), d.toExplicit2D().Y());
}

inline int incircle2d_IEEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return incircle_indirect_IEEE(a, b.toExplicit2D().X(), b.toExplicit2D().Y(), c.toExplicit2D().X(), c.toExplicit2D().Y(), d.toExplicit2D().X(), d.toExplicit2D().Y());
}

inline int incircle2d_IIEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return incircle_indirect_IIEE(a, b, c.toExplicit2D().X(), c.toExplicit2D().Y(), d.toExplicit2D().X(), d.toExplicit2D().Y());
}

inline int incircle2d_IIIE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return incircle_indirect_IIIE(a, b, c, d.toExplicit2D().X(), d.toExplicit2D().Y());
}


inline int genericPoint::incircle(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	const int i = a.isExplicit2D() + b.isExplicit2D() + c.isExplicit2D() + d.isExplicit2D();

	if (i == 4) return incircle2d_EEEE(a, b, c, d);

	if (i == 3)
	{
		if (!a.isExplicit2D()) return incircle2d_IEEE(a, b, c, d);
		if (!b.isExplicit2D()) return incircle2d_IEEE(b, c, a, d);
		if (!c.isExplicit2D()) return incircle2d_IEEE(c, d, a, b);
		return incircle2d_IEEE(d, a, c, b);
	}

	if (i == 2)
	{
		if (c.isExplicit2D() && d.isExplicit2D()) return incircle2d_IIEE(a, b, c, d);
		if (b.isExplicit2D() && d.isExplicit2D()) return incircle2d_IIEE(a, c, d, b);
		if (a.isExplicit2D() && d.isExplicit2D()) return incircle2d_IIEE(b, c, a, d);
		if (b.isExplicit2D() && c.isExplicit2D()) return incircle2d_IIEE(d, a, c, b);
		if (a.isExplicit2D() && c.isExplicit2D()) return incircle2d_IIEE(d, b, a, c);
		return incircle2d_IIEE(c, d, a, b);
	}

	if (i == 1)
	{
		if (d.isExplicit2D()) return incircle2d_IIIE(a, b, c, d);
		if (c.isExplicit2D()) return incircle2d_IIIE(d, b, a, c);
		if (b.isExplicit2D()) return incircle2d_IIIE(a, c, d, b);
		return incircle2d_IIIE(b, d, c, a);
	}

	return incircle_indirect_IIII(a, b, c, d);
}

inline int incircle2dxy_EEEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return incircle(a.toExplicit3D().X(), a.toExplicit3D().Y(), b.toExplicit3D().X(), b.toExplicit3D().Y(), c.toExplicit3D().X(), c.toExplicit3D().Y(), d.toExplicit3D().X(), d.toExplicit3D().Y());
}

inline int incircle2dxy_IEEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return incirclexy_indirect_IEEE(a, b.toExplicit3D().X(), b.toExplicit3D().Y(), c.toExplicit3D().X(), c.toExplicit3D().Y(), d.toExplicit3D().X(), d.toExplicit3D().Y());
}

inline int incircle2dxy_IIEE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return incirclexy_indirect_IIEE(a, b, c.toExplicit3D().X(), c.toExplicit3D().Y(), d.toExplicit3D().X(), d.toExplicit3D().Y());
}

inline int incircle2dxy_IIIE(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	return incirclexy_indirect_IIIE(a, b, c, d.toExplicit3D().X(), d.toExplicit3D().Y());
}


inline int genericPoint::incirclexy(const genericPoint& a, const genericPoint& b, const genericPoint& c, const genericPoint& d)
{
	int i = a.isExplicit3D() + b.isExplicit3D() + c.isExplicit3D() + d.isExplicit3D();

	if (i == 4) return incircle2dxy_EEEE(a, b, c, d);

	if (i == 3)
	{
		if (!a.isExplicit3D()) return incircle2dxy_IEEE(a, b, c, d);
		if (!b.isExplicit3D()) return incircle2dxy_IEEE(b, c, a, d);
		if (!c.isExplicit3D()) return incircle2dxy_IEEE(c, d, a, b);
		return incircle2dxy_IEEE(d, a, c, b);
	}

	if (i == 2)
	{
		if (c.isExplicit3D() && d.isExplicit3D()) return incircle2dxy_IIEE(a, b, c, d);
		if (b.isExplicit3D() && d.isExplicit3D()) return incircle2dxy_IIEE(a, c, d, b);
		if (a.isExplicit3D() && d.isExplicit3D()) return incircle2dxy_IIEE(b, c, a, d);
		if (b.isExplicit3D() && c.isExplicit3D()) return incircle2dxy_IIEE(d, a, c, b);
		if (a.isExplicit3D() && c.isExplicit3D()) return incircle2dxy_IIEE(d, b, a, c);
		return incircle2dxy_IIEE(c, d, a, b);
	}

	if (i == 1)
	{
		if (d.isExplicit3D()) return incircle2dxy_IIIE(a, b, c, d);
		if (c.isExplicit3D()) return incircle2dxy_IIIE(d, b, a, c);
		if (b.isExplicit3D()) return incircle2dxy_IIIE(a, c, d, b);
		return incircle2dxy_IIIE(b, d, c, a);
	}

	return incirclexy_indirect_IIII(a, b, c, d);
}


// These functions assume that point is an SSI
inline bool genericPoint::getIntervalLambda(interval_number& lx, interval_number& ly, interval_number& d) const {
	return toSSI().getIntervalLambda(lx, ly, d);
}

inline void genericPoint::getExactLambda(double** lx, int& lxl, double** ly, int& lyl, double** d, int& dl) const {
	toSSI().getExactLambda(lx, lxl, ly, lyl, d, dl);
}

inline void genericPoint::getBigfloatLambda(bigfloat& lx, bigfloat& ly, bigfloat& d) const {
	toSSI().getBigfloatLambda(lx, ly, d);
}

// These functions assume that point is an implicit 3D
inline bool genericPoint::getIntervalLambda(interval_number& lx, interval_number& ly, interval_number& lz, interval_number& d) const {
	if (isLPI()) return toLPI().getIntervalLambda(lx, ly, lz, d);
	else if (isTPI()) return toTPI().getIntervalLambda(lx, ly, lz, d);
	else return toLNC().getIntervalLambda(lx, ly, lz, d);
}

inline void genericPoint::getExactLambda(double** lx, int& lxl, double** ly, int& lyl, double** lz, int& lzl, double** d, int& dl) const {
	if (isLPI()) toLPI().getExactLambda(lx, lxl, ly, lyl, lz, lzl, d, dl);
	else if (isTPI()) toTPI().getExactLambda(lx, lxl, ly, lyl, lz, lzl, d, dl);
	else toLNC().getExactLambda(lx, lxl, ly, lyl, lz, lzl, d, dl);
}

inline void genericPoint::getBigfloatLambda(bigfloat& lx, bigfloat& ly, bigfloat& lz, bigfloat& d) const {
	if (isLPI()) toLPI().getBigfloatLambda(lx, ly, lz, d);
	else if (isTPI()) toTPI().getBigfloatLambda(lx, ly, lz, d);
	else toLNC().getBigfloatLambda(lx, ly, lz, d);
}

// Type-specific lambdas

inline bool implicitPoint2D_SSI::getIntervalLambda(interval_number& lx, interval_number& ly, interval_number &d) const
{
	if (needsIntervalLambda())
	{
		lambda2d_SSI_interval(l1_1.X(), l1_1.Y(), l1_2.X(), l1_2.Y(), l2_1.X(), l2_1.Y(), l2_2.X(), l2_2.Y(), dfilter_lambda_x, dfilter_lambda_y, dfilter_denominator);
		if (dfilter_denominator.isNegative()) {
			dfilter_lambda_x.invert();
			dfilter_lambda_y.invert();
			dfilter_denominator.invert();
		}
	}

	lx = dfilter_lambda_x;
	ly = dfilter_lambda_y;
	d = dfilter_denominator;
	return (dfilter_denominator.signIsReliable());
}

inline bool implicitPoint3D_LPI::getIntervalLambda(interval_number& lx, interval_number& ly, interval_number& lz, interval_number &d) const
{
	if (needsIntervalLambda())
	{
		lambda3d_LPI_interval(P().X(), P().Y(), P().Z(), Q().X(), Q().Y(), Q().Z(), R().X(), R().Y(), R().Z(), S().X(), S().Y(), S().Z(), T().X(), T().Y(), T().Z(), dfilter_lambda_x, dfilter_lambda_y, dfilter_lambda_z, dfilter_denominator);
		if (dfilter_denominator.isNegative()) {
			dfilter_lambda_x.invert();
			dfilter_lambda_y.invert();
			dfilter_lambda_z.invert();
			dfilter_denominator.invert();
		}
	}
	lx = dfilter_lambda_x;
	ly = dfilter_lambda_y;
	lz = dfilter_lambda_z;
	d = dfilter_denominator;
	return (dfilter_denominator.signIsReliable());
}

inline bool implicitPoint3D_TPI::getIntervalLambda(interval_number& lx, interval_number& ly, interval_number& lz, interval_number& d) const
{
	if (needsIntervalLambda())
	{
	 lambda3d_TPI_interval(
		V1().X(), V1().Y(), V1().Z(), V2().X(), V2().Y(), V2().Z(), V3().X(), V3().Y(), V3().Z(),
		W1().X(), W1().Y(), W1().Z(), W2().X(), W2().Y(), W2().Z(), W3().X(), W3().Y(), W3().Z(),
		U1().X(), U1().Y(), U1().Z(), U2().X(), U2().Y(), U2().Z(), U3().X(), U3().Y(), U3().Z(),
		dfilter_lambda_x, dfilter_lambda_y, dfilter_lambda_z, dfilter_denominator);
	 if (dfilter_denominator.isNegative()) {
		dfilter_lambda_x.invert();
		dfilter_lambda_y.invert();
		dfilter_lambda_z.invert();
		dfilter_denominator.invert();
	 }
	}

	lx = dfilter_lambda_x;
	ly = dfilter_lambda_y;
	lz = dfilter_lambda_z;
	d = dfilter_denominator;
	return (dfilter_denominator.signIsReliable());
}

inline bool implicitPoint3D_LNC::getIntervalLambda(interval_number& lx, interval_number& ly, interval_number& lz, interval_number& d) const
{
	if (needsIntervalLambda())
	{
		lambda3d_LNC_interval(P().X(), P().Y(), P().Z(), Q().X(), Q().Y(), Q().Z(), T(), dfilter_lambda_x, dfilter_lambda_y, dfilter_lambda_z, dfilter_denominator);
	}
	lx = dfilter_lambda_x;
	ly = dfilter_lambda_y;
	lz = dfilter_lambda_z;
	d = dfilter_denominator;
	return true;
}


inline void implicitPoint2D_SSI::getExactLambda(double **lx, int& lxl, double **ly, int& lyl, double **d, int& dl) const
{
	lambda2d_SSI_exact(l1_1.X(), l1_1.Y(), l1_2.X(), l1_2.Y(), l2_1.X(), l2_1.Y(), l2_2.X(), l2_2.Y(), lx, lxl, ly, lyl, d, dl);
	if ((*d)[dl - 1] < 0)
	{
		expansionObject::Gen_Invert(lxl, *lx);
		expansionObject::Gen_Invert(lyl, *ly);
		expansionObject::Gen_Invert(dl, *d);
	}
}

// Keeps lambda/d pairs as close to one as possible to avoid under/overflows
inline void normalizeLambda3D(double* lx, int& lxl, double* ly, int& lyl, double* lz, int& lzl, double* d, int& dl)
{
	double maxd, maxsd, ad, aad;
	maxsd = expansionObject::To_Double(lxl, lx);
	maxd = fabs(maxsd);
	if ((aad = fabs((ad = expansionObject::To_Double(lyl, ly)))) > maxd) { maxd = aad; maxsd = ad; }
	if ((aad = fabs((ad = expansionObject::To_Double(lzl, lz)))) > maxd) { maxd = aad; maxsd = ad; }
	if ((aad = fabs((ad = expansionObject::To_Double(dl, d)))) > maxd) { maxd = aad; maxsd = ad; }

	int e;
	frexp(maxsd, &e);
	const double m = ldexp(2, -e);

	expansionObject::ExactScale(lxl, lx, m);
	expansionObject::ExactScale(lyl, ly, m);
	expansionObject::ExactScale(lzl, lz, m);
	expansionObject::ExactScale(dl, d, m);
}

inline void implicitPoint3D_LPI::getExactLambda(double **lx, int& lxl, double **ly, int& lyl, double **lz, int& lzl, double **d, int& dl) const
{
	lambda3d_LPI_exact(P().X(), P().Y(), P().Z(), Q().X(), Q().Y(), Q().Z(), R().X(), R().Y(), R().Z(), S().X(), S().Y(), S().Z(), T().X(), T().Y(), T().Z(), lx, lxl, ly, lyl, lz, lzl, d, dl);
	if ((*d)[dl - 1] < 0)
	{
		expansionObject::Gen_Invert(lxl, *lx);
		expansionObject::Gen_Invert(lyl, *ly);
		expansionObject::Gen_Invert(lzl, *lz);
		expansionObject::Gen_Invert(dl, *d);
	}
	normalizeLambda3D(*lx, lxl, *ly, lyl, *lz, lzl, *d, dl);
}

inline void implicitPoint3D_TPI::getExactLambda(double **lx, int& lxl, double **ly, int& lyl, double **lz, int& lzl, double **d, int& dl) const
{
	lambda3d_TPI_exact(V1().X(), V1().Y(), V1().Z(), V2().X(), V2().Y(), V2().Z(), V3().X(), V3().Y(), V3().Z(),
		W1().X(), W1().Y(), W1().Z(), W2().X(), W2().Y(), W2().Z(), W3().X(), W3().Y(), W3().Z(),
		U1().X(), U1().Y(), U1().Z(), U2().X(), U2().Y(), U2().Z(), U3().X(), U3().Y(), U3().Z(), lx, lxl, ly, lyl, lz, lzl, d, dl);
	if ((*d)[dl - 1] < 0)
	{
		expansionObject::Gen_Invert(lxl, *lx);
		expansionObject::Gen_Invert(lyl, *ly);
		expansionObject::Gen_Invert(lzl, *lz);
		expansionObject::Gen_Invert(dl, *d);
	}
	normalizeLambda3D(*lx, lxl, *ly, lyl, *lz, lzl, *d, dl);
}

inline void implicitPoint3D_LNC::getExactLambda(double** lx, int& lxl, double** ly, int& lyl, double** lz, int& lzl, double** d, int& dl) const
{
	lambda3d_LNC_exact(P().X(), P().Y(), P().Z(), Q().X(), Q().Y(), Q().Z(), T(), lx, lxl, ly, lyl, lz, lzl, d, dl);
	normalizeLambda3D(*lx, lxl, *ly, lyl, *lz, lzl, *d, dl);
}

inline void implicitPoint2D_SSI::getBigfloatLambda(bigfloat& lx, bigfloat& ly, bigfloat& d) const
{
	lambda2d_SSI_bigfloat(l1_1.X(), l1_1.Y(), l1_2.X(), l1_2.Y(), l2_1.X(), l2_1.Y(), l2_2.X(), l2_2.Y(), lx, ly, d);
	if (sgn(d) < 0)
	{
		lx = -lx;
		ly = -ly;
		d = -d;
	}
}

inline void implicitPoint3D_LPI::getBigfloatLambda(bigfloat& lx, bigfloat& ly, bigfloat& lz, bigfloat& d) const
{
	lambda3d_LPI_bigfloat(P().X(), P().Y(), P().Z(), Q().X(), Q().Y(), Q().Z(), R().X(), R().Y(), R().Z(), S().X(), S().Y(), S().Z(), T().X(), T().Y(), T().Z(), lx, ly, lz, d);
	if (sgn(d) < 0)
	{
		lx = -lx;
		ly = -ly;
		lz = -lz;
		d = -d;
	}
}

inline void implicitPoint3D_TPI::getBigfloatLambda(bigfloat& lx, bigfloat& ly, bigfloat& lz, bigfloat& d) const
{
	lambda3d_TPI_bigfloat(V1().X(), V1().Y(), V1().Z(), V2().X(), V2().Y(), V2().Z(), V3().X(), V3().Y(), V3().Z(),
		W1().X(), W1().Y(), W1().Z(), W2().X(), W2().Y(), W2().Z(), W3().X(), W3().Y(), W3().Z(),
		U1().X(), U1().Y(), U1().Z(), U2().X(), U2().Y(), U2().Z(), U3().X(), U3().Y(), U3().Z(), lx, ly, lz, d);
	if (sgn(d) < 0)
	{
		lx = -lx;
		ly = -ly;
		lz = -lz;
		d = -d;
	}
}

inline void implicitPoint3D_LNC::getBigfloatLambda(bigfloat& lx, bigfloat& ly, bigfloat& lz, bigfloat& d) const
{
	lambda3d_LNC_bigfloat(P().X(), P().Y(), P().Z(), Q().X(), Q().Y(), Q().Z(), T(), lx, ly, lz, d);
}


inline bool genericPoint::apapExplicit(explicitPoint2D& e) const
{
	if (isExplicit2D()) e = toExplicit2D();
	else {
		double l1x_p[128], * l1x = l1x_p, l1y_p[128], * l1y = l1y_p, d1_p[128], * d1 = d1_p;
		int l1x_len, l1y_len, d1_len;
		getExactLambda(&l1x, l1x_len, &l1y, l1y_len, &d1, d1_len);
		const double lambda_x = expansionObject::To_Double(l1x_len, l1x);
		const double lambda_y = expansionObject::To_Double(l1y_len, l1y);
		const double lambda_d = expansionObject::To_Double(d1_len, d1);
		if (l1x_p != l1x) free(l1x);
		if (l1y_p != l1y) free(l1y);
		if (d1_p != d1) free(d1);
		if (lambda_d == 0) return false;

		e = explicitPoint2D(lambda_x / lambda_d, lambda_y / lambda_d);
	}
	return true;
}

inline bool genericPoint::approxExplicit(explicitPoint2D& e) const
{
	if (isExplicit2D()) e = toExplicit2D();
	else {
		double lambda_x, lambda_y, lambda_d;
		interval_number ilx, ily, id;
		if (!getIntervalLambda(ilx, ily, id)) return apapExplicit(e);
		else
		{
			lambda_x = ilx.sup() + ilx.inf();
			lambda_y = ily.sup() + ily.inf();
			lambda_d = id.sup() + id.inf();
		}
		e = explicitPoint2D(lambda_x / lambda_d, lambda_y / lambda_d);
	}
	return true;
}

inline bool genericPoint::apapExplicit(explicitPoint3D& e) const
{
	if (isExplicit3D()) e = toExplicit3D();
	else {
		double l1z_p[128], * l1z = l1z_p, l1x_p[128], * l1x = l1x_p, l1y_p[128], * l1y = l1y_p, d1_p[128], * d1 = d1_p;
		int l1z_len, l1x_len, l1y_len, d1_len;
		getExactLambda(&l1x, l1x_len, &l1y, l1y_len, &l1z, l1z_len, &d1, d1_len);
		const double lambda_x = expansionObject::To_Double(l1x_len, l1x);
		const double lambda_y = expansionObject::To_Double(l1y_len, l1y);
		const double lambda_z = expansionObject::To_Double(l1z_len, l1z);
		const double lambda_d = expansionObject::To_Double(d1_len, d1);
		if (l1z_p != l1z) free(l1z);
		if (l1x_p != l1x) free(l1x);
		if (l1y_p != l1y) free(l1y);
		if (d1_p != d1) free(d1);
		if (lambda_d == 0) return false;
		e = explicitPoint3D(lambda_x / lambda_d, lambda_y / lambda_d, lambda_z / lambda_d);
	}
	return true;
}


inline bool genericPoint::approxExplicit(explicitPoint3D& e) const
{
	if (isExplicit3D()) e = toExplicit3D();
	else {
		double lambda_x, lambda_y, lambda_z, lambda_d;
		interval_number ilx, ily, ilz, id;
		if (!getIntervalLambda(ilx, ily, ilz, id)) return apapExplicit(e);
		else
		{
			lambda_x = ilx.sup() + ilx.inf();
			lambda_y = ily.sup() + ily.inf();
			lambda_z = ilz.sup() + ilz.inf();
			lambda_d = id.sup() + id.inf();
		}
		e = explicitPoint3D(lambda_x / lambda_d, lambda_y / lambda_d, lambda_z / lambda_d);
	}
	return true;
}


inline bool genericPoint::getApproxXYCoordinates(double& x, double& y, bool apap) const
{
	if (is2D())
	{
		explicitPoint2D op;
		if (apap && !apapExplicit(op)) return false;
		if (!apap && !approxExplicit(op)) return false;
		x = op.X(); y = op.Y();
		return true;
	}
	if (is3D())
	{
		explicitPoint3D op;
		if (apap && !apapExplicit(op)) return false;
		if (!apap && !approxExplicit(op)) return false;
		x = op.X(); y = op.Y();
		return true;
	}
	ip_error("genericPoint::getApproxXYCoordinates - should not happen\n");
	return false;
}

inline bool genericPoint::getApproxXYZCoordinates(double& x, double& y, double& z, bool apap) const
{
	if (is3D())
	{
		explicitPoint3D op;
		if (apap && !apapExplicit(op)) return false;
		if (!apap && !approxExplicit(op)) return false;
		x = op.X(); y = op.Y(); z = op.Z();
		return true;
	}
	ip_error("genericPoint::getApproxXYZCoordinates - should not happen\n");
	return false;
}

inline bool genericPoint::getExactXYCoordinates(bigrational& x, bigrational& y) const
{
	if (isExplicit2D()) return toExplicit2D().getExactXYCoordinates(x, y);
	else if (isSSI()) return toSSI().getExactXYCoordinates(x, y);
	else ip_error("genericPoint::getExactXYCoordinates - should not happen\n");
	return false;
}

inline bool genericPoint::getExactXYZCoordinates(bigrational& x, bigrational& y, bigrational& z) const
{
	if (isExplicit3D()) return toExplicit3D().getExactXYZCoordinates(x, y, z);
	else if (isLPI()) return toLPI().getExactXYZCoordinates(x, y, z);
	else if (isTPI()) return toTPI().getExactXYZCoordinates(x, y, z);
	else if (isLNC()) return toLNC().getExactXYZCoordinates(x, y, z);
	else if (isExplicit2D()) { z = bigfloat(0); return toExplicit2D().getExactXYCoordinates(x, y); }
	else if (isSSI()) { z = bigfloat(0); return toSSI().getExactXYCoordinates(x, y); }
	else ip_error("genericPoint::getExactXYZCoordinates - should not happen\n");
	return false;
}

inline bool implicitPoint2D_SSI::getExactXYCoordinates(bigrational& x, bigrational& y) const
{
	bigfloat lx, ly, d;
	getBigfloatLambda(lx, ly, d);
	if (sgn(d) == 0) return false;
	const bigrational rd(d);
	x = bigrational(lx) / rd;
	y = bigrational(ly) / rd;
	return true;
}

inline bool implicitPoint3D_LPI::getExactXYZCoordinates(bigrational& x, bigrational& y, bigrational& z) const
{
	bigfloat lx, ly, lz, d;
	getBigfloatLambda(lx, ly, lz, d);
	if (sgn(d) == 0) return false;
	const bigrational rd(d);
	x = bigrational(lx) / rd;
	y = bigrational(ly) / rd;
	z = bigrational(lz) / rd;
	return true;
}

inline bool implicitPoint3D_TPI::getExactXYZCoordinates(bigrational& x, bigrational& y, bigrational& z) const
{
	bigfloat lx, ly, lz, d;
	getBigfloatLambda(lx, ly, lz, d);
	if (sgn(d) == 0) return false;
	const bigrational rd(d);
	x = bigrational(lx) / rd;
	y = bigrational(ly) / rd;
	z = bigrational(lz) / rd;
	return true;
}

inline bool implicitPoint3D_LNC::getExactXYZCoordinates(bigrational& x, bigrational& y, bigrational& z) const
{
	bigfloat lx, ly, lz, d;
	getBigfloatLambda(lx, ly, lz, d);
	x = bigrational(lx);
	y = bigrational(ly);
	z = bigrational(lz);
	return true;
}

inline ostream& operator<<(ostream& os, const genericPoint& p)
{
	if (p.isExplicit2D()) return os << p.toExplicit2D();
	else if (p.isExplicit3D()) return os << p.toExplicit3D();
	else if (p.isSSI()) return os << p.toSSI();
	else if (p.isLPI()) return os << p.toLPI();
	else if (p.isTPI()) return os << p.toTPI();
	else if (p.isLNC()) return os << p.toLNC();
	else ip_error("genericPoint::operator<< - should not happen\n");
	return os;
}

inline int maxComponentInTriangleNormal_filtered(double ov1x, double ov1y, double ov1z, double ov2x, double ov2y, double ov2z, double ov3x, double ov3y, double ov3z)
{
	double v3x = ov3x - ov2x;
	double v3y = ov3y - ov2y;
	double v3z = ov3z - ov2z;
	double v2x = ov2x - ov1x;
	double v2y = ov2y - ov1y;
	double v2z = ov2z - ov1z;
	double nvx1 = v2y * v3z;
	double nvx2 = v2z * v3y;
	double nvx = nvx1 - nvx2;
	double nvy1 = v3x * v2z;
	double nvy2 = v3z * v2x;
	double nvy = nvy1 - nvy2;
	double nvz1 = v2x * v3y;
	double nvz2 = v2y * v3x;
	double nvz = nvz1 - nvz2;

	double _tmp_fabs, max_var = 0;
	if ((_tmp_fabs = fabs(v3x)) > max_var) max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(v3y)) > max_var) max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(v3z)) > max_var) max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(v2x)) > max_var) max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(v2y)) > max_var) max_var = _tmp_fabs;
	if ((_tmp_fabs = fabs(v2z)) > max_var) max_var = _tmp_fabs;
	double epsilon = 8.88395e-016 * max_var * max_var;

	double nvxc = fabs(nvx);
	double nvyc = fabs(nvy);
	double nvzc = fabs(nvz);
	double nv = nvxc;
	if (nvyc > nv) nv = nvyc;
	if (nvzc > nv) nv = nvzc;

	if (nv > epsilon)
	{
		if (nv == nvxc) return 0;
		if (nv == nvyc) return 1;
		if (nv == nvzc) return 2;
	}
	return -1;
}

inline int maxComponentInTriangleNormal_exact(double ov1x, double ov1y, double ov1z, double ov2x, double ov2y, double ov2z, double ov3x, double ov3y, double ov3z)
{
	expansionObject o;
	double v3x[2];
	o.two_Diff(ov3x, ov2x, v3x);
	double v3y[2];
	o.two_Diff(ov3y, ov2y, v3y);
	double v3z[2];
	o.two_Diff(ov3z, ov2z, v3z);
	double v2x[2];
	o.two_Diff(ov2x, ov1x, v2x);
	double v2y[2];
	o.two_Diff(ov2y, ov1y, v2y);
	double v2z[2];
	o.two_Diff(ov2z, ov1z, v2z);
	double nvx1[8];
	o.Two_Two_Prod(v2y, v3z, nvx1);
	double nvx2[8];
	o.Two_Two_Prod(v2z, v3y, nvx2);
	double nvx[16];
	int nvx_len = o.Gen_Diff(8, nvx1, 8, nvx2, nvx);
	double nvy1[8];
	o.Two_Two_Prod(v3x, v2z, nvy1);
	double nvy2[8];
	o.Two_Two_Prod(v3z, v2x, nvy2);
	double nvy[16];
	int nvy_len = o.Gen_Diff(8, nvy1, 8, nvy2, nvy);
	double nvz1[8];
	o.Two_Two_Prod(v2x, v3y, nvz1);
	double nvz2[8];
	o.Two_Two_Prod(v2y, v3x, nvz2);
	double nvz[16];
	int nvz_len = o.Gen_Diff(8, nvz1, 8, nvz2, nvz);

	double nvxc = fabs(nvx[nvx_len - 1]);
	double nvyc = fabs(nvy[nvy_len - 1]);
	double nvzc = fabs(nvz[nvz_len - 1]);
	double nv = nvxc;
	if (nvyc > nv) nv = nvyc;
	if (nvzc > nv) return 2;
	if (nv == nvxc) return 0;
	return 1;
}

inline int genericPoint::maxComponentInTriangleNormal(double ov1x, double ov1y, double ov1z, double ov2x, double ov2y, double ov2z, double ov3x, double ov3y, double ov3z)
{
	int ret;
	if ((ret = maxComponentInTriangleNormal_filtered(ov1x, ov1y, ov1z, ov2x, ov2y, ov2z, ov3x, ov3y, ov3z)) >= 0) return ret;
	return maxComponentInTriangleNormal_exact(ov1x, ov1y, ov1z, ov2x, ov2y, ov2z, ov3x, ov3y, ov3z);
}


/////////////////////////////////////////////
//
// Derived predicates
//
/////////////////////////////////////////////

inline bool genericPoint::innerSegmentsCross(const genericPoint& A, const genericPoint& B, const genericPoint& P, const genericPoint& Q)
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

inline bool genericPoint::segmentsCross(const genericPoint& A, const genericPoint& B, const genericPoint& P, const genericPoint& Q)
{
	int o11, o12, o21, o22;

	o11 = orient2Dxy(P, A, B);
	o12 = orient2Dxy(Q, B, A);
	o21 = orient2Dxy(A, P, Q);
	o22 = orient2Dxy(B, Q, P);
	if ((o11 || o12) && (o11 * o12 >= 0) && (o21 || o22) && (o21 * o22 >= 0)) return true;

	o11 = orient2Dyz(P, A, B);
	o12 = orient2Dyz(Q, B, A);
	o21 = orient2Dyz(A, P, Q);
	o22 = orient2Dyz(B, Q, P);
	if ((o11 || o12) && (o11 * o12 >= 0) && (o21 || o22) && (o21 * o22 >= 0)) return true;

	o11 = orient2Dzx(P, A, B);
	o12 = orient2Dzx(Q, B, A);
	o21 = orient2Dzx(A, P, Q);
	o22 = orient2Dzx(B, Q, P);
	if ((o11 || o12) && (o11 * o12 >= 0) && (o21 || o22) && (o21 * o22 >= 0)) return true;

	return false;
}

inline bool genericPoint::innerSegmentCrossesInnerTriangle(const genericPoint& s1, const genericPoint& s2, const genericPoint& v1, const genericPoint& v2, const genericPoint& v3)
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

inline bool genericPoint::pointInInnerSegment(const genericPoint& p, const genericPoint& v1, const genericPoint& v2)
{
	if (misaligned(p, v1, v2)) return false;

	int lt2, lt3;
	lt2 = lessThanOnX(v1, p);
	lt3 = lessThanOnX(p, v2);
	if (lt2) return (lt2 == lt3);
	lt2 = lessThanOnY(v1, p);
	lt3 = lessThanOnY(p, v2);
	if (lt2) return (lt2 == lt3);
	lt2 = lessThanOnZ(v1, p);
	lt3 = lessThanOnZ(p, v2);
	if (lt2) return (lt2 == lt3);
	return false;
}

inline bool genericPoint::pointInSegment(const genericPoint& p, const genericPoint& v1, const genericPoint& v2)
{
	if (misaligned(p, v1, v2)) return false;

	int lt2x = lessThanOnX(v1, p);
	int lt3x = lessThanOnX(p, v2);
	if (lt2x && lt3x) return (lt2x == lt3x);
	int lt2y = lessThanOnY(v1, p);
	int lt3y = lessThanOnY(p, v2);
	if (lt2y && lt3y) return (lt2y == lt3y);
	int lt2z = lessThanOnZ(v1, p);
	int lt3z = lessThanOnZ(p, v2);
	if (lt2z && lt3z) return (lt2z == lt3z);

	return ((lt2x == 0 && lt2y == 0 && lt2z == 0) || (lt3x == 0 && lt3y == 0 && lt3z == 0));
}

inline bool genericPoint::pointInTriangle(const genericPoint& P, const genericPoint& A, const genericPoint& B, const genericPoint& C)
{
	int o1, o2, o3;
	o1 = orient2Dxy(P, A, B);
	o2 = orient2Dxy(P, B, C);
	o3 = orient2Dxy(P, C, A);
	if (o1 || o2 || o3) return ((o1 >= 0 && o2 >= 0 && o3 >= 0) || (o1 <= 0 && o2 <= 0 && o3 <= 0));
	o1 = orient2Dyz(P, A, B);
	o2 = orient2Dyz(P, B, C);
	o3 = orient2Dyz(P, C, A);
	if (o1 || o2 || o3) return ((o1 >= 0 && o2 >= 0 && o3 >= 0) || (o1 <= 0 && o2 <= 0 && o3 <= 0));
	o1 = orient2Dzx(P, A, B);
	o2 = orient2Dzx(P, B, C);
	o3 = orient2Dzx(P, C, A);
	return ((o1 >= 0 && o2 >= 0 && o3 >= 0) || (o1 <= 0 && o2 <= 0 && o3 <= 0));
}


inline bool genericPoint::pointInTriangle(const genericPoint& P, const genericPoint& A, const genericPoint& B, const genericPoint& C, int& o1, int& o2, int& o3)
{
	o1 = orient2Dxy(P, A, B);
	o2 = orient2Dxy(P, B, C);
	o3 = orient2Dxy(P, C, A);
	if (o1 || o2 || o3) return ((o1 >= 0 && o2 >= 0 && o3 >= 0) || (o1 <= 0 && o2 <= 0 && o3 <= 0));
	o1 = orient2Dyz(P, A, B);
	o2 = orient2Dyz(P, B, C);
	o3 = orient2Dyz(P, C, A);
	if (o1 || o2 || o3) return ((o1 >= 0 && o2 >= 0 && o3 >= 0) || (o1 <= 0 && o2 <= 0 && o3 <= 0));
	o1 = orient2Dzx(P, A, B);
	o2 = orient2Dzx(P, B, C);
	o3 = orient2Dzx(P, C, A);
	return ((o1 >= 0 && o2 >= 0 && o3 >= 0) || (o1 <= 0 && o2 <= 0 && o3 <= 0));
}

inline bool genericPoint::pointInInnerTriangle(const genericPoint& P, const genericPoint& A, const genericPoint& B, const genericPoint& C)
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

inline bool genericPoint::lineCrossesInnerTriangle(const genericPoint& s1, const genericPoint& s2, const genericPoint& v1, const genericPoint& v2, const genericPoint& v3)
{
	const int o1 = genericPoint::orient3D(s1, s2, v1, v2);
	const int o2 = genericPoint::orient3D(s1, s2, v2, v3);
	if ((o1 >= 0 && o2 <= 0) || (o1 <= 0 && o2 >= 0)) return false;
	const int o3 = genericPoint::orient3D(s1, s2, v3, v1);
	if ((o1 >= 0 && o3 <= 0) || (o1 <= 0 && o3 >= 0)) return false;
	if ((o2 >= 0 && o3 <= 0) || (o2 <= 0 && o3 >= 0)) return false;
	return true;
}

inline bool genericPoint::lineCrossesTriangle(const genericPoint& s1, const genericPoint& s2, const genericPoint& v1, const genericPoint& v2, const genericPoint& v3)
{
	const int o1 = genericPoint::orient3D(s1, s2, v1, v2);
	const int o2 = genericPoint::orient3D(s1, s2, v2, v3);
	if ((o1 > 0 && o2 < 0) || (o1 < 0 && o2 > 0)) return false;
	const int o3 = genericPoint::orient3D(s1, s2, v3, v1);
	if ((o1 > 0 && o3 < 0) || (o1 < 0 && o3 > 0)) return false;
	if ((o2 > 0 && o3 < 0) || (o2 < 0 && o3 > 0)) return false;
	return true;
}

inline bool genericPoint::innerSegmentCrossesTriangle(const genericPoint& s1, const genericPoint& s2, const genericPoint& v1, const genericPoint& v2, const genericPoint& v3)
{
	const int o1 = genericPoint::orient3D(s1, v1, v2, v3); if (o1 == 0) return false;
	const int o2 = genericPoint::orient3D(s2, v1, v2, v3); if (o2 == 0) return false;

	if ((o1 > 0 && o2 > 0) || (o1 < 0 && o2 < 0)) return false;
	return lineCrossesTriangle(s1, s2, v1, v2, v3);
}



inline bool genericPoint::pointInInnerSegment(const genericPoint& p, const genericPoint& v1, const genericPoint& v2, int xyz)
{
	int lt2, lt3;
	if (xyz == 0)
	{
		if (orient2Dyz(p, v1, v2)) return false;
		lt2 = lessThanOnY(v1, p);
		lt3 = lessThanOnY(p, v2);
		if (lt2) return (lt2 == lt3);
		lt2 = lessThanOnZ(v1, p);
		lt3 = lessThanOnZ(p, v2);
	}
	else if (xyz == 1)
	{
		if (orient2Dzx(p, v1, v2)) return false;
		lt2 = lessThanOnX(v1, p);
		lt3 = lessThanOnX(p, v2);
		if (lt2) return (lt2 == lt3);
		lt2 = lessThanOnZ(v1, p);
		lt3 = lessThanOnZ(p, v2);
	}
	else
	{
		if (orient2Dxy(p, v1, v2)) return false;
		lt2 = lessThanOnX(v1, p);
		lt3 = lessThanOnX(p, v2);
		if (lt2) return (lt2 == lt3);
		lt2 = lessThanOnY(v1, p);
		lt3 = lessThanOnY(p, v2);
	}
	return (lt2 && (lt2 == lt3));
}

inline bool genericPoint::pointInSegment(const genericPoint& p, const genericPoint& v1, const genericPoint& v2, int xyz)
{
	int lt2, lt3, lt4, lt5;
	if (xyz == 0)
	{
		if (orient2Dyz(p, v1, v2)) return false;
		lt2 = lessThanOnY(v1, p);
		lt3 = lessThanOnY(p, v2);
		if (lt2 && lt3) return (lt2 == lt3);
		lt4 = lessThanOnZ(v1, p);
		lt5 = lessThanOnZ(p, v2);
	}
	else if (xyz == 1)
	{
		if (orient2Dzx(p, v1, v2)) return false;
		lt2 = lessThanOnX(v1, p);
		lt3 = lessThanOnX(p, v2);
		if (lt2 && lt3) return (lt2 == lt3);
		lt4 = lessThanOnZ(v1, p);
		lt5 = lessThanOnZ(p, v2);
	}
	else
	{
		if (orient2Dxy(p, v1, v2)) return false;
		lt2 = lessThanOnX(v1, p);
		lt3 = lessThanOnX(p, v2);
		if (lt2 && lt3) return (lt2 == lt3);
		lt4 = lessThanOnY(v1, p);
		lt5 = lessThanOnY(p, v2);
	}
	return ((lt2 == 0 && lt4 == 0) || (lt3 == 0 && lt5 == 0));
}

inline bool genericPoint::pointInInnerTriangle(const genericPoint& P, const genericPoint& A, const genericPoint& B, const genericPoint& C, int xyz)
{
	int o1, o2, o3;
	if (xyz == 2)
	{
		o1 = genericPoint::orient2Dxy(P, A, B);
		o2 = genericPoint::orient2Dxy(P, B, C);
		o3 = genericPoint::orient2Dxy(P, C, A);
	}
	else if (xyz == 0)
	{
		o1 = genericPoint::orient2Dyz(P, A, B);
		o2 = genericPoint::orient2Dyz(P, B, C);
		o3 = genericPoint::orient2Dyz(P, C, A);
	}
	else
	{
		o1 = genericPoint::orient2Dzx(P, A, B);
		o2 = genericPoint::orient2Dzx(P, B, C);
		o3 = genericPoint::orient2Dzx(P, C, A);
	}
	return ((o1 > 0 && o2 > 0 && o3 > 0) || (o1 < 0 && o2 < 0 && o3 < 0));
}

inline bool genericPoint::pointInTriangle(const genericPoint& P, const genericPoint& A, const genericPoint& B, const genericPoint& C, int xyz)
{
	int o1, o2, o3;
	if (xyz == 2)
	{
		o1 = genericPoint::orient2Dxy(P, A, B);
		o2 = genericPoint::orient2Dxy(P, B, C);
		o3 = genericPoint::orient2Dxy(P, C, A);
	}
	else if (xyz == 0)
	{
		o1 = genericPoint::orient2Dyz(P, A, B);
		o2 = genericPoint::orient2Dyz(P, B, C);
		o3 = genericPoint::orient2Dyz(P, C, A);
	}
	else
	{
		o1 = genericPoint::orient2Dzx(P, A, B);
		o2 = genericPoint::orient2Dzx(P, B, C);
		o3 = genericPoint::orient2Dzx(P, C, A);
	}
	return ((o1 >= 0 && o2 >= 0 && o3 >= 0) || (o1 <= 0 && o2 <= 0 && o3 <= 0));
}

inline bool genericPoint::innerSegmentsCross(const genericPoint& A, const genericPoint& B, const genericPoint& P, const genericPoint& Q, int xyz)
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

inline bool genericPoint::segmentsCross(const genericPoint& A, const genericPoint& B, const genericPoint& P, const genericPoint& Q, int xyz)
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

	return ((o11 || o12) && (o11 * o12 >= 0) && (o21 || o22) && (o21 * o22 >= 0));
}
