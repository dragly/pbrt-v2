
/*
    pbrt source code Copyright(c) 1998-2009 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef PBRT_CORE_TEXTURE_H
#define PBRT_CORE_TEXTURE_H

// core/texture.h*
#include "pbrt.h"
#include "spectrum.h"
#include "geometry.h"
#include "transform.h"
#include "memory.h"

// Texture Declarations
class TextureMapping2D {
public:
    // TextureMapping2D Interface
    virtual ~TextureMapping2D() { }
    virtual void Map(const DifferentialGeometry &dg,
                     float *s, float *t, float *dsdx, float *dtdx,
                     float *dsdy, float *dtdy) const = 0;
};


class UVMapping2D : public TextureMapping2D {
public:
    // UVMapping2D Public Methods
    UVMapping2D(float su = 1, float sv = 1,
        float du = 0, float dv = 0);
    void Map(const DifferentialGeometry &dg, float *s, float *t,
        float *dsdx, float *dtdx,
        float *dsdy, float *dtdy) const;
private:
    float su, sv, du, dv;
};


class SphericalMapping2D : public TextureMapping2D {
public:
    // SphericalMapping2D Public Methods
    SphericalMapping2D(const Transform &toSph)
        : WorldToTexture(toSph) {
    }
    void Map(const DifferentialGeometry &dg, float *s, float *t,
        float *dsdx, float *dtdx,
        float *dsdy, float *dtdy) const;
private:
    void sphere(const Point &P, float *s, float *t) const;
    Transform WorldToTexture;
};


class CylindricalMapping2D : public TextureMapping2D {
public:
    // CylindricalMapping2D Public Methods
    CylindricalMapping2D(const Transform &toCyl)
        : WorldToTexture(toCyl) {
    }
    void Map(const DifferentialGeometry &dg, float *s, float *t,
        float *dsdx, float *dtdx,
        float *dsdy, float *dtdy) const;
private:
    void cylinder(const Point &P, float *s, float *t) const;
    Transform WorldToTexture;
};


class PlanarMapping2D : public TextureMapping2D {
public:
    // PlanarMapping2D Public Methods
    PlanarMapping2D(const Vector &v1, const Vector &v2,
        float du = 0, float dv = 0);
    void Map(const DifferentialGeometry &dg, float *s, float *t,
        float *dsdx, float *dtdx,
        float *dsdy, float *dtdy) const;
private:
    Vector vs, vt;
    float ds, dt;
};


class TextureMapping3D {
public:
    // TextureMapping3D Interface
    virtual ~TextureMapping3D() { }
    virtual Point Map(const DifferentialGeometry &dg,
                      Vector *dpdx, Vector *dpdy) const = 0;
};


class IdentityMapping3D : public TextureMapping3D {
public:
    IdentityMapping3D(const Transform &x)
        : WorldToTexture(x) { }
    Point Map(const DifferentialGeometry &dg, Vector *dpdx,
              Vector *dpdy) const;
private:
    Transform WorldToTexture;
};


template <typename T> class Texture : public ReferenceCounted {
public:
    // Texture Interface
    virtual T Evaluate(const DifferentialGeometry &) const = 0;
    virtual ~Texture() { }
};


float Lanczos(float, float tau=2);
float Noise(float x, float y = .5f, float z = .5f);
float Noise(const Point &P);
float FBm(const Point &P, const Vector &dpdx, const Vector &dpdy,
    float omega, int octaves);
float Turbulence(const Point &P, const Vector &dpdx, const Vector &dpdy,
    float omega, int octaves);

#endif // PBRT_CORE_TEXTURE_H