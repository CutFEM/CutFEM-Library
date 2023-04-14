/*
This file is part of CutFEM-Library.

CutFEM-Library is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

CutFEM-Library is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
CutFEM-Library. If not, see <https://www.gnu.org/licenses/>
*/
#ifndef COMMON_POINT_HPP
#define COMMON_POINT_HPP
#include <cmath>
#include <algorithm>
#include <vector>
#include <span>

class R0 {
  public:
    using R = double;

    static const int d = 0;
    R0() {}
};

class R1 {
  public:
    using R = double;

    static const int d = 1;
    static const std::vector<R1> KHat;

    R x;
    R1() : x(0.) {}
    R1(R a) : x(a) {}
    R1(R *a) : x(a[0]) {}
    R1(std::span<R> a) : x(a[0]) {}
    R1(const R1 &a, const R1 &b) : x(b.x - a.x) {}
    operator double() const { return x; }
    operator double *() { return &x; }
    operator std::span<R>() { return std::span<R>(&x, 1); }

    operator const double *() const { return &x; }
    const R X() const { return x; }
    R &X() { return x; }
    R1 &operator=(const R *P) {
        x = P[0];
        return *this;
    }
    R1 &operator+=(const R1 &P) {
        x += P.x;
        return *this;
    }
    R1 &operator-=(const R1 &P) {
        x -= P.x;
        return *this;
    }
    R1 &operator*=(R a) {
        x *= a;
        return *this;
    }
    R1 &operator/=(R a) {
        x /= a;
        return *this;
    }
    R1 operator+(const R1 &P) const { return R1(x + P.x); }
    R1 operator-(const R1 &P) const { return R1(x - P.x); }
    R operator,(const R1 &P) const { return x * P.x; }
    R1 operator*(R c) const { return R1(x * c); }
    R1 operator/(R c) const { return R1(x / c); }
    R1 operator-() const { return R1(-x); }
    R1 operator+() const { return *this; }

    R sum() const { return x; }
    R &operator[](int i) { return (&x)[i]; }
    const R &operator[](int i) const { return (&x)[i]; }
    friend R1 operator*(R c, const R1 &P) { return P * c; }

    friend std::ostream &operator<<(std::ostream &f, const R1 &P) {
        f << P.x;
        return f;
    }
    friend std::istream &operator>>(std::istream &f, R1 &P) {
        f >> P.x;
        return f;
    }
};

class R2 {
  public:
    using R = double;

    static const int d = 2;
    static const std::vector<R2> KHat;
    R x, y;

    R2() : x(0.), y(0.) {}
    R2(R a, R b) : x(a), y(b) {}
    R2(const R *a) : x(a[0]), y(a[1]) {}
    R2(R *a) : x(a[0]), y(a[1]) {}
    R2(std::span<R> a) : x(a[0]), y(a[1]) {}
    R2(const R2 &a, const R2 &b) : x(b.x - a.x), y(b.y - a.y) {}

    operator const double *() const { return &x; }
    operator double *() { return &x; }
    operator std::span<R>() { return std::span<R>(&x, 2); }

    R2 &operator=(const R *P) {
        x = P[0];
        y = P[1];
        return *this;
    }
    R2 &operator+=(const R2 &P) {
        x += P.x;
        y += P.y;
        return *this;
    }
    R2 &operator-=(const R2 &P) {
        x -= P.x;
        y -= P.y;
        return *this;
    }
    R2 &operator*=(R a) {
        x *= a;
        y *= a;
        return *this;
    }
    R2 &operator/=(R a) {
        x /= a;
        y /= a;
        return *this;
    }

    R2 operator+(const R2 &P) const { return R2(x + P.x, y + P.y); }
    R2 operator-(const R2 &P) const { return R2(x - P.x, y - P.y); }
    R operator,(const R2 &P) const { return x * P.x + y * P.y; }
    R operator^(const R2 &P) const { return x * P.y - y * P.x; }
    R2 operator*(R c) const { return R2(x * c, y * c); }
    R2 operator/(R c) const { return R2(x / c, y / c); }
    R2 operator-() const { return R2(-x, -y); }
    R2 operator+() const { return *this; }

    R2 perp() const { return R2(-y, x); }
    R sum() const { return x + y; }

    R &operator[](int i) { return (&x)[i]; }
    const R &operator[](int i) const { return (&x)[i]; }

    R X() const { return x; }
    R Y() const { return y; }
    R Z() const { return 0.; }

    R norme() const { return std::sqrt(x * x + y * y); }
    R norm() const { return std::sqrt(x * x + y * y); }
    R norme2() const { return (x * x + y * y); }

    friend R2 operator*(R c, const R2 &P) { return P * c; }
    friend R2 perp(const R2 &P) { return R2(-P.y, P.x); }
    friend R det(const R2 &A, const R2 &B, const R2 &C) { return R2(A, B) ^ R2(A, C); }

    friend std::ostream &operator<<(std::ostream &f, const R2 &P) {
        f << P.x << ' ' << P.y;
        return f;
    }
    friend std::istream &operator>>(std::istream &f, R2 &P) {
        f >> P.x >> P.y;
        return f;
    }
};

inline double Norme_infty(const R2 &A) { return std::max(std::fabs(A.x), std::fabs(A.y)); }
inline double Norme2_2(const R2 &A) { return (A, A); }
inline double Norme2(const R2 &A) { return std::sqrt((A, A)); }

inline double norm2(const R2 &A) { return std::sqrt((A, A)); }
inline double norm2_2(const R2 &A) { return (A, A); }

class R3 {

  public:
    using R = double;

    static const int d = 3;
    static const std::vector<R3> KHat;

    R x, y, z;

    R3() : x(0), y(0), z(0){};
    R3(R a, R b, R c) : x(a), y(b), z(c) {}
    R3(const R *a) : x(a[0]), y(a[1]), z(a[2]) {}
    R3(std::span<R> a) : x(a[0]), y(a[1]), z(a[2]) {}

    R3(R2 P2) : x(P2.x), y(P2.y), z(0) {}
    R3(R2 P2, R zz) : x(P2.x), y(P2.y), z(zz) {}
    R3(const R3 &A, const R3 &B) : x(B.x - A.x), y(B.y - A.y), z(B.z - A.z) {}
    static R3 diag(R a) { return R3(a, a, a); }
    operator const double *() const { return &x; }
    operator double *() { return &x; }
    operator std::span<R>() { return std::span<R>(&x, 3); }

    R3 &operator=(const R *P) {
        x = P[0];
        y = P[1];
        z = P[2];
        return *this;
    }
    R3 &operator=(const R2 &P2) {
        x = P2.x;
        y = P2.y;
        z = 0;
        return *this;
    }
    R3 operator+(const R3 &P) const { return R3(x + P.x, y + P.y, z + P.z); }
    R3 &operator+=(const R3 &P) {
        x += P.x;
        y += P.y;
        z += P.z;
        return *this;
    }
    R3 operator-(const R3 &P) const { return R3(x - P.x, y - P.y, z - P.z); }
    R3 &operator-=(const R3 &P) {
        x -= P.x;
        y -= P.y;
        z -= P.z;
        return *this;
    }
    R3 &operator*=(R c) {
        x *= c;
        y *= c;
        z *= c;
        return *this;
    }
    R3 &operator/=(R c) {
        x /= c;
        y /= c;
        z /= c;
        return *this;
    }
    R3 operator-() const { return R3(-x, -y, -z); }
    R3 operator+() const { return *this; }
    R operator,(const R3 &P) const { return x * P.x + y * P.y + z * P.z; }
    R3 operator^(const R3 &P) const { return R3(y * P.z - z * P.y, P.x * z - x * P.z, x * P.y - y * P.x); }
    R3 operator*(R c) const { return R3(x * c, y * c, z * c); }
    R3 operator/(R c) const { return R3(x / c, y / c, z / c); }
    R &operator[](int i) { return (&x)[i]; }
    const R &operator[](int i) const { return (&x)[i]; }
    friend R3 operator*(R c, const R3 &P) { return P * c; }
    friend R3 operator/(R c, const R3 &P) { return P / c; }

    R norme() const { return std::sqrt(x * x + y * y + z * z); }
    R norm() const { return std::sqrt(x * x + y * y + z * z); }
    R norme2() const { return (x * x + y * y + z * z); }
    R sum() const { return x + y + z; }

    R X() const { return x; }
    R Y() const { return y; }
    R Z() const { return z; }

    friend std::ostream &operator<<(std::ostream &f, const R3 &P) {
        f << P.x << ' ' << P.y << ' ' << P.z;
        return f;
    }
    friend std::istream &operator>>(std::istream &f, R3 &P) {
        f >> P.x >> P.y >> P.z;
        return f;
    }

    friend R det(R3 A, R3 B, R3 C) {
        R s = 1.;
        if (fabs(A.x) < fabs(B.x))
            std::swap(A, B), s = -s;
        if (fabs(A.x) < fabs(C.x))
            std::swap(A, C), s = -s;
        if (fabs(A.x) > 1e-50) {
            s *= A.x;
            A.y /= A.x;
            A.z /= A.x;
            B.y -= A.y * B.x;
            B.z -= A.z * B.x;
            C.y -= A.y * C.x;
            C.z -= A.z * C.x;
            return s * (B.y * C.z - B.z * C.y);
        } else
            return 0.;
    }

    friend R det(const R3 &A, const R3 &B, const R3 &C, const R3 &D) { return det(R3(A, B), R3(A, C), R3(A, D)); }

    R2 p2() const { return R2(x, y); }
};

inline double Norme_infty(const R3 &A) { return std::max({std::fabs(A.x), std::fabs(A.y), std::fabs(A.z)}); }
inline double Norme2_2(const R3 &A) { return (A, A); }
inline double Norme2(const R3 &A) { return sqrt((A, A)); }

struct lessRd {
    bool operator()(const R1 &s1, const R1 &s2) const { return s1.X() < s2.X(); }
    bool operator()(const R2 &s1, const R2 &s2) const { return s1.x == s2.x ? (s1.y < s2.y) : s1.x < s2.x; }
    bool operator()(const R3 &s1, const R3 &s2) const {
        return s1.x == s2.x ? (s1.y == s2.y ? (s1.z < s2.z) : s1.y < s2.y) : s1.x < s2.x;
    }
};

template <int d> struct typeRd {
    typedef R0 Rd;
};
template <> struct typeRd<1> {
    typedef R1 Rd;
};
template <> struct typeRd<2> {
    typedef R2 Rd;
};
template <> struct typeRd<3> {
    typedef R3 Rd;
};

#endif