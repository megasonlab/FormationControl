#ifndef _VECTOR_H
#define _VECTOR_H

//Vector Calculation Operator
//http://www3.u-toyama.ac.jp/akiyama/

#include <math.h>

typedef struct {
    double x;
    double y;
} Vector2D;

const static Vector2D Zero2D = {0.0, 0.0};

typedef struct {
    double x;
    double y;
    double z;
} Vector3D;

const static Vector3D Zero3D = {0.0, 0.0, 0.0};

//3D Vector Calculation Operators
inline Vector3D operator + (Vector3D u,Vector3D v) {
	Vector3D w = { u.x + v.x, u.y + v.y, u.z + v.z };
	return w;
}

inline Vector3D& operator += (Vector3D& u,Vector3D v) {
    u.x += v.x; u.y += v.y; u.z += v.z;
    return u;
}

inline Vector3D operator - (Vector3D u,Vector3D v) {
	Vector3D w = { u.x - v.x, u.y - v.y, u.z - v.z, };
	return w;
}

inline Vector3D& operator -= (Vector3D& u,Vector3D v) {
    u.x -= v.x; u.y -= v.y; u.z -= v.z;
    return u;
}

inline Vector3D operator * (double a,Vector3D u) {
	Vector3D w = { a * u.x, a * u.y, a * u.z };
	return w;
}

inline Vector3D operator * (Vector3D u,double a) {
	Vector3D w = { a * u.x, a * u.y, a * u.z };
	return w;
}

inline Vector3D& operator - (Vector3D& u) {
	u.x *= - 1; u.y *= - 1; u.z *= - 1;
	return u;
}

inline Vector3D operator / (Vector3D u,double a) {
	Vector3D w = { u.x / a, u.y / a, u.z / a };
	return w;
}

inline Vector3D& operator /= (Vector3D& u,double a) {
	u.x /= a; u.y /= a;  u.z /= a;
	return u;
}

inline double operator * (Vector3D u,Vector3D v) {
	return u.x * v.x + u.y * v.y + u.z * v.z;
}

inline Vector3D operator ^ (Vector3D u,Vector3D v) {
	Vector3D w =
	{ u.y * v.z - v.y * u.z, u.z * v.x - v.z * u.x, u.x * v.y - v.x * u.y };
	return w;
}

inline double operator ~ (Vector3D u) {
	return sqrt(u.x * u.x + u.y * u.y + u.z * u.z);
}

//2D Vector Calculation Operators
inline Vector2D operator + (Vector2D u,Vector2D v) {
	Vector2D w = { u.x + v.x, u.y + v.y };
	return w;
}

inline Vector2D& operator += (Vector2D& u,Vector2D v) {
    u.x += v.x; u.y += v.y;
    return u;
}

inline Vector2D operator - (Vector2D u,Vector2D v) {
	Vector2D w =  { u.x - v.x, u.y - v.y };
	return w;
}

inline Vector2D& operator -= (Vector2D& u,Vector2D v) {
    u.x -= v.x; u.y -= v.y;
    return u;
}

inline Vector2D operator * (double a,Vector2D u) {
	Vector2D w = { a * u.x, a * u.y };
	return w;
}

inline Vector2D operator * (Vector2D u,double a) {
	Vector2D w = { a * u.x, a * u.y };
	return w;
}

inline Vector2D& operator - (Vector2D& u) {
	u.x *= - 1; u.y *= - 1;
	return u;
}

inline Vector2D operator / (Vector2D u,double a) {
	Vector2D w = { u.x / a, u.y / a };
	return w;
}

inline Vector2D& operator /= (Vector2D& u,double a) {
	u.x /= a; u.y /= a;
	return u;
}

inline double operator * (Vector2D u,Vector2D v) {
	return u.x * v.x + u.y * v.y;
}

inline double operator ^ (Vector2D u,Vector2D v) {
	double  w;
    w = u.x * v.y - v.x * u.y;
	return w;
}

inline double operator ~ (Vector2D u) {
	return sqrt(u.x * u.x + u.y * u.y);
}

#endif

