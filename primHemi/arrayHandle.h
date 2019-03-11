// A simple 4D vector with a variety of operators and other member functions 
// for doing 4D vector arithmetic. This code uses HEMI to be portable between
// standard host code (compiled with any C++ compiler, including NVCC) and 
// CUDA device code (compiled with NVCC).
#ifndef __ddouble_H__
#define __ddouble_H__

#include "hemi/hemi.h"


struct HEMI_ALIGN(16) ddouble
{
  double val;

  HEMI_DEV_CALLABLE_INLINE_MEMBER
  ddouble(){}

  HEMI_DEV_CALLABLE_INLINE_MEMBER
  ddouble(double app) : val(app) {}

  HEMI_DEV_CALLABLE_INLINE_MEMBER
  ddouble(const ddouble& vv) : val(vv.val) {}

  HEMI_DEV_CALLABLE_INLINE_MEMBER
  ddouble& operator=(const ddouble& v) { 
    val = v.val; 
    return *this;
  }

  HEMI_DEV_CALLABLE_INLINE_MEMBER
  ddouble operator+(const ddouble& v) const {
    return ddouble(val+v.val);
  }
/*
 HEMI_DEV_CALLABLE_INLINE_MEMBER
  ddouble operator-(const Vec4f& v) const {
    return Vec4f(x-v.x, y-v.y, z-v.z, w-v.w); 
  }

  HEMI_DEV_CALLABLE_INLINE_MEMBER
  ddouble operator-() const {
    return Vec4f(-x, -y, -z, -w); 
  }

  HEMI_DEV_CALLABLE_INLINE_MEMBER
  ddouble operator*(const Vec4f& v) const {
    return Vec4f(x*v.x, y*v.y, z*v.z, w*v.w); 
  }

  HEMI_DEV_CALLABLE_INLINE_MEMBER
  ddouble operator*(float s) const {
    return Vec4f(x*s, y*s, z*s, w*s);
  }

  HEMI_DEV_CALLABLE_INLINE_MEMBER
  ddouble& operator+=(const Vec4f& v) {
    x += v.x; y += v.y; z += v.z; w += v.w;
    return *this;
  }

  HEMI_DEV_CALLABLE_INLINE_MEMBER
  ddouble& operator-=(const Vec4f& v) {
    x -= v.x; y -= v.y; z -= v.z; w -= v.w;
    return *this;
  }

  HEMI_DEV_CALLABLE_INLINE_MEMBER
  ddouble& operator*=(const Vec4f& v) {
    x *= v.x; y *= v.y; z *= v.z; w *= v.w;
    return *this;
  }

  HEMI_DEV_CALLABLE_INLINE_MEMBER
  ddouble& operator/=(const Vec4f& v) {
    x /= v.x; y /= v.y; z /= v.z; w /= v.w;
    return *this;
  }

  HEMI_DEV_CALLABLE_INLINE_MEMBER
  ddouble& operator*=(float s) {
    x *= s; y *= s; z *= s; w *= s;
    return *this;
  }

  HEMI_DEV_CALLABLE_INLINE_MEMBER
  ddouble& operator/=(float s) {
    x /= s; y /= s; z /= s; w /= s;
    return *this;
  }*/

  HEMI_DEV_CALLABLE_INLINE_MEMBER
  float dot(const ddouble& v) const {
    return val*v.val;
  }

  HEMI_DEV_CALLABLE_INLINE_MEMBER
  float lengthSqr() const {
    return this->dot(*this);
  }

  HEMI_DEV_CALLABLE_INLINE_MEMBER
  int length() const {
    return sqrt(lengthSqr());
  }

  HEMI_DEV_CALLABLE_INLINE_MEMBER
  float inverseLength(float softening = 0.0f) const {
#ifdef HEMI_DEV_CODE
    return rsqrtf(lengthSqr() + softening);
#else
    return 1.0f / sqrtf(lengthSqr() + softening);
#endif
  }
};

#endif // __ddouble_H__