#ifndef __helper_h__
#define __helper_h__

#include "parser.h"
#include "ppm.h"
using namespace parser;

struct Ray
{
	Vec3f origin, dir;
};

Vec3f multScaler(const Vec3f &a,float s);
Vec3f multVector(const Vec3f &a,const Vec3f &b);
Vec3f cross(const Vec3f &first, const Vec3f &second);

float dotProduct(const Vec3f &a, Vec3f &b);
float length(Vec3f a);
Vec3f makeUnitVector(Vec3f a);
float determinant(float m [3][3]);
Vec3f operator+(const parser::Vec3f& v1, const parser::Vec3f& v2);
Vec3f operator-(const parser::Vec3f& v1, const parser::Vec3f& v2);
Vec3i clamp(Vec3f colors);

// Vec3f operator*(const parser::Vec3f& v1, float alpha){
//     Vec3f result;
//     result.x = v1.x * alpha;
//     result.y = v1.y * alpha;
//     result.z = v1.z * alpha;

//     return result;
// }

#endif // __helper_h__

