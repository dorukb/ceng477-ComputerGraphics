#ifndef __helper_h__
#define __helper_h__

#include "parser.h"
#include "ppm.h"
using namespace parser;

Vec3f multScaler(const Vec3f &a,float s);
Vec3f add(const Vec3f &a, const Vec3f &b);
Vec3f cross(const Vec3f &first, const Vec3f &second);

struct Ray
{
	Vec3f origin, dir;
};
Vec3f substract(const Vec3f &a,const Vec3f &b);

Vec3f multVector(const Vec3f &a,const Vec3f &b);

double dotProduct(const Vec3f &a, Vec3f &b);
double length(Vec3f a);
Vec3f unitVector(Vec3f a);
float determinant(float m [3][3]);




#endif // __helper_h__
//compute eye rays
//ray objesi olucak
//intersection yaz
///color func
//delta fonk
//trişange intersection algo bak
//shade
//mirror


