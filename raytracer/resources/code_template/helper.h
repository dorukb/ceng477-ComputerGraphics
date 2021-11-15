#ifndef __helper_h__
#define __helper_h__

#include "parser.h"
using namespace parser;

Vec3f multScaler(const Vec3f &a,float s);
Vec3f add(const Vec3f &a, const Vec3f &b);
Vec3f cross(const Vec3f &first, const Vec3f &second);

struct Ray
{
	Vec3f origin, dir;
};

#endif // __helper_h__
//compute eye rays
//ray objesi olucak
//intersection yaz
///color func
//delta fonk
//trişange intersection algo bak
//shade
//mirror


