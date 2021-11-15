#include "helper.h"

using namespace parser;

Vec3f multScaler(const Vec3f &a,float s)
{
	Vec3f result;
	result.x = a.x*s;
	result.y = a.y*s;
	result.z = a.z*s;
	return result;
}

Vec3f add(const Vec3f &a, const Vec3f &b)
{
	Vec3f result;
	result.x = a.x+b.x;
	result.y = a.y+b.y;
	result.z = a.z+b.z;
	return result;
}

Vec3f cross(const Vec3f &first, const Vec3f &second){
    Vec3f result;
    result.x = first.y * second.z - first.z * second.y;
    result.y = first.z * second.x - first.x * second.z;
    result.z = first.x * second.y - first.y * second.x;
    return result;
}