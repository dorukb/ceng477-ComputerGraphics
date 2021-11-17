#include "helper.h"
#include <math.h>

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

Vec3f substract(const Vec3f &a,const Vec3f &b){

    Vec3f result;
    result.x = b.x - a.x;
    result.y = b.y - a.y;
    result.z = b.z - a.z;

    return result;
}

Vec3f multVector(const Vec3f &a,const Vec3f &b){

    Vec3f result;
    result.x = a.x * b.x;
    result.y = a.y * b.y;
    result.z = a.z * b.z;

    return result;
}
double dotProduct(const Vec3f &a, Vec3f &b)
{
    double result = a.x*b.x + a.y*b.y + a.z*b.z;
    return result;
}

double length(Vec3f a)
{
    return sqrt((a.x*a.x)+(a.y*a.y)+(a.z*a.z));
}

Vec3f unitVector(Vec3f a)
{
    Vec3f result;
    double l;
    l = length(a);
    result.x = a.x/l;
    result.y = a.y/l;
    result.z = a.z/l;
    return result;
}
