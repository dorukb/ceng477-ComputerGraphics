#include "helper.h"
#include "ppm.h"
#include <math.h>
#include <iostream>

using namespace parser;


Vec3f multScaler(const Vec3f &a,float s)
{
	Vec3f result;
	result.x = a.x*s;
	result.y = a.y*s;
	result.z = a.z*s;
	return result;
}

Vec3f cross(const Vec3f &first, const Vec3f &second){
    Vec3f result;
    result.x = first.y * second.z - first.z * second.y;
    result.y = first.z * second.x - first.x * second.z;
    result.z = first.x * second.y - first.y * second.x;
    return result;
}


Vec3f operator+(const parser::Vec3f& v1, const parser::Vec3f& v2){
    Vec3f result;
    result.x = v1.x + v2.x;
    result.y = v1.y + v2.y;
    result.z = v1.z + v2.z;

    return result;
}

Vec3f operator-(const parser::Vec3f& v1, const parser::Vec3f& v2){
    Vec3f result;
    result.x = v1.x - v2.x;
    result.y = v1.y - v2.y;
    result.z = v1.z - v2.z;

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

Vec3f makeUnitVector(Vec3f a)
{
    Vec3f result;
    double l;
    l = length(a);
    result.x = a.x/l;
    result.y = a.y/l;
    result.z = a.z/l;
    return result;
}

float determinant(float m [3][3])
{
    float firstTerm =  m[0][0] * (m[1][1] * m[2][2] - m[1][2]* m[2][1]);
    float secondTerm = m[1][0] * (m[0][2] * m[2][1] - m[0][1]* m[2][2]);
    float thirdTerm =  m[2][0] * (m[0][1] * m[1][2] - m[1][1]* m[0][2]);
    return firstTerm + secondTerm + thirdTerm;
}

Vec3f clamp(Vec3f colors)
{
    Vec3f result;
    if(colors.x > 255){ 
        result.x = 255;
         }
    else if(colors.x < 0) { 
        result.x = 0;
         }
    else {
         result.x = (int) round(colors.x);
          }
    if(colors.y > 255){ 
        result.y = 255;
         }
    else if(colors.y < 0) { 
        result.y = 0; 
        }
    else { 
        result.y = (int) round(colors.y);
         }
    if(colors.z > 255){ 
        result.z = 255; 
        }
    else if(colors.z < 0) { 
        result.z = 0;
         }
    else { 
        result.z = (int) round(colors.z); 
        }
    return result;
}