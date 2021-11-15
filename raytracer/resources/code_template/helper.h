#ifndef __helper_h__
#define __helper_h__

typedef struct
{
	float x,y,z;
} vec3f;

typedef struct
{
	vec3f center;
	float radius;
	vec3f color;
} sphere;


typedef struct
{
	int x,y,z;
} vec3i;
typedef struct
{
	vec3f o,d;
} ray;

vec3i **image;

#endif // __helper_h__
//compute eye rays
//ray objesi olucak
//intersection yaz
///color func
//delta fonk
//trişange intersection algo bak
//shade
//mirror


