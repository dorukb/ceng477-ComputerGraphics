#include <iostream>
// #include "parser.h"
#include "ppm.h"
#include "helper.h"
#include <math.h>

typedef unsigned char RGB[3];

using namespace parser;
using namespace std;

Ray generateRay(int i, int j, Camera cam);
float intersectSphere(Ray r, Sphere s);
Vec3f computeColor(Ray r);
float determinant(float m [3][3]);
bool intersectTriangle(Ray ray, Triangle tri);

parser::Scene scene;


int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    // parser::Scene scene;
    scene.loadFromXml(argv[1]);

//compute eye rays
//ray objesi olucak
//intersection yaz
///color func
//delta fonk
//trişange intersection algo bak
//shade
//mirror
    
    // for each camera do
    //  for each pixel in image plane, do:
    //      compute viewing (eye, primary) rays
    //      find the first object hit by ray and its surface normal n
    //      set pixel color to value computed from hit point, light, and n


    // ·NearPlane attribute defines the coordinates of the image plane with Left, Right, Bottom,
    // Top floating point parameters, respectively.
    // ·NearDistance defines the distance of the image plane to the camera.
    vector<Camera> cameras = scene.cameras;

    // Debug
    Vec3f sCent = scene.vertex_data[scene.spheres[0].center_vertex_id - 1];
    // cout << "sphere center: " << sCent.x << " " << sCent.y << " " << sCent.z << endl;

    for(int camIndex = 0; camIndex < cameras.size(); camIndex++)
    {
        Camera currCam = cameras[camIndex];
        unsigned char* image = new unsigned char [currCam.image_width * currCam.image_height * 3];

     
        for(int j = 0; j < currCam.image_height; j++)
        {
            for(int i = 0; i < currCam.image_width; i++)
            {
                Ray ray = generateRay(i, j, currCam); 

                // std::cout << "Ray dir: "<< ray.dir.x << ray.dir.y << ray.dir.z << " o: " << ray.origin.x << ray.origin.y << ray.origin.z << std::endl;
                Vec3f rayColor = computeColor(ray);
                if(rayColor.x > 0.01 || rayColor.y > 0.01|| rayColor.z > 0.01){
                    // cout << "Raycolor:" << rayColor.x << " "<< rayColor.y << " "<< rayColor.z<<endl;
                }
                int imgIndex = 3 * (i + j * currCam.image_width);
                image[imgIndex] = (unsigned char)(rayColor.x*255);
                image[imgIndex + 1] = (unsigned char)(rayColor.y*255);
                image[imgIndex + 2]= (unsigned char)(rayColor.z*255);
			
            }            
        }
        write_ppm(currCam.image_name.c_str(), image, currCam.image_width, currCam.image_height);

    }



}


bool intersectTriangle(Ray ray, Triangle tri)
{
    Vec3f a = scene.vertex_data[tri.indices.v0_id - 1];
    Vec3f b = scene.vertex_data[tri.indices.v1_id - 1];
    Vec3f c = scene.vertex_data[tri.indices.v2_id - 1];

    float matrixA [3][3] = {a.x - b.x, a.x - c.x, ray.dir.x,
                            a.y - b.y, a.y - c.y, ray.dir.y,
                            a.z - b.z, a.z - c.z, ray.dir.z};
    float detA = determinant(matrixA);
    if(detA == 0) return false;

    // Cramers Rule
    float matrixBeta [3][3] = {a.x - ray.origin.x, a.x - c.x, ray.dir.x,
                               a.y - ray.origin.y, a.y - c.y, ray.dir.y,
                               a.z - ray.origin.z, a.z - c.z, ray.dir.z};
    float beta = determinant(matrixBeta) / detA;                      
    if(beta < 0) return false;

    float matrixGama [3][3] = {a.x - b.x, a.x - ray.origin.x, ray.dir.x,
                               a.y - b.y, a.y - ray.origin.y, ray.dir.y,
                               a.z - b.z, a.z - ray.origin.z, ray.dir.z};
    float gama = determinant(matrixGama) / detA;
    if(gama < 0 || gama + beta > 1) return false;

    float matrixT[3][3] = {a.x - b.x, a.x - c.x, a.x - ray.origin.x,
                           a.y - b.y, a.y - c.y, a.y - ray.origin.y,
                           a.z - b.z, a.z - c.z, a.z - ray.origin.z};
    float t = determinant(matrixT) / detA;
    if(t <= 0 || t > 99999)
    {
        return false;
    }
    // all conditions satisfied, it does intersect.
    return true;
}
float determinant(float m [3][3])
{
    float firstTerm =  m[0][0] * (m[1][1] * m[2][2] - m[1][2]* m[2][1]);
    float secondTerm = m[1][0] * (m[0][2] * m[2][1] - m[0][1]* m[2][2]);
    float thirdTerm =  m[2][0] * (m[0][1] * m[1][2] - m[1][1]* m[0][2]);
    return firstTerm + secondTerm + thirdTerm;
}
Vec3f computeColor(Ray r)
{
	int i;
	Vec3f c;
	float minT = 90000; // some large number
	float t;
	Vec3f L,N;
	Vec3f P;
	int minI;
	
	c.x=c.y=c.z=0;
	minI = -1;

    vector<Sphere> spheres = scene.spheres;
	for (i = 0; i < spheres.size(); i++)
	{
		t = intersectSphere(r,spheres[i]);
		if (t<minT && t>=0)
		{
            // Material mat = scene.materials[spheres[i].material_id - 1];
            Vec3f tempColor = {1,0,0};
			c = tempColor;// can be replaced with any material property
			minI = i;
			minT = t;
		}
	}

    // triangle intersections
    vector<Triangle> triangles = scene.triangles;
    for (i = 0; i < triangles.size(); i++)
	{
		if (intersectTriangle(r, triangles[i]))
		{
            // Material mat = scene.materials[spheres[i].material_id - 1];
            Vec3f tempColor = {1,0,0};
			c = tempColor;// can be replaced with any material property
			minI = i;
			minT = t;
		}
	}
	// if (minI!=-1)
	// {
	// 	P = add(r.o,multS(r.d,minT));
	// 	L = add(light,multS(P,-1));
	// 	N = add(P,multS(spheres[minI].center,-1));
	// 	L = normalize(L);
	// 	N = normalize(N);
	// 	c = multS(c,dot(L,N));
	// }
	return c;
}

float intersectSphere(Ray r, Sphere s)
{
	float A,B,C; //constants for the quadratic equation
	float delta;
	Vec3f c = scene.vertex_data[s.center_vertex_id - 1];
    float t,t1,t2;
	
	C = (r.origin.x-c.x)*(r.origin.x-c.x)+(r.origin.y-c.y)*(r.origin.y-c.y)+(r.origin.z-c.z)*(r.origin.z-c.z)-s.radius*s.radius;
	B = 2*r.dir.x*(r.origin.x-c.x)+2*r.dir.y*(r.origin.y-c.y)+2*r.dir.z*(r.origin.z-c.z);
	A = r.dir.x*r.dir.x+r.dir.y*r.dir.y+r.dir.z*r.dir.z;
	delta = B*B-4*A*C;
	
	if (delta<0){
        return -1;
    } 
	else if (delta==0)
	{
		t = -B / (2*A);
	}
	else
	{
		delta = sqrt(delta);
		A = 2*A;
		t1 = (-B + delta) / A;
		t2 = (-B - delta) / A;
				
		if (t1<t2) t=t1; else t=t2;
	}
	
	return t;
}
Ray generateRay(int i, int j, Camera cam)
{
	Ray result;
	float su,sv;
	Vec3f m,q,s;

    // NearPlane: coords of image plane with Left, Right, Bottom, Top floats respectively.
    float left = cam.near_plane.x;
	float right = cam.near_plane.y;
    float bottom = cam.near_plane.z;
    float top = cam.near_plane.w;

    float nx = cam.image_width;
    float ny = cam.image_height;

	su = (i +0.5)*(right-left)/ nx;
	sv = (j+0.5)*(top-bottom)/ ny;
	
    Vec3f e = cam.position;
    Vec3f gaze = cam.gaze;
    float dist = cam.near_distance;

    // Up = v,  Gaze = −w, u = v ×w
    Vec3f v = cam.up;
    Vec3f w = multScaler(gaze,-1);

    Vec3f u = cross(v, w);

	m = add(e,multScaler(gaze,dist));
	q = add(m,add(multScaler(u,left),multScaler(v,top)));
	s = add(q,add(multScaler(u,su),multScaler(v,-sv)));
	
	result.origin= e;
	result.dir = add(s,multScaler(e,-1));
	
	return result;
}

