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
    cout << "sphere center: " << sCent.x << " " << sCent.y << " " << sCent.z << endl;

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

