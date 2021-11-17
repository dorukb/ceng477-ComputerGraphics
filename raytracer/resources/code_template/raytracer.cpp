#include <iostream>
// #include "parser.h"
#include "ppm.h"
#include "helper.h"
#include <math.h>

#define INF 999999.0

typedef unsigned char RGB[3];

using namespace parser;
using namespace std;

Ray generateRay(int i, int j, Camera cam);
float intersectSphere(Ray r, Sphere s);
Vec3f computeColor(Ray r);
float determinant(float m [3][3]);
float intersectTriangle(Ray ray, Triangle tri);

parser::Scene scene;
int t_min;

//L_a = k_a * l_a
Vec3f L_a(Scene scene,int i){
    return multVector(scene.materials[i].ambient,scene.ambient_light);
}
//L_d  = k_d * cos(the)' * E_i
//cos(the)' = max (0,w_i.n)
Vec3f L_d(Scene scene,int i,Vec3f w_i,Vec3f n,Vec3f E_i){
    Vec3f result;
    float cos_theta = max(0.0,dotProduct(w_i,n));
    Vec3f k_d = scene.materials[i].diffuse;
    result = multVector(multScaler(k_d,cos_theta),E_i);
    return result;
}
//L_s =k_s * cos(alp)' ^p * E_i
//cos(alp)' = max (0,n.h)
//vec3f h = unitVector(add(w_i, w_o));
Vec3f L_s(Scene scene,int i,Vec3f h,Vec3f n,Vec3f E_i){
    Vec3f result;
    Vec3f k_s = scene.materials[i].specular;
    float cos_alpha = max(0.0, dotProduct(n, h));
    float p = scene.materials[i].phong_exponent;
    double exp = pow(cos_alpha,p);
    result = multVector(multScaler(k_s,exp),E_i);
    return result;
}



int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    // parser::Scene scene;
    scene.loadFromXml(argv[1]);
    vector<Camera> cameras = scene.cameras;
    int camSize = cameras.size();
    int triangelsSize = scene.triangles.size();
    int meshesSize = scene.meshes.size();
    int spheresSize = scene.spheres.size();
    int lightSize = scene.point_lights.size();
    int depth = scene.max_recursion_depth;
    Vec3i backround = scene.background_color;
    
     

    
    // for each camera do
    //  for each pixel in image plane, do:
    //      compute viewing (eye, primary) rays
    //      find the first object hit by ray and its surface normal n
    //      set pixel color to value computed from hit point, light, and n
    
    // ·NearPlane attribute defines the coordinates of the image plane with Left, Right, Bottom,
    // Top floating point parameters, respectively.
    // ·NearDistance defines the distance of the image plane to the camera.
    
    
    // Debug
   // Vec3f sCent = scene.vertex_data[scene.spheres[0].center_vertex_id - 1];
    // cout << "sphere center: " << sCent.x << " " << sCent.y << " " << sCent.z << endl;
   
    for(int camIndex = 0; camIndex < camSize; camIndex++)
    {
        Camera currCam = cameras[camIndex];
        int imageWidth = cameras[camIndex].image_height;
        int imageHeight = cameras[camIndex].image_width;
        
        
        unsigned char* image = new unsigned char [imageWidth * imageHeight * 3];
        

     //image produces
        for(int j = 0; j < imageHeight; j++)
        {
            for(int i = 0; i < imageWidth; i++)
            {
                //RAY TRACING
                //compute viewving ray from e to s
                Ray ray = generateRay(i, j, currCam);
                int closest_S_T_M[3]={-1,-1,-1};// [sphere_index,triange_index,mesh_index]
                Vec3f x;
                t_min = INF;
                double t ;
                Vec3f shade ;
                /*
                
            //for sphere
                //for each object o:
                for(int o=0;o<spheresSize;o++){
                    //if ray intersects o at point x:
                    t = intersectSphere(ray,scene.spheres[o]);
                    // if t<tmin:
                    //  tmin= t,obj = o;
                    if(t>=0){
                        if(t<t_min){
                            t_min = t;
                            closest_S_T_M[0]=o;
                            closest_S_T_M[1]=-1;
                            closest_S_T_M[2]=-1;
                            //calculate intersection point x;
                            x =  add(ray.origin ,multScaler(ray.dir,t_min));
                        }
                    }
                    
                    // if t<tmin:
                    //  tmin= t,obj = o;
                    
                }
            // for triange
                for(int o=0;o<triangelsSize;o++){
                    //if ray intersects o at point x:
                    //t = intersectTriangle(ray,scene.spheres[o]);
                    // if t<tmin:
                    //  tmin= t,obj = o;
                    if(t>=0){
                        if(t<t_min){
                            t_min = t;
                            closest_S_T_M[0]=-1;
                            closest_S_T_M[1]=o;
                            closest_S_T_M[2]=-1;
                            //calculate intersection point x;
                            x =  add(ray.origin ,multScaler(ray.dir,t_min));
                        }
                    }
                    // if t<tmin:
                    //  tmin= t,obj = o;
                    
                }
                
            // for meshes
                for(int o=0;o<meshesSize;o++){
                    //if ray intersects o at point x:
                   // t = intersectMesh(ray,scene.spheres[o]);
                    // if t<tmin:
                    //  tmin= t,obj = o;
                    if(t>=0){
                        if(t<t_min){
                            t_min = t;
                            closest_S_T_M[0]=-1;
                            closest_S_T_M[1]=-1;
                            closest_S_T_M[2]=o;
                            //calculate intersection point x;
                            x = add(ray.origin ,multScaler(ray.dir,t_min));
                        }
                    }
                    
                        
                   
                    
                }
                Vec3f pixel_color ;
                //if object is not null:
                if(triangelsSize > 0 || spheresSize > 0 || meshesSize > 0 ){
                    //  pixel color = La ->ambiant shading is not effected by shadosw
                    if( closest_S_T_M[0] != -1){//closest is a sphere
                        pixel_color = L_a(scene,scene.spheres[closest_S_T_M[0]].material_id);
                        //  for easch light l:
                        for(int l=0;l<lightSize;l++){
                            
                        }
                        
                    }
                   
                    
                    
                    //      compute shadow ray s from x to l:
                    //      for each object p:
                    //          if s intersects p before the light source:
                    //              continue the light loop
                    //      pixel color = pixel color + Ld +Ls ->diffuse and specular components
                }
                
            
                
                else{
                    //  pixel color = color of backround
                }
                
                
              */
                
                
                // do not forget clamping the pixel value to [0,255] range and rounding it to the nearest integer
                
                
               // Vec3f shaded_pixel =
                
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
        write_ppm(currCam.image_name.c_str(), image, imageWidth, imageHeight);

    }

    




    

    
    


}

float intersectFace(Ray ray, Face face)
{
    Vec3f a = scene.vertex_data[face.v0_id - 1];
    Vec3f b = scene.vertex_data[face.v1_id - 1];
    Vec3f c = scene.vertex_data[face.v2_id - 1];

    float matrixA [3][3] = {a.x - b.x, a.x - c.x, ray.dir.x,
                            a.y - b.y, a.y - c.y, ray.dir.y,
                            a.z - b.z, a.z - c.z, ray.dir.z};
    float detA = determinant(matrixA);
    if(detA == 0) return -1;

    // Cramers Rule
    float matrixBeta [3][3] = {a.x - ray.origin.x, a.x - c.x, ray.dir.x,
                               a.y - ray.origin.y, a.y - c.y, ray.dir.y,
                               a.z - ray.origin.z, a.z - c.z, ray.dir.z};
    float beta = determinant(matrixBeta) / detA;                      
    if(beta < 0) return -1;

    float matrixGama [3][3] = {a.x - b.x, a.x - ray.origin.x, ray.dir.x,
                               a.y - b.y, a.y - ray.origin.y, ray.dir.y,
                               a.z - b.z, a.z - ray.origin.z, ray.dir.z};
    float gama = determinant(matrixGama) / detA;
    if(gama < 0 || gama + beta > 1) return -1;

    float matrixT[3][3] = {a.x - b.x, a.x - c.x, a.x - ray.origin.x,
                           a.y - b.y, a.y - c.y, a.y - ray.origin.y,
                           a.z - b.z, a.z - c.z, a.z - ray.origin.z};
    float t = determinant(matrixT) / detA;
    return t;  
}


float intersectTriangle(Ray ray, Triangle tri)
{
    Vec3f a = scene.vertex_data[tri.indices.v0_id - 1];
    Vec3f b = scene.vertex_data[tri.indices.v1_id - 1];
    Vec3f c = scene.vertex_data[tri.indices.v2_id - 1];

    float matrixA [3][3] = {a.x - b.x, a.x - c.x, ray.dir.x,
                            a.y - b.y, a.y - c.y, ray.dir.y,
                            a.z - b.z, a.z - c.z, ray.dir.z};
    float detA = determinant(matrixA);
    if(detA == 0) return -1;

    // Cramers Rule
    float matrixBeta [3][3] = {a.x - ray.origin.x, a.x - c.x, ray.dir.x,
                               a.y - ray.origin.y, a.y - c.y, ray.dir.y,
                               a.z - ray.origin.z, a.z - c.z, ray.dir.z};
    float beta = determinant(matrixBeta) / detA;                      
    if(beta < 0) return -1;

    float matrixGama [3][3] = {a.x - b.x, a.x - ray.origin.x, ray.dir.x,
                               a.y - b.y, a.y - ray.origin.y, ray.dir.y,
                               a.z - b.z, a.z - ray.origin.z, ray.dir.z};
    float gama = determinant(matrixGama) / detA;
    if(gama < 0 || gama + beta > 1) return -1;

    float matrixT[3][3] = {a.x - b.x, a.x - c.x, a.x - ray.origin.x,
                           a.y - b.y, a.y - c.y, a.y - ray.origin.y,
                           a.z - b.z, a.z - c.z, a.z - ray.origin.z};
    float t = determinant(matrixT) / detA;
    return t;  
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
        t = intersectTriangle(r, triangles[i]);
		if (t > 0 && t<minT)
		{
            // Material mat = scene.materials[spheres[i].material_id - 1];
            Vec3f tempColor = {1,0,0};
			c = tempColor;// can be replaced with any material property
			minI = i;
			minT = t;
		}
	}
    
    
    vector<Mesh> meshes = scene.meshes;
    for(i = 0; i < meshes.size(); i++)
    {
        vector<Face> faces = meshes[i].faces;
        for(int j = 0; j < faces.size(); j++)
        {
            t = intersectFace(r, faces[j]);
            if (t > 0 && t<minT)
            {
                // Material mat = scene.materials[spheres[i].material_id - 1];
                Vec3f tempColor = {1,0,0};
                c = tempColor;// can be replaced with any material property
                minI = i;
                minT = t;
            }
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
	
	if (delta<0) return -1;
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




//compute eye rays
//ray objesi olucak
//intersection yaz
///color func
//delta fonk
//trişange intersection algo bak
//shade
//mirror
