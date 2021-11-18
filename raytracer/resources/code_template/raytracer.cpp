#include <iostream>
// #include "parser.h"
#include "ppm.h"
#include "helper.h"
#include <math.h>

#define INF 999999.0

typedef unsigned char RGB[3];

using namespace parser;
using namespace std;

parser::Scene scene;
int t_min;

Ray generateRay(int i, int j, Camera cam)
{
    Ray result;
    float su, sv;
    Vec3f m, q, s;

    // NearPlane: coords of image plane with Left, Right, Bottom, Top floats respectively.
    float left = cam.near_plane.x;
    float right = cam.near_plane.y;
    float bottom = cam.near_plane.z;
    float top = cam.near_plane.w;

    float nx = cam.image_width;
    float ny = cam.image_height;

    su = (i + 0.5) * (right - left) / nx;
    sv = (j + 0.5) * (top - bottom) / ny;

    Vec3f e = cam.position;
    Vec3f gaze = cam.gaze;
    float dist = cam.near_distance;

    // Up = v,  Gaze = −w, u = v ×w
    Vec3f v = cam.up;
    Vec3f w = multScaler(gaze, -1);

    Vec3f u = cross(v, w);

    m = add(e, multScaler(gaze, dist));
    q = add(m, add(multScaler(u, left), multScaler(v, top)));
    s = add(q, add(multScaler(u, su), multScaler(v, -sv)));

    result.origin = e;
    result.dir = add(s, multScaler(e, -1));

    return result;
}

float intersectFace(Ray ray, Face face)
{
    Vec3f a = scene.vertex_data[face.v0_id - 1];
    Vec3f b = scene.vertex_data[face.v1_id - 1];
    Vec3f c = scene.vertex_data[face.v2_id - 1];

    float matrixA[3][3] = {a.x - b.x, a.x - c.x, ray.dir.x,
                           a.y - b.y, a.y - c.y, ray.dir.y,
                           a.z - b.z, a.z - c.z, ray.dir.z};
    float detA = determinant(matrixA);
    if (detA == 0)
        return -1;

    // Cramers Rule
    float matrixBeta[3][3] = {a.x - ray.origin.x, a.x - c.x, ray.dir.x,
                              a.y - ray.origin.y, a.y - c.y, ray.dir.y,
                              a.z - ray.origin.z, a.z - c.z, ray.dir.z};
    float beta = determinant(matrixBeta) / detA;
    if (beta < 0)
        return -1;

    float matrixGama[3][3] = {a.x - b.x, a.x - ray.origin.x, ray.dir.x,
                              a.y - b.y, a.y - ray.origin.y, ray.dir.y,
                              a.z - b.z, a.z - ray.origin.z, ray.dir.z};
    float gama = determinant(matrixGama) / detA;
    if (gama < 0 || gama + beta > 1)
        return -1;

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

    float matrixA[3][3] = {a.x - b.x, a.x - c.x, ray.dir.x,
                           a.y - b.y, a.y - c.y, ray.dir.y,
                           a.z - b.z, a.z - c.z, ray.dir.z};
    float detA = determinant(matrixA);
    if (detA == 0)
        return -1;

    // Cramers Rule
    float matrixBeta[3][3] = {a.x - ray.origin.x, a.x - c.x, ray.dir.x,
                              a.y - ray.origin.y, a.y - c.y, ray.dir.y,
                              a.z - ray.origin.z, a.z - c.z, ray.dir.z};
    float beta = determinant(matrixBeta) / detA;
    if (beta < 0)
        return -1;

    float matrixGama[3][3] = {a.x - b.x, a.x - ray.origin.x, ray.dir.x,
                              a.y - b.y, a.y - ray.origin.y, ray.dir.y,
                              a.z - b.z, a.z - ray.origin.z, ray.dir.z};
    float gama = determinant(matrixGama) / detA;
    if (gama < 0 || gama + beta > 1)
        return -1;

    float matrixT[3][3] = {a.x - b.x, a.x - c.x, a.x - ray.origin.x,
                           a.y - b.y, a.y - c.y, a.y - ray.origin.y,
                           a.z - b.z, a.z - c.z, a.z - ray.origin.z};
    float t = determinant(matrixT) / detA;
    return t;
}
float intersectSphere(Ray r, Sphere s)
{
    float A, B, C; //constants for the quadratic equation
    float delta;
    Vec3f c = scene.vertex_data[s.center_vertex_id - 1];
    float t, t1, t2;

    C = (r.origin.x - c.x) * (r.origin.x - c.x) + (r.origin.y - c.y) * (r.origin.y - c.y) + (r.origin.z - c.z) * (r.origin.z - c.z) - s.radius * s.radius;
    B = 2 * r.dir.x * (r.origin.x - c.x) + 2 * r.dir.y * (r.origin.y - c.y) + 2 * r.dir.z * (r.origin.z - c.z);
    A = r.dir.x * r.dir.x + r.dir.y * r.dir.y + r.dir.z * r.dir.z;
    delta = B * B - 4 * A * C;

    if (delta < 0)
        return -1;
    else if (delta == 0)
    {
        t = -B / (2 * A);
    }
    else
    {
        delta = sqrt(delta);
        A = 2 * A;
        t1 = (-B + delta) / A;
        t2 = (-B - delta) / A;

        if (t1 < t2)
            t = t1;
        else
            t = t2;
    }

    return t;
}


Vec3f computeColor(Ray r)
{
    int i;
    Vec3f c;
    float minT = 90000; // some large number
    float t;
    Vec3f L, N;
    Vec3f P;
    int minI;

    c.x = c.y = c.z = 0;
    minI = -1;

    vector<Sphere> spheres = scene.spheres;
    for (i = 0; i < spheres.size(); i++)
    {
        t = intersectSphere(r, spheres[i]);
        if (t < minT && t >= 0)
        {
            // Material mat = scene.materials[spheres[i].material_id - 1];
            Vec3f tempColor = {1, 0, 0};
            c = tempColor; // can be replaced with any material property
            minI = i;
            minT = t;
        }
    }

    // triangle intersections

    vector<Triangle> triangles = scene.triangles;
    for (i = 0; i < triangles.size(); i++)
    {
        t = intersectTriangle(r, triangles[i]);
        if (t > 0 && t < minT)
        {
            // Material mat = scene.materials[spheres[i].material_id - 1];
            Vec3f tempColor = {1, 0, 0};
            c = tempColor; // can be replaced with any material property
            minI = i;
            minT = t;
        }
    }

    vector<Mesh> meshes = scene.meshes;
    for (i = 0; i < meshes.size(); i++)
    {
        vector<Face> faces = meshes[i].faces;
        for (int j = 0; j < faces.size(); j++)
        {
            t = intersectFace(r, faces[j]);
            if (t > 0 && t < minT)
            {
                // Material mat = scene.materials[spheres[i].material_id - 1];
                Vec3f tempColor = {1, 0, 0};
                c = tempColor; // can be replaced with any material property
                minI = i;
                minT = t;
            }
        }
    }
    // if (minI!=-1)
    // {
    //     P = add(r.o,multS(r.d,minT));
    //     L = add(light,multS(P,-1));
    //     N = add(P,multS(spheres[minI].center,-1));
    //     L = normalize(L);
    //     N = normalize(N);
    //     c = multS(c,dot(L,N));
    // }
    return c;
}


//L_a = k_a * l_a
Vec3f L_a(Scene scene, int i)
{
    return multVector(scene.materials[i].ambient, scene.ambient_light);
}
//L_d  = k_d * cos(the)' * E_i
//cos(the)' = max (0,w_i.n)
Vec3f L_d(Scene scene, int i, Vec3f w_i, Vec3f n, Vec3f E_i)
{
    Vec3f result;
    float cos_theta = max((double)0, dotProduct(w_i, n));
    Vec3f k_d = scene.materials[i].diffuse;
    result = multVector(multScaler(k_d, cos_theta), E_i);
    return result;
}
//L_s =k_s * cos(alp)' ^p * E_i
//cos(alp)' = max (0,n.h)
//vec3f h = unitVector(add(w_i, w_o));
Vec3f L_s(Scene scene, int i, Vec3f h, Vec3f n, Vec3f E_i)
{
    Vec3f result;
    Vec3f k_s = scene.materials[i].specular;
    float cos_alpha = max(0.0, dotProduct(n, h));
    float p = scene.materials[i].phong_exponent;
    double exp = pow(cos_alpha, p);
    result = multVector(multScaler(k_s, exp), E_i);
    return result;
}
//E_i = I/(d^2)
Vec3f E_i(Vec3f I, double d)
{
    d = 1 / (d * d);
    return multScaler(I, d);
}

//L_m = k_m * L_i(x,w_r)
Vec3f L_m()
{
    Vec3f result;
    return result;
}

Vec3f clamb(Vec3f colors){
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
int main(int argc, char *argv[])
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

    for (int camIndex = 0; camIndex < camSize; camIndex++)
    {
        Camera currCam = cameras[camIndex];
        int imageWidth = cameras[camIndex].image_height;
        int imageHeight = cameras[camIndex].image_width;

        unsigned char *image = new unsigned char[imageWidth * imageHeight * 3];

        //image produces
        for (int j = 0; j < imageHeight; j++)
        {
            for (int i = 0; i < imageWidth; i++)
            {
                //RAY TRACING
                //compute viewving ray from e to s
                Ray ray = generateRay(i, j, currCam);
                int closest_S_T_M[3] = {-1, -1, -1}; // [sphere_index,triange_index,mesh_index]
                Vec3f x;
                t_min = INF;
                double t;
                Vec3f shade;
                Vec3f pixel_color;
                pixel_color.x = scene.background_color.x;
                pixel_color.y = scene.background_color.y;
                pixel_color.z = scene.background_color.z;
                Vec3f n;        //normal vector
                Vec3f w_i, w_o; //in_vector,out_vector
                Vec3f I, E;

                //for sphere
                //for each object o:
                for (int o = 0; o < spheresSize; o++)
                {
                    //if ray intersects o at point x:
                    t = intersectSphere(ray, scene.spheres[o]);
                    // if t<tmin:
                    if (t >= 0 && t < t_min)
                    {
                        //tmin= t
                        t_min = t;
                        // obj = o;
                        closest_S_T_M[0] = o;
                        closest_S_T_M[1] = -1;
                        closest_S_T_M[2] = -1;
                        //calculate intersection point x;

                        x = add(ray.origin, multScaler(ray.dir, t_min));
                    }
                }

                // for triange

                // for meshes
                vector<Mesh> meshes = scene.meshes;
                for (int o = 0; o < meshesSize; o++)
                {
                    vector<Face> faces = meshes[o].faces;
                    for (int j = 0; j < faces.size(); j++)
                    {
                        t = intersectFace(ray, faces[j]);
                        if (t >= 0)
                        {
                            if (t < t_min)
                            {
                                //  tmin= t,
                                t_min = t;
                                //obj = o;
                                closest_S_T_M[0] = -1;
                                closest_S_T_M[1] = -1;
                                closest_S_T_M[2] = o;
                                //calculate intersection point x;
                                x = add(ray.origin, multScaler(ray.dir, t_min));
                                // Material mat = scene.materials[spheres[i].material_id - 1];

                               
                            }
                        }
                    }
                }

                //////////////////buraya kadar C ye color atayıp verebilyoruz şimdi color belirle/////

                //if object is not null:
                if (triangelsSize > 0 || spheresSize > 0 || meshesSize > 0)
                {
                    //  pixel color = La ->ambiant shading is not effected by shadosw

                    if (closest_S_T_M[0] != -1)
                    { //closest is a sphere
                        Sphere curr_sphere = scene.spheres[closest_S_T_M[0]];
                        //pixel_color = L_a(scene, curr_sphere.material_id - 1);

                        //  for easch light l:
                        for (int l = 0; l < lightSize; l++)
                        {
                            //compute shadow rays from x to l;
                            w_i = substract(x, scene.point_lights[l].position);
                            w_o = substract(x, currCam.position);
                            I = scene.point_lights[l].intensity;
                            E = E_i(I, length(w_i));
                            Vec3f center = scene.vertex_data[curr_sphere.center_vertex_id - 1];
                            n = substract(center, x);
                            n = multScaler(n, 1 / curr_sphere.radius);
                            Vec3f L_dif = L_d(scene, curr_sphere.material_id - 1,unitVector(w_i), n, E);
                            pixel_color = add(pixel_color,L_dif);
                            pixel_color = L_dif;

                            // foreach object p:
                            //     if s intersects p before the light source:
                            //     continue the light loop; // point is in shadow – no contribution from this light
                            //     pixel color += Ld + Ls // add diffuse and specular components for this light source
                        }
                    }
                    else if (closest_S_T_M[1] != -1)
                    { //closest is a triangle
                        pixel_color = L_a(scene, scene.triangles[closest_S_T_M[1]].material_id - 1);
                        for (int l = 0; l < lightSize; l++)
                        {
                            //compute shadow rays from x to l;
                            /*
                            w_i = substract(x, scene.point_lights[l].position);
                            w_o = substract(x, currCam.position);
                            I = scene.point_lights[l].intensity;
                            E = E_i(I, length(w_i));
                            int first_edge_vertex_id = scene.triangles[closest_S_T_M[1]].indices.v0_id;

                            n = substract(center, x);
                            n = multScaler(n, 1 / curr_sphere.radius);
                            Vec3f L_dif = L_d(scene, curr_sphere.material_id - 1,unitVector(w_i), n, E);
                            pixel_color = add(pixel_color,L_dif);
                            pixel_color = L_dif;
                            */
                            // foreach object p:
                            //     if s intersects p before the light source:
                            //     continue the light loop; // point is in shadow – no contribution from this light
                            //     pixel color += Ld + Ls // add diffuse and specular components for this light source
                        }

                    }
                    else if (closest_S_T_M[2] != -1)
                    { //closest is a mesh
                        pixel_color = L_a(scene, scene.meshes[closest_S_T_M[2]].material_id - 1);
                    }

                    //      compute shadow ray s from x to l:
                    //      for each object p:
                    //          if s intersects p before the light source:
                    //              continue the light loop
                    //      pixel color = pixel color + Ld +Ls ->diffuse and specular components
                }

                else
                {
                    //  pixel color = color of backround
                }

                //////////////////////////////////////

                // do not forget clamping the pixel value to [0,255] range and rounding it to the nearest integer

                // Vec3f shaded_pixel =

                // std::cout << "Ray dir: "<< ray.dir.x << ray.dir.y << ray.dir.z << " o: " << ray.origin.x << ray.origin.y << ray.origin.z << std::endl;
                //Vec3f rayColor = computeColor(ray);
                Vec3f rayColor = clamb(pixel_color);

                if (rayColor.x > 0.01 || rayColor.y > 0.01 || rayColor.z > 0.01)
                {
                    // cout << "Raycolor:" << rayColor.x << " "<< rayColor.y << " "<< rayColor.z<<endl;
                }
                int imgIndex = 3 * (i + j * currCam.image_width);
                image[imgIndex] = (unsigned char)(rayColor.x);
                image[imgIndex + 1] = (unsigned char)(rayColor.y);
                image[imgIndex + 2] = (unsigned char)(rayColor.z);
            }
        }
        write_ppm(currCam.image_name.c_str(), image, imageWidth, imageHeight);
    }
}

//compute eye rays
//ray objesi olucak
//intersection yaz
///color func
//delta fonk
//trişange intersection algo bak
//shade
//mirror
