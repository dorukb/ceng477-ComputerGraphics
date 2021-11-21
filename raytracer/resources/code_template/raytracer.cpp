#include <iostream>
#include "ppm.h"
#include "helper.h"
#include <math.h>

#define INF 999999.0

typedef unsigned char RGB[3];

using namespace parser;
using namespace std;

parser::Scene scene;

Vec3f shadeTriangle(int objIndex, Vec3f intersectionPoint, Vec3f camPos, PointLight currLight,int depth,Camera currCam);
Vec3f shadeMesh(int objIndex, Vec3f intersectionPoint, Vec3f camPos,int meshTriangleIndex, PointLight currLight,int depth,Camera currCam);
Vec3f shadeSphere(int objIndex, Vec3f intersectionPoint, Vec3f camPos, PointLight currLight,int depth,Camera currCam);
bool doesIntersectWithMesh(Ray ray, float tLight);
Vec3f colorOfPixel(Ray ray, int recursionDepth,Camera currCam);
int numOfSpheres;
int numOfTriangles;
int numOfMeshes;
int numOfLights;
int numOfCams ;
Vec3i bgColor;
vector<Vec3f> triNormals;

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

    m = e + multScaler(gaze, dist);
    q = m + multScaler(u, left)+ multScaler(v, top);
    s = q + multScaler(u, su) + multScaler(v, -sv);

    result.origin = e;
    result.dir = s + multScaler(e, -1);

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


enum ObjectType
{
    NONE,
    TRIANGLE,
    SPHERE,
    MESH,
};
//L_a = k_a * l_a
Vec3f L_a(Material mat)
{
    return multVector(mat.ambient, scene.ambient_light);
}
//L_d  = k_d * cos(the)' * E_i
//cos(the)' = max (0,w_i.n)
Vec3f L_d(Material mat, Vec3f w_i, Vec3f n, Vec3f E_i)
{
    float cos_theta = max(0.0f, dotProduct(w_i, n));
    Vec3f k_d = mat.diffuse;
    Vec3f result = multVector(multScaler(k_d, cos_theta), E_i);
    return result;
}
//L_s =k_s * cos(alp)' ^p * E_i
//cos(alp)' = max (0,n.h)
//vec3f h = unitVector(add(w_i, w_o));
Vec3f L_s(Material mat, Vec3f h, Vec3f n, Vec3f E_i)
{
    Vec3f result;
    Vec3f k_s = mat.specular;
    double cos_alpha = max(0.0f, dotProduct(n, h));
    float p = mat.phong_exponent;
    float exp = pow(cos_alpha, p);
    result = multVector(multScaler(k_s, exp), E_i);
    return result;
}

//L_m = k_m * L_i(x,w_r)
// eğer k_m mirror 0 sa her yerden 0 l_m olcak
    // değilse recursive 
    // yansıyan vektor m olsun m =e+dt e intersection point d de gittiği yön yani w_r yönü olucak
Vec3f L_m(Material mat,Vec3f w_o,Vec3f n,int depth,Vec3f intersection,Camera currCam)
{
    Vec3f result;
    Vec3f k_m = mat.mirror;
    if ((k_m.x <= 0.0001 && k_m.y <= 0.0001 && k_m.z <= 0.0001)|| depth == 0){
        return {0,0,0};
    }
    else {
        //float cos_theta = max(0.0f, );
        Vec3f w_r =  multScaler(multScaler(n,2),dotProduct(n, w_o))-w_o; 
        w_r = makeUnitVector(w_r);
        Ray reflection_ray ;
        reflection_ray.origin = intersection + multScaler(n,scene.shadow_ray_epsilon);
       // reflection_ray.origin = intersection;
        reflection_ray.dir = w_r;
        Vec3f L_i = colorOfPixel(reflection_ray, depth-1,currCam);
        result = multVector(k_m,L_i);
        return result;
    }

}
//E_i = I/(d^2)
Vec3f E_i(Vec3f I, float d)
{
    return multScaler(I,  1.0 / (d * d));
}

Vec3f colorOfPixel(Ray ray, int recursionDepth,Camera currCam){


                int closestObjIndex = -1;
                ObjectType closestObjType = NONE;
                float t_min = INF;
                double t;

                int meshTriangleIndex;
                // Intersect with all spheres in the scene
                for (int sphereIndex = 0; sphereIndex < numOfSpheres; sphereIndex++)
                {
                    t = intersectSphere(ray, scene.spheres[sphereIndex]);
                    if (t > 0 && t < t_min)
                    {
                        t_min = t;
                        closestObjIndex = sphereIndex;
                        closestObjType = SPHERE;
                    }
                }

                // Intersect with all triangles in the scene
                for (int triIndex = 0; triIndex < numOfTriangles; triIndex++)
                {
                    t = intersectTriangle(ray, scene.triangles[triIndex]);
                    if (t > 0 && t < t_min)
                    {
                        t_min = t;
                        closestObjIndex = triIndex;
                        closestObjType = TRIANGLE;
                    }
                }
                // Intersect with all Meshes in the scene
                vector<Mesh> meshes = scene.meshes;
                for (int meshIndex = 0; meshIndex < numOfMeshes; meshIndex++)
                {
                    vector<Face> faces = meshes[meshIndex].faces;
                    for (int j = 0; j < faces.size(); j++)
                    {
                        t = intersectFace(ray, faces[j]);
                        if (t > 0 && t < t_min)
                        {
                            t_min = t;
                            closestObjIndex = meshIndex;
                            closestObjType = MESH;
                            meshTriangleIndex = j;

                        }
                    }
                }
                // SHADING , pixel color = La -> ambient shading is not effected by shadows
               
                Vec3f shadedColor = {0,0,0};
                Vec3f intersectionPoint = ray.origin + multScaler(ray.dir, t_min);
    
                switch (closestObjType)
                {
                    case SPHERE:    
                        shadedColor = L_a(scene.materials[scene.spheres[closestObjIndex].material_id-1]);
                        for(int l = 0; l < numOfLights; l++)
                        {    
                            PointLight currLight = scene.point_lights[l];
                            Ray shadowRay;
                            shadowRay.dir = makeUnitVector(currLight.position - intersectionPoint);
                            shadowRay.origin = intersectionPoint + multScaler(shadowRay.dir, scene.shadow_ray_epsilon);
                            
                            float tLight = length(currLight.position - intersectionPoint);
                            bool inShadow = false;
                            // do intersection test with all objects, if found one before the light source, skip next calculation
                            // Intersect with all spheres in the scene
                            for (int sphereIndex = 0; sphereIndex < numOfSpheres; sphereIndex++)
                            {
                                t = intersectSphere(shadowRay, scene.spheres[sphereIndex]);
                                if (t > 0 && t < tLight)
                                {
                                    // found one obj that intersects before the light source. quit
                                    inShadow = true;
                                    break;
                                }
                            }
                            if(inShadow) continue;

                            // Intersect with all triangles in the scene
                            for (int triIndex = 0; triIndex < numOfTriangles; triIndex++)
                            {
                                t = intersectTriangle(shadowRay, scene.triangles[triIndex]);
                                if (t > 0 && t < tLight)
                                {
                                    // found one obj that intersects before the light source. quit
                                    inShadow = true;
                                    break;
                                }
                            }
                            if(inShadow) continue;
                            
                            // Intersect with all Meshes in the scene
                            if(doesIntersectWithMesh(shadowRay, tLight)){
                                inShadow = true;
                                continue;
                            }

                            shadedColor = shadedColor + shadeSphere(closestObjIndex, intersectionPoint, currCam.position, currLight, recursionDepth, currCam);
                        }
                        break;

                    case TRIANGLE:
                        shadedColor = L_a(scene.materials[scene.triangles[closestObjIndex].material_id-1]);
                        for(int l = 0; l < numOfLights; l++)
                        {    
                            PointLight currLight = scene.point_lights[l];
                            Ray shadowRay;
                            shadowRay.dir = makeUnitVector(currLight.position - intersectionPoint);
                            shadowRay.origin = intersectionPoint + multScaler(shadowRay.dir, scene.shadow_ray_epsilon);
                            
                            float tLight = length(currLight.position - intersectionPoint);
                            bool inShadow = false;
                            // do intersection test with all objects, if found one before the light source, skip next calculation
                            // Intersect with all spheres in the scene
                            for (int sphereIndex = 0; sphereIndex < numOfSpheres; sphereIndex++)
                            {
                                t = intersectSphere(shadowRay, scene.spheres[sphereIndex]);
                                if (t > 0 && t < tLight)
                                {
                                    // found one obj that intersects before the light source. quit
                                    inShadow = true;
                                    break;
                                }
                            }
                            if(inShadow) continue;

                            // Intersect with all triangles in the scene
                            for (int triIndex = 0; triIndex < numOfTriangles; triIndex++)
                            {
                                t = intersectTriangle(shadowRay, scene.triangles[triIndex]);
                                if (t > 0 && t < tLight)
                                {
                                    // found one obj that intersects before the light source. quit
                                    inShadow = true;
                                    break;
                                }
                            }
                            if(inShadow) continue;
                            
                            // Intersect with all Meshes in the scene
                            if(doesIntersectWithMesh(shadowRay, tLight)){
                                inShadow = true;
                                continue;
                            }

                            shadedColor = shadedColor + shadeTriangle(closestObjIndex, intersectionPoint, currCam.position, currLight,recursionDepth,currCam);
                        }
                        break;

                    case MESH:
                        shadedColor = L_a(scene.materials[scene.meshes[closestObjIndex].material_id-1]);
                        for(int l = 0; l < numOfLights; l++)
                        {    
                            PointLight currLight = scene.point_lights[l];
                            Ray shadowRay;
                            shadowRay.dir = makeUnitVector(currLight.position - intersectionPoint);
                            shadowRay.origin = intersectionPoint + multScaler(shadowRay.dir, scene.shadow_ray_epsilon);
                            
                            float tLight = length(currLight.position - intersectionPoint);
                            bool inShadow = false;
                            // do intersection test with all objects, if found one before the light source, skip next calculation
                            // Intersect with all spheres in the scene
                            for (int sphereIndex = 0; sphereIndex < numOfSpheres; sphereIndex++)
                            {
                                t = intersectSphere(shadowRay, scene.spheres[sphereIndex]);
                                if (t > 0 && t < tLight)
                                {
                                    // found one obj that intersects before the light source. quit
                                    inShadow = true;
                                    break;
                                }
                            }
                            if(inShadow) continue;

                            // Intersect with all triangles in the scene
                            for (int triIndex = 0; triIndex < numOfTriangles; triIndex++)
                            {
                                t = intersectTriangle(shadowRay, scene.triangles[triIndex]);
                                if (t > 0 && t < tLight)
                                {
                                    // found one obj that intersects before the light source. quit
                                    inShadow = true;
                                    break;
                                }
                            }
                            if(inShadow) continue;
                            
                            // Intersect with all Meshes in the scene
                            if(doesIntersectWithMesh(shadowRay, tLight)){
                                inShadow = true;
                                continue;
                            }

                            shadedColor = shadedColor + shadeMesh(closestObjIndex, intersectionPoint, currCam.position, meshTriangleIndex, currLight,recursionDepth,currCam);
                        }
                        break;

                    case NONE:
                        // use bg color
                        shadedColor.x = bgColor.x;
                        shadedColor.y = bgColor.y;
                        shadedColor.z = bgColor.z;
                        break;

                    default:
                        cout << "ERROR -- Closest object type is not set." << endl;
                        break;
                }

                return shadedColor;
}

int main(int argc, char *argv[])
{
    scene.loadFromXml(argv[1]);
    vector<Camera> cameras = scene.cameras;
    numOfCams = cameras.size();
    numOfTriangles = scene.triangles.size();
    numOfMeshes = scene.meshes.size();
    numOfSpheres = scene.spheres.size();
    int recursionDepth = scene.max_recursion_depth;
    numOfLights = scene.point_lights.size();
    bgColor = scene.background_color;

    // ·NearPlane attribute defines the coordinates of the image plane with Left, Right, Bottom,
    // Top floating point parameters, respectively.
    // ·NearDistance defines the distance of the image plane to the camera.

    for (int camIndex = 0; camIndex < numOfCams; camIndex++)
    {
        Camera currCam = cameras[camIndex];
        int imageWidth = currCam.image_width;
        int imageHeight = currCam.image_height;

        unsigned char *image = new unsigned char[imageWidth * imageHeight * 3];

        
        // Iterate over the image plane
        for (int j = 0; j < imageHeight; j++)
        {
            for (int i = 0; i < imageWidth; i++)
            {
                //RAY TRACING
                //compute viewing ray from e to s(i,j)
                Ray ray = generateRay(i, j, currCam);
                Vec3f shadedColor = colorOfPixel(ray,recursionDepth, currCam);

                //      compute shadow ray s from x to l:
                //      for each object p:
                //          if s intersects p before the light source:
                //              continue the light loop
                //      pixel color = pixel color + Ld +Ls ->diffuse and specular components

                // do not forget clamping the pixel value to [0,255] range and rounding it to the nearest integer
                
                Vec3i rayColor = clamp(shadedColor);

                int imgIndex = 3 * (i + j * currCam.image_width);
                image[imgIndex] = (unsigned char)(rayColor.x);
                image[imgIndex + 1] = (unsigned char)(rayColor.y);
                image[imgIndex + 2] = (unsigned char)(rayColor.z);
            }
        }
        write_ppm(currCam.image_name.c_str(), image, imageWidth, imageHeight);
    }
}

bool doesIntersectWithMesh(Ray ray, float tLight)
{
    vector<Mesh> meshes = scene.meshes;
    int numOfMeshes = scene.meshes.size();
    float t = 0;
    for (int meshIndex = 0; meshIndex < numOfMeshes; meshIndex++)
    {
        vector<Face> faces = meshes[meshIndex].faces;
        for (int j = 0; j < faces.size(); j++)
        {
            t = intersectFace(ray, faces[j]);
            if (t > 0 && t < tLight)
            {
                return true;
            }
        }
    }
    return false;
}
Vec3f shadeSphere(int objIndex, Vec3f intersectionPoint, Vec3f camPos, PointLight currLight,int depth,Camera currCam)
{
    Sphere curr_sphere = scene.spheres[objIndex];
    Material mat = scene.materials[curr_sphere.material_id -1];

    //compute shadow rays from x to l;
    Vec3f w_i = currLight.position - intersectionPoint;
    Vec3f w_o = camPos - intersectionPoint; 

    Vec3f I = currLight.intensity;
    Vec3f E = E_i(I, length(w_i));
    Vec3f center = scene.vertex_data[curr_sphere.center_vertex_id - 1];
    Vec3f n = intersectionPoint - center;
    n = multScaler(n, 1.0 / curr_sphere.radius);
    
    Vec3f res = L_d(mat, makeUnitVector(w_i), n, E);    
    
    float theta = acos(dotProduct(makeUnitVector(w_i),n) );
    if (theta<90) {
        Vec3f h = makeUnitVector(w_i+w_o);
        Vec3f L_spe = L_s(mat, h, n, E);
        res = res + L_spe;
    }
    Vec3f L_mir = L_m(mat,makeUnitVector(w_o),n,depth,intersectionPoint, currCam);
    res = res + L_mir;
    // foreach object p:
    //     if s intersects p before the light source:
    //     continue the light loop; // point is in shadow – no contribution from this light
    //     pixel color += Ld + Ls // add diffuse and specular components for this light source
    return res;
}
Vec3f shadeTriangle(int objIndex, Vec3f intersectionPoint, Vec3f camPos, PointLight currLight,int depth,Camera currCam)
{
    Triangle currTri = scene.triangles[objIndex];
    Material mat = scene.materials[currTri.material_id -1];

    //compute shadow rays from x to l;
    
    Vec3f w_i = currLight.position - intersectionPoint;
    Vec3f w_o = camPos - intersectionPoint;

    Vec3f I = currLight.intensity;
    Vec3f E = E_i(I, length(w_i));
    Vec3f a_vertex = scene.vertex_data[currTri.indices.v0_id-1];
    Vec3f b_vertex = scene.vertex_data[currTri.indices.v1_id-1];
    Vec3f c_vertex = scene.vertex_data[currTri.indices.v2_id-1];
    Vec3f ab_edge = b_vertex - a_vertex;
    Vec3f ac_edge = c_vertex - a_vertex;
    Vec3f n = cross(ab_edge,ac_edge);
    // n = multScaler(n,1/length(n));

    return L_d(mat, makeUnitVector(w_i), makeUnitVector(n), E);
}
Vec3f shadeMesh(int objIndex, Vec3f intersectionPoint, Vec3f camPos,int meshTriangleIndex, PointLight currLight,int depth,Camera currCam )
{
    Mesh currMesh = scene.meshes[objIndex];
    Material mat = scene.materials[currMesh.material_id - 1];
    Face face = currMesh.faces[meshTriangleIndex];
    
    Vec3f w_i = currLight.position - intersectionPoint;
    Vec3f w_o = camPos - intersectionPoint;

    Vec3f I = currLight.intensity;
    Vec3f E = E_i(I, length(w_i));
    Vec3f a_vertex = scene.vertex_data[face.v0_id -1];
    Vec3f b_vertex = scene.vertex_data[face.v1_id -1];
    Vec3f c_vertex = scene.vertex_data[face.v2_id -1];
    Vec3f ab_edge = b_vertex - a_vertex;
    Vec3f ac_edge = c_vertex - a_vertex;

    Vec3f n = cross(ab_edge, ac_edge);
    n = multScaler(n,1.0 / length(n));

    return L_d(mat, makeUnitVector(w_i), n, E) + L_m(mat, makeUnitVector(w_o), n, depth, intersectionPoint, currCam);
}
