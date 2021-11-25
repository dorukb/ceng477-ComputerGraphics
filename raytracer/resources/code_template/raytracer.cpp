#include <iostream>
#include "ppm.h"
#include "helper.h"
#include <math.h>
#include <thread>

#define INF 999999.0

typedef unsigned char RGB[3];

using namespace parser;
using namespace std;

struct traceArguments{
    int cameraIndex;
    int threadIndex;
    unsigned char *imagePtr;
};
parser::Scene scene;

Vec3f shadeTriangle(int objIndex,Material currMat, Vec3f intersectionPoint, Vec3f camPos, PointLight currLight,int depth,Camera currCam);
Vec3f shadeMesh(int objIndex, Material currMat, Vec3f intersectionPoint, Vec3f camPos,int meshTriangleIndex, PointLight currLight,int depth,Camera currCam);
Vec3f shadeSphere(int objIndex, Material currMat, Vec3f intersectionPoint, Vec3f camPos, PointLight currLight,int depth,Camera currCam);
bool doesIntersectWithMesh(Ray ray, float tLight);
Vec3f calculateColor(Ray ray, int recursionDepth,Camera currCam);
int numOfSpheres;
int numOfTriangles;
int numOfMeshes;
int numOfLights;
int numOfCams ;
Vec3i bgColor;
vector<Vec3f> triNormals;
vector<vector<Vec3f>> meshFaceNormals;

static int threadCount = 8;
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
    if (fabs(detA) < 0)
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
    // float delta;
    Vec3f center = scene.vertex_data[s.center_vertex_id - 1];
    float t, t1, t2;

    Vec3f oMinusC = r.origin - center;

    float c = dotProduct(oMinusC,oMinusC) - (pow(s.radius, 2));
    float b =  2 * dotProduct(r.dir, oMinusC);
    float a = dotProduct(r.dir, r.dir);
    float delta = pow(b,2) - (4 * a * c);

    if(delta < 0.0)
    {
        return -1;
    }
    else
    {
        delta = sqrtf(delta);
        a = 2.0 * a;
        float t1 = (-b + delta) / a;
        float t2 = (-b - delta) / a;
        t = t1 < t2 ? t1 : t2;
    }
    return t;
}


void calculateTriNormals()
{
    Vec3f a,b,c,ab,ac,n;
    for(int i = 0; i < numOfTriangles; i++)
    {
        a = scene.vertex_data[scene.triangles[i].indices.v0_id -1];
        b = scene.vertex_data[scene.triangles[i].indices.v1_id -1];
        c = scene.vertex_data[scene.triangles[i].indices.v2_id -1];
        ab = b - a;
        ac = c - a;

        n = makeUnitVector(cross(ab, ac));
        triNormals[i] = n;
    }
}
void calculateMeshFaceNormals()
{
    Vec3f a,b,c,ab,ac,n;
    Mesh mesh;
    for(int i = 0; i < numOfMeshes; i++)
    {   
        mesh = scene.meshes[i];
        int numOfFaces = mesh.faces.size();    
        meshFaceNormals[i].resize(numOfFaces);
        for(int j = 0; j < numOfFaces; j++)
        {
            // calculate normal for this face, store it in meshFaceNormals[i][j];
            a = scene.vertex_data[mesh.faces[j].v0_id -1];
            b = scene.vertex_data[mesh.faces[j].v1_id -1];
            c = scene.vertex_data[mesh.faces[j].v2_id -1];
            ab = b - a;
            ac = c - a;

            n = makeUnitVector(cross(ab, ac));
            meshFaceNormals[i][j] = n;
        }
    }
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
    if (!mat.is_mirror || depth == 0){
        return {0,0,0};
    }
    else {
        //float cos_theta = max(0.0f, );
        Vec3f w_r =  multScaler(multScaler(n,2.0f),dotProduct(n, w_o))-w_o; 
        w_r = makeUnitVector(w_r);

        Ray reflection_ray;
        reflection_ray.origin = intersection+ multScaler(n,scene.shadow_ray_epsilon);
        reflection_ray.dir = w_r;
        return multVector(k_m, calculateColor(reflection_ray, depth-1,currCam));
    }

}
//E_i = I/(d^2)
// Vec3f E_i(Vec3f I, float d)
// {
//     return multScaler(I,  1.0 / (d * d));
// }

Vec3f calculateColor(Ray ray, int recursionDepth,Camera currCam)
{
    int closestObjIndex = -1;
    ObjectType closestObjType = NONE;
    float t_min = INF;
    double t = -1;

    int meshTriangleIndex = -1;
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
    Sphere currSphere ;
    Vec3f w_o, n , center;
    Material currMat;

    switch (closestObjType)
    {
        case SPHERE:    
            currMat = scene.materials[scene.spheres[closestObjIndex].material_id-1];
            shadedColor = L_a(currMat);
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

                shadedColor = shadedColor + shadeSphere(closestObjIndex, currMat, intersectionPoint, ray.origin, currLight, recursionDepth, currCam);
            }
            currSphere = scene.spheres[closestObjIndex];
            w_o = ray.origin - intersectionPoint; 
            w_o = makeUnitVector(w_o);
            center = scene.vertex_data[currSphere.center_vertex_id - 1];
            n = intersectionPoint - center;
            n.x /= currSphere.radius;
            n.y /= currSphere.radius;
            n.z /= currSphere.radius;
            n = makeUnitVector(n);
            shadedColor = shadedColor + L_m(scene.materials[currSphere.material_id -1],w_o,n,recursionDepth,intersectionPoint, currCam);

            break;

        case TRIANGLE:
            currMat = scene.materials[scene.triangles[closestObjIndex].material_id-1];
            shadedColor = L_a(currMat);
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

                shadedColor = shadedColor + shadeTriangle(closestObjIndex, currMat, intersectionPoint, ray.origin, currLight,recursionDepth,currCam);
            }
            w_o = ray.origin - intersectionPoint; 
            w_o = makeUnitVector(w_o);

            // get the precalculated normal for this triangle.
            n = triNormals[closestObjIndex];
            shadedColor = shadedColor + L_m(currMat,w_o,n,recursionDepth,intersectionPoint, currCam);

            break;

        case MESH:
            currMat = scene.materials[scene.meshes[closestObjIndex].material_id-1];
            shadedColor = L_a(currMat);
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

                shadedColor = shadedColor + shadeMesh(closestObjIndex, currMat, intersectionPoint, ray.origin, meshTriangleIndex, currLight,recursionDepth,currCam);
            }
            w_o = ray.origin - intersectionPoint; 
            w_o = makeUnitVector(w_o);

            // get the precalculated normal for this face.
            n = meshFaceNormals[closestObjIndex][meshTriangleIndex];
            shadedColor = shadedColor + L_m(currMat,w_o,n,recursionDepth,intersectionPoint, currCam);

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

void raytraceMain(traceArguments args)
{
    int threadIndex = args.threadIndex;
    int cameraIndex = args.cameraIndex;

    Camera currCam = scene.cameras[cameraIndex];
    int recursionDepth = scene.max_recursion_depth;
    // calculate the starting and ending indices for width & height using these parameters.
    int width = currCam.image_width;
    int height = currCam.image_height;

    int startingHeight = threadIndex * (height / 8);
    int endingHeight = startingHeight + (height / 8);

    // unsigned char *image = ;
    // Iterate over the image plane
    for (int j = startingHeight; j <= endingHeight; j++)
    {
        for (int i = 0; i < width; i++)
        {
            //RAY TRACING
            Ray ray = generateRay(i, j, currCam);
            Vec3f shadedColor = calculateColor(ray,recursionDepth, currCam);
            Vec3i rayColor = clamp(shadedColor);

            int imgIndex = 3 * (i + j * currCam.image_width);
            args.imagePtr[imgIndex] = (unsigned char)(rayColor.x);
            args.imagePtr[imgIndex + 1] = (unsigned char)(rayColor.y);
            args.imagePtr[imgIndex + 2] = (unsigned char)(rayColor.z);
        }
    }
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


    // Precalculate mesh face normals here
    meshFaceNormals.resize(numOfMeshes);
    calculateMeshFaceNormals();

    triNormals.resize(numOfTriangles);
    calculateTriNormals();
    for (int camIndex = 0; camIndex < numOfCams; camIndex++)
    {
        Camera currCam = cameras[camIndex];
        int imageWidth = currCam.image_width;
        int imageHeight = currCam.image_height;

        unsigned char *image = new unsigned char[imageWidth * imageHeight * 3];

        vector<thread> traceThreads;
        traceThreads.resize(threadCount);

        vector<traceArguments> args;
        args.resize(threadCount);
        for(int t = 0; t < threadCount; t++)
        {
            args[t].cameraIndex = camIndex;
            args[t].threadIndex = t;
            args[t].imagePtr = image;
            traceThreads[t] = thread(raytraceMain, args[t]);
        }
        // // Iterate over the image plane
        // for (int j = 0; j < imageHeight; j++)
        // {
        //     for (int i = 0; i < imageWidth; i++)
        //     {
        //         //RAY TRACING
        //         Ray ray = generateRay(i, j, currCam);
        //         Vec3f shadedColor = calculateColor(ray,recursionDepth, currCam);
        //         Vec3i rayColor = clamp(shadedColor);

        //         int imgIndex = 3 * (i + j * currCam.image_width);
        //         image[imgIndex] = (unsigned char)(rayColor.x);
        //         image[imgIndex + 1] = (unsigned char)(rayColor.y);
        //         image[imgIndex + 2] = (unsigned char)(rayColor.z);
        //     }
        // }
         // Wait for all threads to finish before terminating main.
        // cout <<" threads created" <<  endl;
        for(int i = 0; i < threadCount; i++)
        {
            traceThreads[i].join();
        }
        // cout <<"all threads joined" << endl;
        write_ppm(currCam.image_name.c_str(), image, imageWidth, imageHeight);
    }
}

bool doesIntersectWithMesh(Ray ray, float tLight)
{
    vector<Mesh> meshes = scene.meshes;
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
Vec3f shadeSphere(int objIndex, Material currMat, Vec3f intersectionPoint, Vec3f rayOrigin, PointLight currLight,int depth,Camera currCam)
{
    Sphere curr_sphere = scene.spheres[objIndex];

    //compute shadow rays from x to l;
    Vec3f w_i = currLight.position - intersectionPoint;
    Vec3f w_o = rayOrigin - intersectionPoint; 

    Vec3f I = currLight.intensity;

    float dist = length(w_i);
    Vec3f E = multScaler(I,  1.0 / pow(dist,2));

    Vec3f center = scene.vertex_data[curr_sphere.center_vertex_id - 1];
    Vec3f n = intersectionPoint - center;

    n = makeUnitVector(n);
    w_i = makeUnitVector(w_i);
    w_o = makeUnitVector(w_o);
    Vec3f spec = {0,0,0};
    float theta = acos(dotProduct(w_i,n));
    if (theta<90) {
        Vec3f h = makeUnitVector(w_i+w_o);
        spec = L_s(currMat, h, n, E);
    }
    return spec + L_d(currMat, w_i, n, E);
}
Vec3f shadeTriangle(int objIndex, Material currMat, Vec3f intersectionPoint, Vec3f rayOrigin, PointLight currLight,int depth,Camera currCam)
{
    //compute shadow rays from x to l;   
    Vec3f w_i = currLight.position - intersectionPoint;
    Vec3f w_o = rayOrigin - intersectionPoint;

    Vec3f I = currLight.intensity;
    float dist = length(w_i);
    Vec3f E = multScaler(I,  1.0 / pow(dist,2));

    Vec3f n = triNormals[objIndex];
    w_i = makeUnitVector(w_i);
    w_o = makeUnitVector(w_o);
    float theta = acos(dotProduct(w_i,n));
    Vec3f spec = {0,0,0};
    if (theta<90) {
        Vec3f h = makeUnitVector(w_i+w_o);
        spec = L_s(currMat, h, n, E);
    }
    return spec + L_d(currMat, w_i, n, E);
}
Vec3f shadeMesh(int objIndex, Material currMat, Vec3f intersectionPoint, Vec3f rayOrigin,int meshTriangleIndex, PointLight currLight,int depth,Camera currCam )
{   
    Vec3f w_i = currLight.position - intersectionPoint;
    Vec3f w_o = rayOrigin - intersectionPoint;

    Vec3f I = currLight.intensity;
    float dist = length(w_i);
    Vec3f E = multScaler(I,  1.0 / pow(dist,2));

    Vec3f n = meshFaceNormals[objIndex][meshTriangleIndex];
    w_i = makeUnitVector(w_i);
    w_o = makeUnitVector(w_o);

    float theta = acos(dotProduct(w_i,n));
    Vec3f spec = {0,0,0};
    if (theta<90) {
        Vec3f h = makeUnitVector(w_i+w_o);
        spec = L_s(currMat, h, n, E);
    }
    return spec + L_d(currMat, w_i, n, E);
}
