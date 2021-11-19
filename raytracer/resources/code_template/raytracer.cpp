#include <iostream>
#include "ppm.h"
#include "helper.h"
#include <math.h>

#define INF 999999.0

typedef unsigned char RGB[3];

using namespace parser;
using namespace std;

parser::Scene scene;
Vec3f shadeSphere(int objIndex, Vec3f intersectionPoint, Vec3f camPos);
Vec3f shadeTriangle(int objIndex, Vec3f intersectionPoint, Vec3f camPos);
Vec3f shadeMesh(int objIndex, Vec3f intersectionPoint, Vec3f camPos, int TriIndex);

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
    float cos_alpha = max(0.0f, dotProduct(n, h));
    float p = mat.phong_exponent;
    float exp = pow(cos_alpha, p);
    result = multVector(multScaler(k_s, exp), E_i);
    return result;
}
//E_i = I/(d^2)
Vec3f E_i(Vec3f I, float d)
{
    return multScaler(I,  1.0 / (d * d));
}

//L_m = k_m * L_i(x,w_r)
Vec3f L_m()
{
    Vec3f result;
    return result;
}


int main(int argc, char *argv[])
{
    scene.loadFromXml(argv[1]);
    vector<Camera> cameras = scene.cameras;
    int numOfCams = cameras.size();
    int numOfTriangles = scene.triangles.size();
    int numOfMeshes = scene.meshes.size();
    int numOfSpheres = scene.spheres.size();
    int recursionDepth = scene.max_recursion_depth;
    Vec3i bgColor = scene.background_color;

    // ·NearPlane attribute defines the coordinates of the image plane with Left, Right, Bottom,
    // Top floating point parameters, respectively.
    // ·NearDistance defines the distance of the image plane to the camera.

    vector<Vec3f> triNormals;


    for (int camIndex = 0; camIndex < numOfCams; camIndex++)
    {
        Camera currCam = cameras[camIndex];
        int imageWidth = currCam.image_height;
        int imageHeight = currCam.image_width;

        unsigned char *image = new unsigned char[imageWidth * imageHeight * 3];

        float t_min;
        // Iterate over the image plane
        for (int j = 0; j < imageHeight; j++)
        {
            for (int i = 0; i < imageWidth; i++)
            {
                //RAY TRACING
                //compute viewing ray from e to s(i,j)
                Ray ray = generateRay(i, j, currCam);

                int closestObjIndex = -1;
                ObjectType closestObjType = NONE;

                t_min = INF;
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
               
                Vec3f shadedColor;
                Vec3f intersectionPoint = ray.origin + multScaler(ray.dir, t_min);
                switch (closestObjType)
                {
                    case SPHERE:
                        shadedColor = shadeSphere(closestObjIndex, intersectionPoint, currCam.position);
                        break;

                    case TRIANGLE:
                        shadedColor = shadeTriangle(closestObjIndex, intersectionPoint, currCam.position);
                        break;

                    case MESH:
                        shadedColor = shadeMesh(closestObjIndex, intersectionPoint, currCam.position,meshTriangleIndex);
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

Vec3f shadeSphere(int objIndex, Vec3f intersectionPoint, Vec3f camPos)
{
    Sphere curr_sphere = scene.spheres[objIndex];
    Material mat = scene.materials[curr_sphere.material_id -1];
    Vec3f pixel_color = L_a(mat);

    //  Shading
    //calculate intersection point x
    //     P = add(r.o,multS(r.d,minT));  -- p = x
    //     L = add(light,multS(P,-1));   -- w_i
    //     N = add(P,multS(spheres[minI].center,-1));  -- n

    //     L = normalize(L);  N = normalize(N);
    //     c = multS(c,dot(L,N));
    int numOfLights = scene.point_lights.size();
    for (int i = 0; i < numOfLights; i++)
    {
        PointLight currLight = scene.point_lights[i];
        //compute shadow rays from x to l;
        Vec3f w_i = currLight.position - intersectionPoint;
        Vec3f w_o = camPos - intersectionPoint; 

        Vec3f I = currLight.intensity;
        Vec3f E = E_i(I, length(w_i));
        Vec3f center = scene.vertex_data[curr_sphere.center_vertex_id - 1];
        Vec3f n = intersectionPoint - center;
        n = multScaler(n, 1.0 / curr_sphere.radius);

        Vec3f L_dif = L_d(mat, makeUnitVector(w_i), n, E);
        pixel_color = pixel_color + L_dif;

        float theta = acos(dotProduct(makeUnitVector(w_i),n) );
        if (theta<90) {
            Vec3f h = makeUnitVector(w_i+w_o);
        
            Vec3f L_spe = L_s(mat, h, n, E);
            pixel_color = pixel_color + L_spe;
        }
        

        
        // foreach object p:
        //     if s intersects p before the light source:
        //     continue the light loop; // point is in shadow – no contribution from this light
        //     pixel color += Ld + Ls // add diffuse and specular components for this light source
    }
    return pixel_color;
}
Vec3f shadeTriangle(int objIndex, Vec3f intersectionPoint, Vec3f camPos)
{
    Triangle currTri = scene.triangles[objIndex];
    Material mat = scene.materials[currTri.material_id -1];

    Vec3f pixel_color = L_a(mat);

    int numOfLights = scene.point_lights.size();
    for (int l = 0; l < numOfLights; l++)
    {
        //compute shadow rays from x to l;
        
        Vec3f w_i = scene.point_lights[l].position - intersectionPoint;
        Vec3f w_o = camPos - intersectionPoint;

        Vec3f I = scene.point_lights[l].intensity;
        Vec3f E = E_i(I, length(w_i));
        Triangle curr_triangle = scene.triangles[objIndex];
        Vec3f a_vertex = scene.vertex_data[curr_triangle.indices.v0_id];
        Vec3f b_vertex = scene.vertex_data[curr_triangle.indices.v1_id];
        Vec3f c_vertex = scene.vertex_data[curr_triangle.indices.v2_id];
        Vec3f ab_edge = b_vertex - a_vertex;
        Vec3f ac_edge = c_vertex - a_vertex;
        Vec3f n = cross(ab_edge,ac_edge);
        n = multScaler(n,1/length(n));

        Vec3f L_dif = L_d(mat, makeUnitVector(w_i), n, E);
        pixel_color = pixel_color + L_dif;
    }
    return pixel_color;
}
Vec3f shadeMesh(int objIndex, Vec3f intersectionPoint, Vec3f camPos,int meshTriangleIndex)
{
    Mesh currMesh = scene.meshes[objIndex];
    Material mat = scene.materials[currMesh.material_id - 1];
    Face face = currMesh.faces[meshTriangleIndex];
    Vec3f pixel_color = L_a(mat);
    int numOfLights = scene.point_lights.size();
    for (int l = 0; l < numOfLights; l++)
    {
        //compute shadow rays from x to l;
        
        Vec3f w_i = scene.point_lights[l].position - intersectionPoint;
        Vec3f w_o = camPos - intersectionPoint;

        Vec3f I = scene.point_lights[l].intensity;
        Vec3f E = E_i(I, length(w_i));
        Vec3f a_vertex = scene.vertex_data[face.v0_id];
        Vec3f b_vertex = scene.vertex_data[face.v1_id];
        Vec3f c_vertex = scene.vertex_data[face.v2_id];
        Vec3f ab_edge = b_vertex - a_vertex;
        Vec3f ac_edge = c_vertex - a_vertex;

        Vec3f n = cross(makeUnitVector(ab_edge),makeUnitVector(ac_edge));
        n = multScaler(n,1/length(n));

        Vec3f L_dif = L_d(mat, makeUnitVector(w_i), n, E);
        pixel_color = pixel_color + L_dif;
    }

    return pixel_color;
}
