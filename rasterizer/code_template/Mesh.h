#ifndef __MESH_H__
#define __MESH_H__

#include <vector>
#include "Triangle.h"
#include <iostream>
#include "Matrix4.h"
#include <map>
#include "Vec4.h"
using namespace std;

class Mesh
{

public:
  class Line{
        public:
            Vec4 v1;
            Vec4 v2;
            Line();
            Line(Vec3 *vert1, Vec3 *vert2)
            {
                v1.x = vert1->x;
                v1.y = vert1->y;
                v1.z = vert1->z;
                v1.t = 1.0; // homogenous w coord.
                v1.colorId = vert1->colorId;

                v2.x = vert2->x;
                v2.y = vert2->y;
                v2.z = vert2->z;
                v2.t = 1.0; // homogenous w coord.
                v2.colorId = vert2->colorId;
            }
            Line(Vec4 *vert1, Vec4 *vert2){
                v1.x = vert1->x;
                v1.y = vert1->y;
                v1.z = vert1->z;
                v1.t = vert1->t; // homogenous w coord.
                v1.colorId = vert1->colorId;

                v2.x = vert2->x;
                v2.y = vert2->y;
                v2.z = vert2->z;
                v2.t = vert2->t; // homogenous w coord.
                v2.colorId = vert2->colorId;
            }
    };

    int meshId;
    int type; // 0 for wireframe, 1 for solid
    int numberOfTransformations;
    vector<int> transformationIds;
    vector<char> transformationTypes;
    int numberOfTriangles;
    vector<Triangle> triangles;


    // custom fields
    Matrix4 modelM;
    std::map<int, Vec3> sharedVerticesCopy;
    vector<Vec4> vertexDataCopy;

    // For clipping and rendering in wireframe mode
    vector<Line> lines;


    Mesh();
    Mesh(int meshId, int type, int numberOfTransformations,
          vector<int> transformationIds,
          vector<char> transformationTypes,
          int numberOfTriangles,
          vector<Triangle> triangles);

    friend ostream &operator<<(ostream &os, const Mesh &m);

  
    struct Tri{
        vector<Line> lines;
    };

   
};

#endif