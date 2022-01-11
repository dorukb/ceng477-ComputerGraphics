#ifndef _SCENE_H_
#define _SCENE_H_

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "Camera.h"
#include "Color.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "Vec4.h"
#include <unordered_map>

using namespace std;

class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;

	vector< vector<Color> > image;
	vector< Camera* > cameras;
	vector< Vec3* > vertices;
	vector< Color* > colorsOfVertices;
	vector< Scaling* > scalings;
	vector< Rotation* > rotations;
	vector< Translation* > translations;
	vector< Mesh* > meshes;

	// std::unordered_map<int,int> vertexRefCounts;
	// std::set<int> sharedVertexIndices;

	Scene(const char *xmlPath);

	void initializeImage(Camera* camera);
	void forwardRenderingPipeline(Camera* camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera* camera);
	void convertPPMToPNG(string ppmFileName, int osType);

	void calculateModelingTransformations();
	bool clippingTest(Mesh::Line *line);
	void clipAndAddToLinesList(Mesh::Line *line, Mesh *mesh, Camera *camera, Matrix4 &Mvp);

	
	void draw(int x,int y,Color* c);
    double f01(double x,double y,double x0,double x1,double y0,double y1);
    double f12(double x,double y,double x1,double x2,double y1,double y2);
    double f20(double x,double y,double x2,double x0,double y2,double y0);
    void rasterization(Mesh* object);
    void midpoint_algorithm(double x_0,double y_0,Color* c_0,double x_1,double y_1,Color* c_1);

};

#endif
