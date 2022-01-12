#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>

#include "Scene.h"
#include "Camera.h"
#include "Color.h"
#include "Mesh.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "tinyxml2.h"
#include "Helpers.h"
#include "Matrix4.h"

using namespace tinyxml2;
using namespace std;

/*
	Transformations, clipping, culling, rasterization are done here.
	You may define helper functions.
*/

void Scene::draw(int x, int y, Color *c, Camera *cam)
{

	if (x < 0 || y < 0)
	{
		return;
	}
	if (x >= cam->horRes || y >= cam->verRes)
		return;
	/*image[x][y].r = 0;
	image[x][y].b = 0;
	image[x][y].g = 0;*/
	if(image[x][y].r==backgroundColor.r && image[x][y].g==backgroundColor.g && image[x][y].b==backgroundColor.b ){
        image[x][y].r=round(c->r);
        image[x][y].b=round(c->b);
        image[x][y].g =round(c->g);
	}

}

double Scene ::f01(double x, double y, double x0, double x1, double y0, double y1)
{
	return x * (y0 - y1) + y * (x1 - x0) + x0 * y1 - y0 * x1;
}
double Scene ::f12(double x, double y, double x1, double x2, double y1, double y2)
{
	return x * (y1 - y2) + y * (x2 - x1) + x1 * y2 - y1 * x2;
}
double Scene ::f20(double x, double y, double x2, double x0, double y2, double y0)
{
	return x * (y2 - y0) + y * (x0 - x2) + x2 * y0 - y2 * x0;
}

void Scene::midpoint_algorithm(double x_0, double y_0, Color c_0, double x_1, double y_1, Color c_1, Camera *camera)
{
	double m = (y_1 - y_0) / (x_1 - x_0);
	double x, y, d;
	double alpha = ABS(x_1 - x_0);

	int minx = min(x_0, x_1);
	int miny = min(y_0, y_1);

	int maxx = max(x_0, x_1);
	int maxy = max(y_0, y_1);
	if (minx < 0 || miny < 0 || maxx >= camera->horRes || maxy >= camera->verRes)
	{
		return;
	}
	if(0<m && m<=1){
	    if(x_1>=x_0){
            y=y_0;
            d = 2*(y_0-y_1)+(x_1-x_0);
            Color* c = &c_0;
            Color* dc = divColor(subColor(&c_1,&c_0),(x_1-x_0)) ;
            cout<<"r:"<<dc->r<<"g:"<<dc->g<<"b:"<<dc->b<<endl;
            for(x=x_0; x<x_1;x++){
                draw(x,y,c, camera);
                if(d<0){
                    y=y+1;
                    d+=2*((y_0-y_1)+(x_1-x_0));
                }else{
                    d+=2*(y_0-y_1);
                }
                c=addColor(c,dc);
            }
	    }
	    else{
            y=y_1;
            d = 2*(y_1-y_0)+(x_0-x_1);
            Color* c = &c_1;
            Color* dc = divColor(subColor(&c_0,&c_1),(x_0-x_1)) ;
            cout<<"infr:"<<dc->r<<"g:"<<dc->g<<"b:"<<dc->b<<endl;
            for(x=x_1; x<x_0;x++){
                draw(x,y,c, camera);
                if(d<0){
                    y=y+1;
                    d+=2*((y_1-y_0)+(x_0-x_1));
                }else{
                    d+=2*(y_1-y_0);
                }
                c=addColor(c,dc);
            }

        }

	}
	else if(m>1){
	    if(y_1>=y_0){
            cout<<m<<endl;
            //change x and y
            double x=x_0;
            double d = 2*(x_1-x_0)+(y_1-y_0);
            Color* c = &c_0;
            Color* dc = divColor(subColor(&c_1,&c_0),(y_1-y_0)) ;
            for(int y = y_0; y<y_1;y++){
                draw(x,y,c, camera);
                if(d<0){
                    x=x+1;
                    if(x >= camera->horRes){
                        break;
                    }
                    d+=2*((x_0-x_1)+(y_1-y_0));
                }else{
                    d+=2*(x_0-x_1);
                }
                c=addColor(c,dc);
            }
	    }
	    else{
            cout<<m<<endl;
            //change x and y
            double x=x_1;
            double d = 2*(x_0-x_1)+(y_0-y_1);
            Color* c = &c_1;
            Color* dc = divColor(subColor(&c_0,&c_1),(y_0-y_1)) ;
            for(int y = y_1; y<y_0;y++){
                draw(x,y,c, camera);
                if(d<0){
                    x=x+1;
                    if(x >= camera->horRes){
                        break;
                    }
                    d+=2*((x_1-x_0)+(y_0-y_1));
                }else{
                    d+=2*(x_1-x_0);
                }
                c=addColor(c,dc);
            }

        }

	}
	else if(-1<=m && m<0){
        if(x_1>=x_0){
            y=y_0;
            d = 2*(y_0-y_1)-(x_1-x_0);
            Color* c = &c_0;
            Color* dc = divColor(subColor(&c_1,&c_0),(x_1-x_0)) ;
            for(x=x_0; x<x_1;x++){
                draw(x,y,c, camera);
                if(d<0){
                    d+=2*(y_0-y_1);
                }else{

                    y=y-1;
                    d+=2*((y_0-y_1)-(x_1-x_0));
                }
                c=addColor(c,dc);
            }
        }
        else{
            y=y_1;
            d = 2*(y_1-y_0)-(x_0-x_1);
            Color* c = &c_1;
            Color* dc = divColor(subColor(&c_0,&c_1),(x_0-x_1)) ;
            for(x=x_1; x<x_0;x++){
                draw(x,y,c, camera);
                if(d<0){
                    d+=2*(y_1-y_0);
                }else{

                    y=y-1;
                    d+=2*((y_1-y_0)-(x_0-x_1));
                }
                c=addColor(c,dc);
            }

        }
    }
	else if (-1.0 > m)
	{
        if(y_1>=y_0){
            x=x_0;
            d = 2*(x_0-x_1)-(y_1-y_0);
            Color* c = &c_0;
            Color* dc = divColor(subColor(&c_1,&c_0),(y_1-y_0)) ;
            for(y=y_0; y<y_1;y++){
                draw(x,y,c, camera);
                if(d<0){
                    d+=2*(x_0-x_1);
                }else{

                    x=x-1;
                    d+=2*((x_0-x_1)-(y_1-y_0));
                }
                c=addColor(c,dc);
            }
        }

        else{
            x=x_1;
            d = 2*(x_1-x_0)-(y_0-y_1);
            Color* c = &c_1;
            Color* dc = divColor(subColor(&c_0,&c_1),(y_0-y_1)) ;
            for(y=y_1; y<y_0;y++){
                draw(x,y,c, camera);
                if(d<0){
                    d+=2*(x_1-x_0);
                }else{

                    x=x-1;
                    d+=2*((x_1-x_0)-(y_0-y_1));
                }
                c=addColor(c,dc);
            }

        }
	}
}

void Scene::rasterization(Mesh *object, Camera *camera)
{
	//solid
	if (object->type == 1)
	{
		double alpha, beta, gama;
		Color *c;
		Color *c_0;
		Color *c_1;
		Color *c_2;

		int numOfVert = object->vertexDataCopy.size();

		for (int i = 0; i < numOfVert-2; i += 3)
		{
			//
			// int v0_id = object->triangles[i].getFirstVertexId() - 1;
			// int v1_id = object->triangles[i].getSecondVertexId() - 1;
			// int v2_id = object->triangles[i].getThirdVertexId() - 1;


			int c0_id = (object->vertexDataCopy[i].colorId) - 1;
			int c1_id = (object->vertexDataCopy[i+1].colorId) - 1;
			int c2_id = (object->vertexDataCopy[i+2].colorId) - 1;

			double x_0 = object->vertexDataCopy[i].x;
			double y_0 = object->vertexDataCopy[i].y;
			c_0 = colorsOfVertices[c0_id];

			double x_1 = object->vertexDataCopy[i+1].x;
			double y_1 = object->vertexDataCopy[i+1].y;
			c_1 = colorsOfVertices[c1_id];

			double x_2 = object->vertexDataCopy[i+2].x;
			double y_2 = object->vertexDataCopy[i+2].y;
			c_2 = colorsOfVertices[c2_id];

			double y_min = min(y_0, min(y_1, y_2));
			double x_min = min(x_0, min(x_1, x_2));
			double y_max = max(y_0, max(y_1, y_2));
			double x_max = max(x_0, max(x_1, x_2));


			for (int y = y_min; y < y_max; y++)
			{
				for (int x = x_min; x < x_max; x++)
				{
					alpha = f12(x, y, x_1, x_2, y_1, y_2) / f12(x_0, y_0, x_1, x_2, y_1, y_2);
					beta = f20(x, y, x_2, x_0, y_2, y_0) / f20(x_1, y_1, x_2, x_0, y_2, y_0);
					gama = f01(x, y, x_0, x_1, y_0, y_1) / f01(x_2, y_2, x_0, x_1, y_0, y_1);

					if (alpha >= 0 && beta >= 0 && gama >= 0)
					{
						c = addColor(addColor(multColor(c_0,alpha),multColor(c_1,beta)),multColor(c_2,gama));
						draw(x,y,c, camera);
					}
				}
			}
		}
	}
	else
	{
		double alpha, beta, gama;
		Color *c;
		Color *c_0;
		Color *c_1;
		Color *c_2;

		int lineCount = object->lines.size();

		for (int i = 0; i < lineCount; i++)
		{
			Mesh::Line line = object->lines[i];
			double x_0 = line.v1.x;
			double y_0 = line.v1.y;
			int id0 = line.v1.colorId;
			// double x_0 = object->vertexDataCopy[3*i].x;
			// double y_0 = object->vertexDataCopy[i].y;

			// int id0 = object->vertexDataCopy[3*i+2].colorId;
			c_0 = colorsOfVertices[id0 - 1];

			double x_1 = line.v2.x;
			double y_1 = line.v2.y;
			int id1 = line.v2.colorId;
			c_1 = colorsOfVertices[id1 - 1];

			// cout<<"x_0:"<<x_0<<",y_0:"<<y_1<<"x_1:"<<x_1<<",y_1:"<<y_1<<endl;

			//draw for 0->1 1->2 2->0
			midpoint_algorithm(x_0, y_0, *c_0, x_1, y_1, *c_1,camera);

		}
	}

	//cout <<" rasterization doen for mesh: " << object->meshId << endl;
}
void Scene::forwardRenderingPipeline(Camera *camera)
{
	// TODO: Implement this function.

	// After this call, all Meshes have their Model matrix field(.modelM) assigned.

	this->calculateModelingTransformations();

	// Viewing Transformations
	Matrix4 Mcam;
	Mcam = Mcam.GetMcam(camera);

	// 1 for perspective, 0 for orthographic
	Matrix4 Mproj;
	if (camera->projectionType == 0)
	{
		Mproj = Mproj.GetMortho(camera);
	}
	else if (camera->projectionType == 1)
	{
		Mproj = Mproj.GetMpers(camera);
	}
	Matrix4 Mvp;
	Mvp = Mvp.GetMvp(camera);
	// remember to do Perspective divide if mode is perspective!

	// Do: Trans'edVertex = Mvp * PersDiv * Mper * Mcam * Mmodel * Vertex
	// for each vertex of a mesh. remember, Mmodel is specific to a Mesh, while others are generic.

	int meshCount = meshes.size();
	for (int i = 0; i < meshCount; i++)
	{
		Mesh *mesh = meshes[i];

		Matrix4 McamMmodel = multiplyMatrixWithMatrix(Mcam, mesh->modelM);
		Matrix4 MprojMcamMmodel = multiplyMatrixWithMatrix(Mproj, McamMmodel);
		// only M_vp mulp is left, do that after perspective divide if cam is perspective.

		Matrix4 MprojMcam = multiplyMatrixWithMatrix(Mproj, Mcam);
		for (auto tri : mesh->triangles)
		{
			if (mesh->type == 0) // Wireframe mode
			{
				Vec3 *v1data = vertices[tri.getFirstVertexId() - 1];
				Vec3 *v2data = vertices[tri.getSecondVertexId() - 1];
				Vec3 *v3data = vertices[tri.getThirdVertexId() - 1];

				// do also BFC here??

				// Adds w component, convert to homogenous coords.
				Vec4 v1(v1data);
				Vec4 v2(v2data);
				Vec4 v3(v3data);

				// // apply the transformation
				// v1 = multiplyMatrixWithVec4(MprojMcamMmodel, v1);
				// v2 = multiplyMatrixWithVec4(MprojMcamMmodel, v2);
				// v3 = multiplyMatrixWithVec4(MprojMcamMmodel, v3);
				// Apply modelling trasnform

				v1 = multiplyMatrixWithVec4(mesh->modelM, v1);
				v2 = multiplyMatrixWithVec4(mesh->modelM, v2);
				v3 = multiplyMatrixWithVec4(mesh->modelM, v3);

				if(this->cullingEnabled){
					// Backface culling

					// calculate normal
					Vec3 ab, ac, n;
					
					Vec3 a,b,c;
					a.x = v1.x;
					a.y = v1.y;
					a.z = v1.z;

					b.x = v2.x;
					b.y = v2.y;
					b.z = v2.z;
					c.x = v3.x;
					c.y = v3.y;
					c.z =v3.z;

					ab = subtractVec3(b, a);
					ac = subtractVec3(c, a);
					n = crossProductVec3(ab, ac);
					n = normalizeVec3(n);

					Vec3 midpoint;
					midpoint.x = (v1.x + v2.x + v3.x)/3.0;
					midpoint.y = (v1.y + v2.y + v3.y)/3.0;
					midpoint.z = (v1.z + v2.z + v3.z)/3.0;

					// Vec3 v1test;
					// v1test.x = v1.x;
					// v1test.y = v1.y;
					// v1test.z = v1.z;

					Vec3 camToMidpoint = subtractVec3(midpoint, camera->pos);
					if(dotProductVec3(n, camToMidpoint) >= 0){
						// cull this;
						continue;
                	}
				}

				v1 = multiplyMatrixWithVec4(MprojMcam, v1);
				v2 = multiplyMatrixWithVec4(MprojMcam, v2);
				v3 = multiplyMatrixWithVec4(MprojMcam, v3);

				Mesh::Line line1(&v1, &v2); // formed by: v1 and v2
				Mesh::Line line2(&v1, &v3); // formed by: v1 and v3
				Mesh::Line line3(&v2, &v3); // formed by: v2 and v3

				//do clipping test for each line, then add to lines list of our mesh.
				clipAndAddToLinesList(&line1, mesh, camera, Mvp);
				clipAndAddToLinesList(&line2, mesh, camera, Mvp);
				clipAndAddToLinesList(&line3, mesh, camera, Mvp);
				//cout <<"[CLIPPING DONE]"<< endl;

				//cout << "line count: " << mesh->lines.size() << endl;

			}

			else // solid mode
			{
				// perform BFC first?

				// transform all vertices, no need to form lines. keep the same triangle structure.
				// just copy vertex data.
				Vec3 *v1data = vertices[tri.getFirstVertexId() - 1];
				Vec3 *v2data = vertices[tri.getSecondVertexId() - 1];
				Vec3 *v3data = vertices[tri.getThirdVertexId() - 1];

				Vec4 v1(v1data);
				Vec4 v2(v2data);
				Vec4 v3(v3data);

				v1 = multiplyMatrixWithVec4(mesh->modelM, v1);
				v2 = multiplyMatrixWithVec4(mesh->modelM, v2);
				v3 = multiplyMatrixWithVec4(mesh->modelM, v3);

				if(this->cullingEnabled){
					// Backface culling

					// calculate normal
					Vec3 ab, ac, n;
					
					Vec3 a,b,c;
					a.x = v1.x;
					a.y = v1.y;
					a.z = v1.z;

					b.x = v2.x;
					b.y = v2.y;
					b.z = v2.z;
					c.x = v3.x;
					c.y = v3.y;
					c.z =v3.z;

					ab = subtractVec3(b, a);
					ac = subtractVec3(c, a);
					n = crossProductVec3(ab, ac);
					n = normalizeVec3(n);

					Vec3 midpoint;
					midpoint.x = (v1.x + v2.x + v3.x)/3.0;
					midpoint.y = (v1.y + v2.y + v3.y)/3.0;
					midpoint.z = (v1.z + v2.z + v3.z)/3.0;

					// Vec3 v1test;
					// v1test.x = v1.x;
					// v1test.y = v1.y;
					// v1test.z = v1.z;

					Vec3 camToMidpoint = subtractVec3(midpoint, camera->pos);
					if(dotProductVec3(n, camToMidpoint) >= 0){
						// cull this;
						continue;
                	}
				}
				v1 = multiplyMatrixWithVec4(MprojMcam, v1);
				v2 = multiplyMatrixWithVec4(MprojMcam, v2);
				v3 = multiplyMatrixWithVec4(MprojMcam, v3);
				
				if (camera->projectionType == 1)
				{
					// perspective cam, do perspective divide
					v1.x /= v1.t;
					v1.y /= v1.t;
					v1.z /= v1.t;
					v1.t = 1.0;

					v2.x /= v2.t;
					v2.y /= v2.t;
					v2.z /= v2.t;
					v2.t = 1.0;

					v3.x /= v3.t;
					v3.y /= v3.t;
					v3.z /= v3.t;
					v3.t = 1.0;
				}

				// complete transformations by multiplying with M_viewport
				v1 = multiplyMatrixWithVec4(Mvp, v1);
				v2 = multiplyMatrixWithVec4(Mvp, v2);
				v3 = multiplyMatrixWithVec4(Mvp, v3);
					// given that we traverse the triangles list in the same order, we can reliably get correct vertex data by indexing without vertex ids.
				mesh->vertexDataCopy.push_back(v1);
				mesh->vertexDataCopy.push_back(v2);
				mesh->vertexDataCopy.push_back(v3);

				// for (int j = 0; j < 3; j++)
				// {
				// 	int ind = tri.vertexIds[j];
				// 	Vec3 *data = vertices[ind - 1];
				// 	Vec4 vert;
				// 	vert.x = data->x;
				// 	vert.y = data->y;
				// 	vert.z = data->z;
				// 	vert.t = 1.0; // homogenous w coord.
				// 	vert.colorId = data->colorId;

				// 	Vec4 transedVert = multiplyMatrixWithVec4(MprojMcamMmodel, vert);

				// 	if (camera->projectionType == 1)
				// 	{
				// 		// perspective cam, do perspective divide
				// 		transedVert.x /= transedVert.t;
				// 		transedVert.y /= transedVert.t;
				// 		transedVert.z /= transedVert.t;
				// 		transedVert.t = 1.0;
				// 	}

				// 	// complete transformations by multiplying with M_viewport
				// 	Vec4 finalVert = multiplyMatrixWithVec4(Mvp, transedVert);

				// 	// given that we traverse the triangles list in the same order, we can reliably get correct vertex data by indexing without vertex ids.
				// 	mesh->vertexDataCopy.push_back(finalVert);

				// 	// but we will drop entire triangles, vertices during clipping? the list will get broken...
				// 	// need to delete the copies of deleted vertices aswell. vector will "move them down" to correct positions after erase...
				// 	// vec.erase( vec.begin() + 3 ); exp code to use for deletion.
				// }
			}



			// do BFC
			// add to list if visible.
		}

        rasterization(&mesh[i],camera);
		if(mesh->type==0){
			mesh->lines.clear();
		    // while(mesh->lines.size()>0){
		    //     mesh->lines.pop_back();
		    // }

		}else{
			mesh->vertexDataCopy.clear();
		    // while(mesh->vertexDataCopy.size()>0){
		    //     mesh->lines.pop_back();
		    // }
		}
	}


}

void Scene::calculateModelingTransformations()
{
	for (int i = 0; i < meshes.size(); i++)
	{
		Mesh *m = meshes[i];
		m->modelM = getIdentityMatrix();

		for (int j = 0; j < m->numberOfTransformations; j++)
		{
			// ids start from 1, again :/
			int tIndex = m->transformationIds[j] - 1;

			char type = m->transformationTypes[j];
			//cout << "index of transformation: " << tIndex << " type: " << type << endl;

			if (type == 't')
			{
				// do translation transformation
				Translation *trans = this->translations[tIndex];
				Matrix4 transMatrix(trans);

				// Accumulate the transformation.
				m->modelM = multiplyMatrixWithMatrix(transMatrix, m->modelM);
			}
			else if (type == 's')
			{
				// do scaling transformation
				Scaling *scaling = this->scalings[tIndex];
				Matrix4 scalingMatrix(scaling);

				// Accumulate the transformation.
				m->modelM = multiplyMatrixWithMatrix(scalingMatrix, m->modelM);
			}
			else if (type == 'r')
			{
				// do rotation transformation
				Rotation *rot = this->rotations[tIndex];
				Matrix4 rotationMatrix;
				rotationMatrix = rotationMatrix.GetRotationMatrix(rot);

				// Accumulate the transformation.
				m->modelM = multiplyMatrixWithMatrix(rotationMatrix, m->modelM);
			}
		}
	}
	// at this point, each mesh has the final modelling matrix info in their "mesh->modelM" field.
}

bool isVisible(double den, double num, double *te, double *tl)
{
	double t;
	if (den > 0)
	{ // potentially entering
		t = num / den;
		if (t > *tl)
		{
			return false;
		}
		if (t > *te)
		{
			*te = t;
		}
	}
	else if (den < 0)
	{ // potentially leaving
		t = num / den;
		if (t < *te)
		{
			return false;
		}
		if (t < *tl)
		{
			*tl = t;
		}
	}
	else if (num > 0)
	{ // line parallel to edge
		return false;
	}
	return true;
}

void Scene::clipAndAddToLinesList(Mesh::Line *line, Mesh *mesh, Camera *camera, Matrix4 &Mvp)
{
	if (camera->projectionType == 1)
	{
		// perspective cam, do perspective divide
		line->v1.x /= line->v1.t;
		line->v1.y /= line->v1.t;
		line->v1.z /= line->v1.t;
		line->v1.t = 1.0;

		line->v2.x /= line->v2.t;
		line->v2.y /= line->v2.t;
		line->v2.z /= line->v2.t;
		line->v2.t = 1.0;
	}
	line->v1 = multiplyMatrixWithVec4(Mvp, line->v1);
	line->v2 = multiplyMatrixWithVec4(Mvp, line->v2);

	//cout << "DO cllipping test for line: v1: " << line->v1 << "v2: "<< line->v2 << endl;
	bool isvisible = clippingTest(line, camera);
	if (isvisible)
	{
		// if (camera->projectionType == 1)
		// {
		// 	// perspective cam, do perspective divide
		// 	line->v1.x /= line->v1.t;
		// 	line->v1.y /= line->v1.t;
		// 	line->v1.z /= line->v1.t;
		// 	line->v1.t = 1.0;

		// 	line->v2.x /= line->v2.t;
		// 	line->v2.y /= line->v2.t;
		// 	line->v2.z /= line->v2.t;
		// 	line->v2.t = 1.0;
		// }
		// complete transformations by multiplying with M_viewport
		// cout << "Before VP: " << line->v1 << endl;
		// line->v1 = multiplyMatrixWithVec4(Mvp, line->v1);
		// line->v2 = multiplyMatrixWithVec4(Mvp, line->v2);
		// cout << "After VP: " << line->v1 << endl;
		//cout << "Visible line: v1: " << line->v1 << "v2: "<< line->v2 << endl;
		mesh->lines.push_back(*line);
	}
	else{
		//cout << "Clipped line: v1: " << line->v1 << "v2: "<< line->v2 << endl;

	}
}
bool Scene::clippingTest(Mesh::Line *line, Camera *cam)
{
	double dx, dy, dz;
	double xmin, xmax, ymin, ymax, zmin, zmax;

	dx = line->v2.x - line->v1.x;
	dy = line->v2.y - line->v1.y;
	dz = line->v2.z - line->v1.z;


	// dx = (line->v1.x - line->v2.x);
	// dy = (line->v1.y - line->v2.y);
	// dz = (line->v1.z - line->v2.z);
	// if (ABS(line->v1.t) > ABS(line->v2.t))
	// {
	// 	xmax = ABS(line->v1.t);
	// }
	// else
	// {
	// 	xmax = ABS(line->v2.t);
	// }
	// ymax = zmax = xmax;

	// if (-ABS(line->v1.t) < -ABS(line->v2.t))
	// {
	// 	xmin = -ABS(line->v1.t);
	// }
	// else
	// {
	// 	xmin = -ABS(line->v2.t);
	// }
	// ymin = zmin = xmin;

	ymin = zmin = xmin = 0;
	ymax = zmax = xmax = cam->horRes;
	ymax = cam->verRes;


	double* te = new double(0);
	double* tl = new double(1);

	*te = 0;
	*tl = 1;
	// double te = 0;
	// double tl = 1;
	bool visible = false;

	double x0 = line->v1.x;
	double y0 = line->v1.y;
	double z0 = line->v1.z;

	// double x0 = line->v2.x;
	// double y0 = line->v2.y;
	// double z0 = line->v2.z;
	Vec4 newV1(line->v1);
	Vec4 newV2(line->v2);
	if (isVisible(dx, xmin - x0, te, tl))
	{
		if (isVisible(-dx, x0 - xmax, te, tl))
		{
			if (isVisible(dy, ymin - y0, te, tl))
			{
				if (isVisible(-dy, y0 - ymax, te, tl))
				{
					if (isVisible(dz, zmin - z0, te, tl))
					{
						if (isVisible(-dz, z0 - zmax, te, tl))
						{
							visible = true;
							// get clipped vertex data
							if (*tl < 1)
							{
								newV2.x = x0 + dx * (*tl);
								newV2.y = y0 + dy * (*tl);
								newV2.z = z0 + dz * (*tl);
							}
							if (*te > 0)
							{
								newV1.x = x0 + dx * (*te);
								newV1.y = y0 + dy * (*te);
								newV1.z = z0 + dz * (*te);
							}
							// update the line.
							line->v1 = newV1;
							line->v2 = newV2;

							//cout <<"[VISIBLE, new vertices] " << line->v1 << " v2: " << line->v2 << endl;
						}
					}
				}
			}
		}
	}

	return visible;
}

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *pElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *pRoot = xmlDoc.FirstChild();

	// read background color
	pElement = pRoot->FirstChildElement("BackgroundColor");
	str = pElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	pElement = pRoot->FirstChildElement("Culling");
	if (pElement != NULL)
	{
		str = pElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			cullingEnabled = true;
		}
		else
		{
			cullingEnabled = false;
		}
	}

	// read cameras
	pElement = pRoot->FirstChildElement("Cameras");
	XMLElement *pCamera = pElement->FirstChildElement("Camera");
	XMLElement *camElement;
	while (pCamera != NULL)
	{
		Camera *cam = new Camera();

		pCamera->QueryIntAttribute("id", &cam->cameraId);

		// read projection type
		str = pCamera->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			cam->projectionType = 0;
		}
		else
		{
			cam->projectionType = 1;
		}

		camElement = pCamera->FirstChildElement("Position");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->pos.x, &cam->pos.y, &cam->pos.z);

		camElement = pCamera->FirstChildElement("Gaze");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->gaze.x, &cam->gaze.y, &cam->gaze.z);

		camElement = pCamera->FirstChildElement("Up");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->v.x, &cam->v.y, &cam->v.z);

		cam->gaze = normalizeVec3(cam->gaze);
		cam->u = crossProductVec3(cam->gaze, cam->v);
		cam->u = normalizeVec3(cam->u);

		cam->w = inverseVec3(cam->gaze);
		cam->v = crossProductVec3(cam->u, cam->gaze);
		cam->v = normalizeVec3(cam->v);

		camElement = pCamera->FirstChildElement("ImagePlane");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &cam->left, &cam->right, &cam->bottom, &cam->top,
			   &cam->near, &cam->far, &cam->horRes, &cam->verRes);

		camElement = pCamera->FirstChildElement("OutputName");
		str = camElement->GetText();
		cam->outputFileName = string(str);

		cameras.push_back(cam);

		pCamera = pCamera->NextSiblingElement("Camera");
	}

	// read vertices
	pElement = pRoot->FirstChildElement("Vertices");
	XMLElement *pVertex = pElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (pVertex != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = pVertex->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = pVertex->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		vertices.push_back(vertex);
		colorsOfVertices.push_back(color);

		pVertex = pVertex->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	pElement = pRoot->FirstChildElement("Translations");
	XMLElement *pTranslation = pElement->FirstChildElement("Translation");
	while (pTranslation != NULL)
	{
		Translation *translation = new Translation();

		pTranslation->QueryIntAttribute("id", &translation->translationId);

		str = pTranslation->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		translations.push_back(translation);

		pTranslation = pTranslation->NextSiblingElement("Translation");
	}

	// read scalings
	pElement = pRoot->FirstChildElement("Scalings");
	XMLElement *pScaling = pElement->FirstChildElement("Scaling");
	while (pScaling != NULL)
	{
		Scaling *scaling = new Scaling();

		pScaling->QueryIntAttribute("id", &scaling->scalingId);
		str = pScaling->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		scalings.push_back(scaling);

		pScaling = pScaling->NextSiblingElement("Scaling");
	}

	// read rotations
	pElement = pRoot->FirstChildElement("Rotations");
	XMLElement *pRotation = pElement->FirstChildElement("Rotation");
	while (pRotation != NULL)
	{
		Rotation *rotation = new Rotation();

		pRotation->QueryIntAttribute("id", &rotation->rotationId);
		str = pRotation->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		rotations.push_back(rotation);

		pRotation = pRotation->NextSiblingElement("Rotation");
	}

	// read meshes
	pElement = pRoot->FirstChildElement("Meshes");

	XMLElement *pMesh = pElement->FirstChildElement("Mesh");
	XMLElement *meshElement;
	while (pMesh != NULL)
	{
		Mesh *mesh = new Mesh();

		pMesh->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = pMesh->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = 0;
		}
		else
		{
			mesh->type = 1;
		}

		// read mesh transformations
		XMLElement *pTransformations = pMesh->FirstChildElement("Transformations");
		XMLElement *pTransformation = pTransformations->FirstChildElement("Transformation");

		while (pTransformation != NULL)
		{
			char transformationType;
			int transformationId;

			str = pTransformation->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			pTransformation = pTransformation->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *clone_str;
		int v1, v2, v3;
		XMLElement *pFaces = pMesh->FirstChildElement("Faces");
		str = pFaces->GetText();
		clone_str = strdup(str);

		row = strtok(clone_str, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			// vertices[v1-1]->refCount++;
			// if(vertices[v1-1]->refCount>1){
			// 	// shared vertice? what if its us? using this in more than one triangle??? FUCK
			// }
			// // increase ref count for v1,v2,v3. if ref count > 1, then add to shared vertices.
			// auto it = this->vertexRefCounts.find(v1);
			// if(it != vertexRefCounts.end()){
			// 	// found it, so shared.
			// 	this->sharedVertexIndices.push_back()
			// }
			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		meshes.push_back(mesh);

		pMesh = pMesh->NextSiblingElement("Mesh");
	}
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
			}

			this->image.push_back(rowOfColors);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				this->image[i][j].r = this->backgroundColor.r;
				this->image[i][j].g = this->backgroundColor.g;
				this->image[i][j].b = this->backgroundColor.b;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFileName.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFileName << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
	os_type == 1 		-> Ubuntu
	os_type == 2 		-> Windows
	os_type == other	-> No conversion
*/
void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
	string command;

	// call command on Ubuntu
	if (osType == 1)
	{
		command = "convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// call command on Windows
	else if (osType == 2)
	{
		command = "magick convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// default action - don't do conversion
	else
	{
	}
}