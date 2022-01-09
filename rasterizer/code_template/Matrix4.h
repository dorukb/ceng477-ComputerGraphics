#ifndef __MATRIX4_H__
#define __MATRIX4_H__

#include <iostream>
#include "Translation.h"
#include "Scaling.h"
#include "Rotation.h"
#include "Camera.h"

using namespace std;

class Matrix4
{
public:
    double val[4][4];

    Matrix4(Translation *trans);
    Matrix4(Scaling *scaling);
    Matrix4 GetRotationMatrix(Rotation *rot);
    Matrix4 GetTranspose(const Matrix4 &m);
    Matrix4 GetMcam(Camera *cam);
    Matrix4 GetMpers(Camera *cam);
    Matrix4 GetMortho(Camera *cam);
    Matrix4 GetMvp(Camera *cam);
    
    Matrix4();
    Matrix4(double val[4][4]);
    Matrix4(const Matrix4 &other);


    friend ostream &operator<<(ostream &os, const Matrix4 &m);
};

#endif