#include "Matrix4.h"
#include <iostream>
#include <iomanip>

#include "Vec3.h"
#include "Helpers.h"
#include "cmath"
#include "Scene.h"
using namespace std;

Matrix4::Matrix4()
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            this->val[i][j] = 0;
        }
    }
}
Matrix4::Matrix4(double val[4][4])
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            this->val[i][j] = val[i][j];
        }
    }
}

Matrix4::Matrix4(const Matrix4 &other)
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            this->val[i][j] = other.val[i][j];
        }
    }
}

// Custom functions
Matrix4::Matrix4(Translation *trans)
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            if(i == j)
            {
                this->val[i][j] = 1;
            }
            else{

                this->val[i][j] = 0;
            }
        }
    }
	this->val[0][3] = trans->tx;
	this->val[1][3] = trans->ty;
	this->val[2][3] = trans->tz;
}

Matrix4::Matrix4(Scaling *scaling)
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            this->val[i][j] = 0;
        }
    }
    this->val[0][0] = scaling->sx;
    this->val[1][1] = scaling->sy;
    this->val[2][2] = scaling->sz;
    this->val[3][3] = 1;
}


// Create worldToCamera matrix, M_cam
Matrix4 Matrix4::GetMcam(Camera *cam)
{
    // Copy these for clarity & easier debugging
    Matrix4 Mcam;
    
    Vec3 u(cam->u);
    Vec3 v(cam->v);
    Vec3 w(cam->w);
    Vec3 cpos(cam->pos);

    Mcam.val[0][0] = u.x;
    Mcam.val[0][1] = u.y;
    Mcam.val[0][2] = u.z;
    Mcam.val[0][3] = -(u.x*cpos.x + u.y*cpos.y + u.z*cpos.z);

    Mcam.val[1][0] = v.x;
    Mcam.val[1][1] = v.y;
    Mcam.val[1][2] = v.z;
    Mcam.val[1][3] = -(v.x*cpos.x + v.y*cpos.y + v.z*cpos.z);

    Mcam.val[2][0] = w.x;
    Mcam.val[2][1] = w.y;
    Mcam.val[2][2] = w.z;
    Mcam.val[2][3] = -(w.x*cpos.x + w.y*cpos.y + w.z*cpos.z);

    // this->val[3][0] = 0.0;
    // this->val[3][1] = 0.0;
    // this->val[3][2] = 0.0;
    Mcam.val[3][3] = 1.0;
    return Mcam;
}

Matrix4 Matrix4::GetMpers(Camera *cam)
{
    Matrix4 Mpers; // already 0's out all cells. only fill non zero ones.
    double twon = 2.0 * cam->near;

    Mpers.val[0][0] = twon / (cam->right-cam->left);
    Mpers.val[0][2] = (cam->right+cam->left) / (cam->right-cam->left);

    Mpers.val[1][1] = twon / (cam->top-cam->bottom);
    Mpers.val[1][2] = (cam->top+cam->bottom) / (cam->top-cam->bottom);
   
    Mpers.val[2][2] = (-cam->far-cam->near) / (cam->far-cam->near);
    Mpers.val[2][3] = (twon*cam->far) / (cam->far-cam->near);

    Mpers.val[3][2] = -1.0;
    return Mpers;
}    
Matrix4 Matrix4::GetMvp(Camera *cam)
{
    Matrix4 Mvp;
    Mvp.val[0][0] = cam->horRes/2.0;
    Mvp.val[0][3] = (cam->horRes-1)/2.0;

    Mvp.val[1][1] = cam->verRes/2.0;
    Mvp.val[1][3] = (cam->verRes-1)/2.0;

    Mvp.val[2][2] = 0.5;
    Mvp.val[2][3] = 0.5;
    return Mvp;
}

Matrix4 Matrix4::GetMortho(Camera *cam)
{
    Matrix4 Mortho;
    Mortho.val[0][0] = 2.0/(cam->right-cam->left);
    Mortho.val[0][3] = (-cam->right -cam->left);

    Mortho.val[1][1] = 2.0/(cam->top-cam->bottom);
    Mortho.val[1][3] = (-cam->top -cam->bottom) / (cam->top - cam->bottom);

    Mortho.val[2][2] = (-2.0)/(cam->far - cam->near);
    Mortho.val[2][3] = (-cam->far -cam->near) / (cam->far - cam->near);

    Mortho.val[3][3] = 1;
    return Mortho;
}
   
Matrix4 Matrix4::GetRotationMatrix(Rotation *rot)
{
    
    Vec3 u(rot->ux, rot->uy, rot->uz,0);
    u = normalizeVec3(u);

// Step 1) of alternate method, calculate V from u
// E.g. if u = (a, b, c) with c
// being the smallest absolute
// value then v = (-b, a, 0)
    Vec3 v;
    

    // handle rot around major axes differently.

    if(u.x == 0.0 && u.y == 0.0 && u.z != 0.0){
        // rot around Z
        cout <<"rot around Z" <<endl;

        Matrix4 rotAroundZ;
        double theta = (rot->angle* M_PI)/180.0;

        rotAroundZ.val[0][0] = cos(theta);
        rotAroundZ.val[0][1] = -sin(theta);

        rotAroundZ.val[1][0] = sin(theta);
        rotAroundZ.val[1][1] = cos(theta);

        rotAroundZ.val[2][2] = 1;
        rotAroundZ.val[3][3] = 1;
        return rotAroundZ;
    }
    else if(u.x == 0.0 && u.z == 0.0 && u.y != 0.0){
        // rot around y
        cout <<"rot around Y" <<endl;

        Matrix4 rotAroundY;
        double theta = (rot->angle* M_PI)/180.0;

        rotAroundY.val[0][0] = cos(theta);
        rotAroundY.val[0][2] = sin(theta);

        rotAroundY.val[1][1] = 1;

        rotAroundY.val[2][0] = -sin(theta);
        rotAroundY.val[2][2] = cos(theta);
        rotAroundY.val[3][3] = 1;
        return rotAroundY;
    }
    else if(u.z == 0.0 && u.y == 0.0 && u.x !=0.0){
        // rot around x
        cout <<"rot around X" <<endl;
        Matrix4 rotAroundX;
        double theta = (rot->angle* M_PI)/180.0;

        rotAroundX.val[0][0] = 1;

        rotAroundX.val[1][1] = cos(theta);
        rotAroundX.val[1][2] = -sin(theta);

        rotAroundX.val[2][1] = sin(theta);
        rotAroundX.val[2][2] = cos(theta);
        rotAroundX.val[3][3] = 1;
        return rotAroundX;
    }
    if(u.x < u.y && u.x < u.z){
        // ux is the smallest;
        v.x = 0;
        if(ABS(u.z - 0.000000) < EPSILON){
            v.y = 0;
        }
        else{
            v.y = -u.z;
        }
        v.z =  u.y;
    }
    else if(u.y < u.x && u.y < u.z){
        // uy is smallest
        v.y = 0;
        if(ABS(u.z - 0.000000) < EPSILON){
            v.x = 0;
        }
        else{
            v.x = -u.z;
        }
        // v.x = -u.z;
        v.z = u.x;
    }
    else{
        // uz is smallest
        v.z = 0;

        if(ABS(u.y - 0.000000) < EPSILON){
            v.x = 0;
        }
        else{
            v.x = -u.y;
        }
        // v.x = -u.y;
        v.y = u.x;
    }


// Step 2)  w = u x v

    Vec3 w = crossProductVec3(u, v);
    cout << "w:" << w << endl;
// Step 3) Normalize v and w

    v = normalizeVec3(v);   
    w = normalizeVec3(w);

// Form M and M^-1 matrices.

    Matrix4 mInv;
    mInv.val[0][0] = u.x;
    mInv.val[0][1] = v.x;
    mInv.val[0][2] = w.x;

    mInv.val[1][0] = u.y;
    mInv.val[1][1] = v.y;
    mInv.val[1][2] = w.y;

    // third row
    mInv.val[2][0] = u.z;
    mInv.val[2][1] = v.z;
    mInv.val[2][2] = w.z;

    mInv.val[3][3] = 1;

    // Final rot trans is: Minv * Rx(theta) * M
    double theta = (rot->angle* M_PI)/180.0;
    Matrix4 rotm;

    rotm.val[0][0] = 1.0;

    rotm.val[1][1] = cos(theta);
    rotm.val[1][2] = (-1.0)*sin(theta);

    rotm.val[2][1] = sin(theta);
    rotm.val[2][2] = cos(theta);

    rotm.val[3][3] = 1.0;

    for(int i =0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            if(rotm.val[i][j] == -0.0){
                cout << "rotm minus zero problem solved at i: " << i << " j: " << j << endl;
                rotm.val[i][j]  = 0.0;
            }
            if(mInv.val[i][j] == -0.0){
                cout << "minv minus zero problem 2 solved\n";
                mInv.val[i][j]  = 0.0;
            }
        }
    }
    cout << "mInv" << mInv << endl;
    cout << "rotm: " << rotm << endl;

    Matrix4 m = GetTranspose(mInv);

    cout << "m: " << m << endl;
    Matrix4 temp = multiplyMatrixWithMatrix(rotm, m);
    cout << "Temp: " << temp << endl;
    Matrix4 result = multiplyMatrixWithMatrix(mInv, temp);
    
    return result;
}

Matrix4 Matrix4::GetTranspose(const Matrix4 &m)
{
    Matrix4 res;
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            res.val[j][i] = m.val[i][j];
        }
    }
    return res;
}
ostream &operator<<(ostream &os, const Matrix4 &m)
{

    os << fixed << setprecision(6) << "|" << m.val[0][0] << "|" << m.val[0][1] << "|" << m.val[0][2] << "|" << m.val[0][3] << "|"
       << endl
       << "|" << m.val[1][0] << "|" << m.val[1][1] << "|" << m.val[1][2] << "|" << m.val[1][3] << "|"
       << endl
       << "|" << m.val[2][0] << "|" << m.val[2][1] << "|" << m.val[2][2] << "|" << m.val[2][3] << "|"
       << endl
       << "|" << m.val[3][0] << "|" << m.val[3][1] << "|" << m.val[3][2] << "|" << m.val[3][3] << "|";

    return os;
}