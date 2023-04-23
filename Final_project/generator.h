#include<vector>
#include<math.h>
#include<cmath>
#include<iostream>
#include"algebra3.h"
#define PI 3.1415926
using namespace std;
/*
    The main function of generating terrain or clouds.
    To ensure the connectivity of each block at (i,j),
    another function to handle this, noted as 'S' below, i.e. smooth step function
    and x,z is the continuous range of the position,
    where a,b,c,d are coefficients.

        f[i,j](x,z) = a[i,j] + (b[i,j]-a[i,j])*S(x-i)
                             + (c[i,j]-a[i,j])*S(z-j)
                             + (a[i,j]-b[i,j]-c[i,j]+d[i,j])*S(x-i)*S(z-j)
        S(a,b,x) = 3l^2 - 2l^3
        where l = min(1,max(0,(x-a)/(b-a)))
*/

vector<double> Generator(vec3 pos_xz); // y-value, dx, dz
// summation of generator, i.e. adding high frequency components.
vector<double> Gen_detail(vec3 pos_xz, int n); // outputs: y-value, dx ,dz
vector<double> Smooth_step(double a, double b, double x); // outputs: value, dx
double Gen_coef(int i, int j); // coefficients for generator

vec3 Mat_mul_vec(mat3 M, vec3 v); // 3x3 matrix multiplication with vec3

// soft-shadow function depends on the difference distance from sun-direction hill
vec3 Shader(vec3 pos_xyz, vec3 v2sun,
            double coef_t=1, double coef_a=1e-3);
mat3 Rot_mat(double radian=0, int power=1); // set a 3x3 rotation matrix
mat3 Mat_mul_mat(mat3 mat1, mat3 mat2); // 3x3 matrix multiplication with another one
vec3 Vec_p2p_mul(vec3 u, vec3 v); // point dot multiplcation of vectors
vec3 Reflect2Eye(vec3 normal_gen, vec3 p2sun,
                 vec3 p2eye, vec3 sun_light_color,
                 double highlight_coef=0.95); // reflection light to viewer(eyes)

