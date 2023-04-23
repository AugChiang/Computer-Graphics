#include<vector>
#include<cmath>
#include<math.h>
#include<iostream>
#include"algebra3.h"
using namespace std;
#define PI 3.1415926


mat3 Rot_mat(double radian, int power){
    // using static_cast to fix the narrowing conversion error
    double cos_theta = cos(power*radian);
    double sin_theta = sin(power*radian);

    // double a{cos(radian)};
    vec3 col0{cos_theta, 0, sin_theta};
    vec3 col1{0,1,0};
    vec3 col2{-sin_theta, 0, cos_theta};
    mat3 rotate_mat{col0,col1,col2};

    return rotate_mat;
}

mat3 Mat_mul_mat(mat3 mat1, mat3 mat2){
    mat3 res{};
    for(size_t i{0}; i<3; i++){
        for(size_t j{0}; j<3; j++){
            for(size_t k{0}; k<3 ; k++){
                res[i][j] = mat1[i][k]*mat2[k][j] +
                            mat1[i][k]*mat2[k][j] +
                            mat1[i][k]*mat2[k][j];
            }
        }
    }
    return res;
}

vec3 Mat_mul_vec(mat3 M, vec3 v){
    vec3 res{};
    for(size_t i{0}; i<3; i++){
        for(size_t j{0}; j<3; j++){
                res[i] += M[i][j]*v[j];
        }
    }
    return res;
}

vec3 Vec_p2p_mul(vec3 u, vec3 v){
    return vec3{ u[0]*v[0], u[1]*v[1], u[2]*v[2] };
}

vector<double> Smooth_step(double a, double b, double x){
    vector<double> res{}; // value, dx
    res.reserve(2);

    double l = min(1.0, max(0.0, (x-a)/(b-a)));
    res.push_back(3*pow(l,2) - 2*pow(l,3)); // value

    double dS = 6*(l-pow(l,2)); // the slope of smooth step at x
    res.push_back(dS); // derivative
    return res; // value, dx
}

double Gen_coef(int i, int j){
    double int_part{};
    double u = 50*modf(i/PI, &int_part);
    double v = 50*modf(j/PI, &int_part);
    // double u = 50*static_cast <double> (rand()) / static_cast <double> (RAND_MAX)+0.01;
    // double v = 50*static_cast <double> (rand()) / static_cast <double> (RAND_MAX)+0.01;
    double a = 2*modf(u*v*(u+v), &int_part)-1;

    return a;
}

vector<double> Generator(vec3 pos_xz){ // pos = (x,z)
    vector<double> res{};
    res.reserve(3); // y-value, dx, dz
    // scale such that one unit of space is nearly 1 meter.
    double x = pos_xz[0];
    double z = pos_xz[2];
    const double scale_factor{600};
    const double bias{600};
    // (i,k) = grid position
    int i = floor(x);
    int j = floor(z);

    double a = Gen_coef(i,j);
    double b = Gen_coef(i+1,j);
    double c = Gen_coef(i,j+1);
    double d = Gen_coef(i+1,j+1);
    // double rand_x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX)+x/2000;
    // double rand_z = static_cast <double> (rand()) / static_cast <double> (RAND_MAX)+z/2000;
    // derivatives of x,z
    double dx{0};
    double dz{0};

    vector<double> smooth_x = Smooth_step(0,1,x-i); // value, dx
    vector<double> smooth_z = Smooth_step(0,1,z-j); // value, dz
    double y = a + (b-a)*smooth_x[0] +
                   (c-a)*smooth_z[0] +
                   (a-b-c+d)*smooth_x[0]*smooth_z[0];

    dx = scale_factor*((b-a)+(a-b-c+d)*smooth_z[0])*smooth_x[1];
    dz = scale_factor*((c-a)+(a-b-c+d)*smooth_x[0])*smooth_z[1];

    // double y = a + (b-a)*Smooth_step(0,1,rand_x) + (c-a)*Smooth_step(0,1,rand_z)
    //             + (a-b-c+d)*Smooth_step(0,1,rand_x)*Smooth_step(0,1,rand_z);
    y = scale_factor*y+bias;
    res.push_back(y);
    res.push_back(dx);
    res.push_back(dz);
    return res; //value ,dx ,dz
}

vector<double> Gen_detail(vec3 pos_xz, int n){
    vector<double> res{};
    res.reserve(3); // y-value, dx, dz
    // k to control the extent of details.
    
    // default rotation matrix:
    // 3-4-5 triangle in radians are  [0.64, 0.93, and 1.57]
    vector<double> y_dx_dz{}; // passed by result of generator
    double y{0};
    double dx{0};
    double dz{0};
    
    for(size_t k{0}; k<n+1; k++){ // add high frequency part
        mat3 rotM = Rot_mat(0.64, int(k));
        vec3 rot = Mat_mul_vec(rotM, pos_xz);
        vec3 xz = rot*pow(2,k);
        y_dx_dz = Generator(xz);

        y += (1/pow(2,k))*y_dx_dz[0]; // y-value
        dx += (1/pow(2,k))*y_dx_dz[1]; // dx
        dz += (1/pow(2,k))*y_dx_dz[2]; // dz
    }
    res.push_back(y); // res[0]
    res.push_back(dx); // res[1]
    res.push_back(dz); // res[2]
    return res; //y-value, dx, dz
}

vec3 Shader(vec3 pos_xyz, vec3 v2sun,
            double coef_t=16, double coef_a=1e-3){
    // check only the pixels in the shadow.
    vec3 new_xz{pos_xyz}; // new position along sun direction
    double new_y{}; // searching height
    double s{1}; // steps record
    double dh{0}; // difference of heights, init as pos_xz's height
    // cout << "res: " << coef_a*(dh/coef_t) << endl;
    for(size_t i{1}; i<=coef_t; i++){
        new_xz += v2sun; // one step towards sun direction per loop.
        new_y = Generator(new_xz)[0];
        if(new_y - pos_xyz[1] > dh){ // update dh, and retain local maximum one.
            dh = new_y - pos_xyz[1];
            s = int(i);
        }
    }
    // cout << "shader coefficient: " << coef_a*(dh/s) << endl;
    return exp(-coef_a*(dh/s));

}

vec3 Reflect2Eye(vec3 normal_gen, vec3 p2sun, vec3 p2eye,
                   vec3 sun_light_color, double highlight_coef=0.95){
    vec3 reflection_ray = (2*normal_gen*(normal_gen*p2sun) - p2sun).normalize();
    double highlight_intensity = 1-5*sqrt(0.5*(1+(p2sun*p2eye)));
    double r_dot_v{}; // reflection_ray dot vector-to-viewer
    vec3 res{0};

    r_dot_v = pow(reflection_ray*p2eye, 8);
    if(r_dot_v > 0 && (p2eye*p2sun) > 0){
        res = r_dot_v*((1-highlight_coef) + highlight_coef*highlight_intensity)
                     *sun_light_color;
        return res;
    }else{
        return res;
    }
}