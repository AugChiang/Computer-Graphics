#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <math.h>
#include <limits> // to set infinity
#include <unordered_map> // just like python dictionary structure
#include <random> // to gen random floats for distributed ray tracing

#define PI 3.1415926
using namespace std;
const vector<double> WHITE{255,255,255};
const vector<double> BLACK{0,0,0};
//const vector<double> RED{255,0,0};
//const vector<double> GREEN{0,255,0};
//const vector<double> BLUE{0,0,255};
const double VIEW_D = 50.0; // the distance from eye to the screen
                            // arbitrary, not affecting output results, when not considering depth
const double APETURE = 15/2.2; // apeture size (diameter in mm)
const int EYE_SAMPLE_N = 1; // number of samples of eye distributed position

void Show_vec_element(vector<double> vec){
    for (auto &v : vec) {
        cout << v << " ";
    }
}

double point_d(vector<double> p1, vector<double> p2){
    double d;
    for(int i=0;i<p1.size();i++){
        d += pow((p1[i]-p2[i]),2);
    }
    return sqrt(d);
}

vector<double> cross_product(vector<double> u, vector<double> v){
    vector<double> res;
    int n = u.size();
    if(n !=v.size()){
        cerr << "The size of vectors are not consistent. (cross product)";
    }
    for(int i=1;i<2*n-2;i++){
        res.push_back(u[i%n]*v[(i+1)%n]-u[(i+1)%n]*v[i%n]);
    }
    return res;
}

double inner_product(vector<double> u, vector<double> v){
    double res;
    int n = u.size();
    if(n !=v.size()){
        cerr << "The size of vectors are not consistent. (inner product)";
    }
    for(int i=1;i<n;i++){
        res += u[i]*v[i];
    }
    return res;
}

double vec_length(vector<double> vec){
    int n=vec.size();
    double length_sq{0.0};
    for(int i=0;i<n;i++){
        length_sq += pow(vec[i],2);
    }
    return sqrt(length_sq);
}

vector<double> unit_vec(vector<double> vec){
    vector<double> u;
    int n=vec.size();
    double l = vec_length(vec);
    u.reserve(n);
    for(int i=0;i<n;i++){
        u.push_back(vec[i]/l);
    }
    return u;
}

vector<double> vec_scalar(double t, vector<double> vec){
    vector<double> u;
    u.reserve(vec.size());
    for(int i=0;i<vec.size();i++){
        u.push_back(t*vec[i]);
    }
    return u;
}

vector<double> vec_p2p_mul(vector<double> u, vector<double> v){
    vector<double> res;
    if(u.size()!=v.size()){
        cerr << "The vector size is not consistent. (sum)";
    }
    for(int i=0;i<u.size();i++){
        res.push_back(u[i]*v[i]);
    }
    return res;
}

vector<double> vec_sum(vector<double> u, vector<double> v){
    if(u.size()!=v.size()){
        cerr << "The vector size is not consistent. (sum)";
    }
    vector<double> w;
    w.reserve(u.size());
    for(int i=0;i<v.size();i++){
        w.push_back(u[i]+v[i]);
    }
    return w;
}

vector<double> vec_subtract(vector<double> u, vector<double> v){
    if(u.size()!=v.size()){
        cerr << "The vector size is not consistent. (subtract)";
    }
    vector<double> w;
    w.reserve(u.size());
    for(int i=0;i<v.size();i++){
        w.push_back(u[i]-v[i]);
    }
    return w;
}

double det33(double A3[3][3]){
    double P;
    double N;
    double res;

    P = A3[0][0]*A3[1][1]*A3[2][2]+
        A3[0][1]*A3[1][2]*A3[2][0]+
        A3[1][0]*A3[2][1]*A3[0][2];

    N = A3[0][2]*A3[1][1]*A3[2][0]+
        A3[0][1]*A3[1][0]*A3[2][2]+
        A3[0][0]*A3[2][1]*A3[1][2];

    res = P - N;
    return res;
}

vector<double> Solve_33_matrix(double A[3][3], vector<double> b){
    vector<double> coef; // to be solved
    double A_inv[3][3]; // inverse matrix of A
    // A_inv = (1/det)*adj(A)
    double c = 1/det33(A);
    if(c==0){
        coef = {0,0,0};
        return coef;
    }
    A_inv[0][0] = c*(A[1][1]*A[2][2]- A[1][2]*A[2][1]);
    A_inv[0][1] = -c*(A[0][1]*A[2][2]- A[0][2]*A[2][1]);
    A_inv[0][2] = c*(A[0][1]*A[1][2]- A[0][2]*A[1][1]);

    A_inv[1][0] = -c*(A[1][0]*A[2][2]- A[1][2]*A[2][0]);
    A_inv[1][1] = c*(A[0][0]*A[2][2]- A[0][2]*A[2][0]);
    A_inv[1][2] = -c*(A[0][0]*A[1][2]- A[1][2]*A[1][0]);

    A_inv[2][0] = c*(A[1][0]*A[2][1]- A[1][1]*A[2][0]);
    A_inv[2][1] = -c*(A[0][0]*A[2][1]- A[0][1]*A[2][0]);
    A_inv[2][2] = c*(A[0][0]*A[1][1]- A[0][1]*A[1][0]);

    // Solve coefficient: A_inv*Av = A_inv*b => v = A_inv*b
    for(int i=0;i<=2;i++){
        coef.push_back(A_inv[i][0]*b[0] +
                       A_inv[i][1]*b[1] +
                       A_inv[i][2]*b[2]  );
    }
    return coef;
}

double Array_mean(vector<double> v){
    double sum{};
    int N = v.size();
    for(size_t i{0}; i < N; i++){
        sum += v[i];
    }
    return sum/N;
}

bool Hit_Sphere(vector<double> start_pos, vector<double> ray, vector<double> sphere){
    double radius = sphere[3];
    // double t; // parameter for ray
    //        t = proj_length / ray_length;
    vector<double> centre(sphere.begin(),sphere.begin()+3); // sphere centre position
    vector<double> to_centre = vec_subtract(centre, start_pos); // vector of source to centre
    vector<double> u_to_centre = unit_vec(to_centre);
    double cos_theta; // projection length of to_ray and to_centre vectors
    double cos_phi; // tangent line and to_centre line
    vector<double> u_ray = unit_vec(ray);
    cos_theta = inner_product(u_to_centre, u_ray);
    if(cos_theta <=0){return false;}

    cos_phi = sqrt(1-pow((radius/vec_length(to_centre)),2));
    if(cos_theta > cos_phi){
        return true;
    }
    else{
        return false;
    }
}

vector<double> Sphere_Diffuse(vector<double> light_pos, vector<double> intersect_pos,
                              vector<double> centre_pos, vector<double> ambient_color){
    vector<double> diffuse_color; // reflected color on pixel (based on light color)
    vector<double> N; // normal vector of the sphere
    vector<double> L; // light vector from intersect point to light source
    diffuse_color.reserve(3);
    N.reserve(3);
    L.reserve(3);

    N = unit_vec(vec_subtract(intersect_pos, centre_pos));
    L = unit_vec(vec_subtract(light_pos, intersect_pos));

    double dot_NL = inner_product(N,L); // cos
    // cout << "dot N*L: " << dot_NL << endl;
    if(dot_NL<=0){
        return BLACK;
    }
    diffuse_color = vec_scalar(dot_NL,ambient_color);
    return diffuse_color;
}

unordered_map<string, vector<double>> Sphere_Specular(vector<double> light_pos, vector<double> intersect_pos,
                              vector<double> centre_pos, vector<double> ray, double power){
    // Specular Formula:= Ii(R dot E)^n
    unordered_map<string, vector<double>> res;
    res["color"] = BLACK;
    vector<double> specular_color; // specular component color on pixel (based on light color)
    vector<double> N; // normal vector of the sphere
    vector<double> L; // light vector from intersect point to light source
    vector<double> R; // reflected light vector
    vector<double> E; // intersected point to eye_pos
    specular_color.reserve(3);
    N.reserve(3);
    L.reserve(3);

    N = unit_vec(vec_subtract(intersect_pos, centre_pos));
    L = unit_vec(vec_subtract(light_pos, intersect_pos));



    double dot_NL = inner_product(N,L);
    // reflected light = -l + 2(l dot n)n
    R = unit_vec(vec_sum(vec_scalar(-1.0, L),
                vec_scalar(2.0, vec_scalar(dot_NL,N))));
    res["reflection"] = R;
    E = unit_vec(vec_scalar(-1.0, ray));
    double dot_RE = inner_product(R,E); // cos
    dot_RE = pow(dot_RE, power);
    // cout << "dot N*L: " << dot_NL << endl;
    if(dot_NL<=0){
        return res;
    }
    specular_color = vec_scalar(dot_RE,WHITE);
    res["color"] = specular_color;
    return res;
}

unordered_map<string, vector<double>> Chk_Sphere(vector<double> light_pos, vector<double> eye_pos,
                          vector<double> ray, vector<double> sphere, vector<vector<double>> tri_p,
                          vector<vector<double>> material, bool &is_sphere){
    // cout << "Checking Sphere Collision..." << endl;
    unordered_map<string, vector<double>> res;
    res["color"] = BLACK; // default color
    // res["reflection"] = {0,0,0};
    // res["intersection"] = {numeric_limits<double>::infinity(),
    //                        numeric_limits<double>::infinity(),
    //                        numeric_limits<double>::infinity()};
    vector<double> mat(material[sphere.back()-1]);
    vector<double> color{0,0,0};
    vector<double> amb_color(mat.begin(),mat.begin()+3); // r,g,b
    amb_color = vec_p2p_mul(WHITE, amb_color); // ambient color
    double t; // parameter for ray
    double radius = sphere[3];
    double cos_theta; // cos value of to_ray and to_centre vectors
    vector<double> centre(sphere.begin(),sphere.begin()+3); // sphere centre position
    vector<double> to_centre; // vector of source to centre
    vector<double> intersect_pos; // intersection point of sphere surface
    vector<double> chk_point; // destination point, for det intersection
    color.reserve(3);
    centre.reserve(3);
    to_centre.reserve(3);
    intersect_pos.reserve(3);
    chk_point.reserve(3);

    to_centre = vec_subtract(centre, eye_pos);

    vector<double> u_ray;
    vector<double> u_to_centre;

    u_to_centre = unit_vec(to_centre);
    u_ray = unit_vec(ray);

    // the numeric.h lib calculate inner product to int,
    // so I write it myself.
    cos_theta = inner_product(u_to_centre,u_ray);
    
    if(cos_theta <=0){ // no intersection with the sphere
        return res;
    }

    t = cos_theta * vec_length(to_centre); // scalar of the u_ray = length of projection
    
    double det{0}; // det = sqrt((x-x0)^2 + (y-y0)^2 + (z-z0)^2) v.s. radius
    chk_point = vec_sum(eye_pos, vec_scalar(t,u_ray));
    det = point_d(centre, chk_point); // length of centre to chk_point

    if(det <= radius){
        is_sphere = true;
        double backward; // used to find the 1st intersection point.
        backward = sqrt(pow(radius,2)-pow(det,2));
        // cout << "Backward length: " << backward << endl;
        // first contact point of the ray onto sphere surface
        intersect_pos = vec_sum(chk_point, vec_scalar(-backward, u_ray)); 
        res["intersection"] = intersect_pos;

        // calculate diffuse color
        vector<double> diffuse_color = Sphere_Diffuse(light_pos, intersect_pos, centre, amb_color);
        
        // ambient + diffuse
        color = vec_sum(vec_scalar(mat[3], color),
                        vec_scalar(mat[4], diffuse_color));
        // adding specular component
        unordered_map<string, vector<double>> cir_spe_res(Sphere_Specular(light_pos, intersect_pos, centre, ray, mat[6]));
        // intersection point reflection vector
        res["reflection"] = cir_spe_res["reflection"];
        // color = vec_sum(color, vec_scalar(mat[5], cir_spe_res["color"]));
        res["color"] = color;
        return res;
    }
    else{
        is_sphere = false;
        return res;
    }
}

vector<double> Tri_diffuse(vector<double> light_pos, vector<double> intersect_pos,
                           vector<double> normal_vec, vector<double> ambient_color){
    vector<double> ray = vec_subtract(intersect_pos, light_pos);
    vector<double> diffuse_color; // reflected color on pixel
    vector<double> L; // light vector from intersect point to light source
    vector<double> N;
    diffuse_color.reserve(3);
    L.reserve(3);
    N.reserve(3);

    L = unit_vec(vec_subtract(light_pos, intersect_pos));
    N=unit_vec(normal_vec);

    double dot_NL = inner_product(N,L);
    // cout << "dot N*L: " << dot_NL << endl;
    if(dot_NL<=0){
        return BLACK;
    }
    diffuse_color = vec_scalar(dot_NL,ambient_color);
    return diffuse_color;
}

unordered_map<string, vector<double>> Tri_specular(vector<double> light_pos, vector<double> intersect_pos,
                            vector<double> normal_vec, vector<double> ray, double power){
    unordered_map<string, vector<double>> res;
    res["color"] = BLACK;
    vector<double> specular_color;
    vector<double> N; // normal vector of the sphere
    vector<double> L; // light vector from intersect point to light source
    vector<double> R; // reflected light vector
    vector<double> E; // intersected point to eye_pos
    specular_color.reserve(3);
    N.reserve(3);
    L.reserve(3);

    N = unit_vec(normal_vec);
    L = unit_vec(vec_subtract(light_pos, intersect_pos));

    double dot_NL = inner_product(N,L);
    // reflected light = -l + 2(l dot n)n
    R = unit_vec(vec_sum(vec_scalar(-1.0, L),
                vec_scalar(2.0, vec_scalar(dot_NL,N))));
    res["reflection"] = R;
    E = unit_vec(vec_scalar(-1.0, ray));
    double dot_RE = inner_product(R,E);
    dot_RE = pow(dot_RE, power);
    // cout << "dot N*L: " << dot_NL << endl;
    if(dot_NL<=0){
        return res;
    }
    specular_color = vec_scalar(dot_RE,WHITE);
    res["color"] = specular_color;
    return res;
}


unordered_map<string, vector<double>> Chk_Triangle(vector<double> light_pos, vector<double> eye_pos,
                            vector<double> ray, vector<double> sphere, vector<vector<double>> tri_p,
                            vector<vector<double>> material){
    unordered_map<string, vector<double>> res;
    res["color"] = BLACK;
    // res["reflection"] = {0,0,0};
    // res["intersection"] = {numeric_limits<double>::infinity(),
    //                        numeric_limits<double>::infinity(),
    //                        numeric_limits<double>::infinity()};

    // tri_p structure: [ [9 points + mat_No.], [9 points + mat_no.], ... ]
    int tri_num = tri_p.size()-1;
    // cout << "Number of triangle: " << tri_num << endl;

    // origin point, take one of the three triangle vertices.
    vector<double> V0; 
    // use for linear combination to represent the triangle plane
    vector<double> v1;
    vector<double> v2;
    vector<double> n; // normal vector of the triangle
    vector<double> chk_point;
    V0.reserve(3);
    v1.reserve(3);
    v2.reserve(3);
    n.reserve(3);
    chk_point.reserve(3);

    for(int i=0;i<tri_num;i++){ // i: index of the triangle.
        // material info
        vector<double> mat = material[tri_p[i].back()-1]; // material No.
        vector<double> amb_color(mat.begin(),mat.begin()+3); // r,g,b ratio
        vector<double> color;
        color.reserve(3);
        amb_color = vec_p2p_mul(WHITE, amb_color); // ambient color

        // init
        V0.clear();
        v1.clear();
        v2.clear();
        n.clear();

        // cout << "i-th: " << i << endl;
        V0.push_back(tri_p[i][0]); // x<i0>
        V0.push_back(tri_p[i][1]); // y<i0>
        V0.push_back(tri_p[i][2]); // z<i0>

        v1.push_back(tri_p[i][3]-tri_p[i][0]);
        v1.push_back(tri_p[i][4]-tri_p[i][1]);
        v1.push_back(tri_p[i][5]-tri_p[i][2]);

        v2.push_back(tri_p[i][6]-tri_p[i][0]);
        v2.push_back(tri_p[i][7]-tri_p[i][1]);
        v2.push_back(tri_p[i][8]-tri_p[i][2]);

        n = cross_product(v1,v2);
        if(n[2]<0){ n = vec_scalar(-1.0,n);} // normal vector points up to +y axis
        n = unit_vec(n); // unit normal vector of the triangle(↑)
        if(inner_product(n,ray) == 0){
            return res;
        }
        // Calculate the intersection
        /*
            Plane: (x,y,z) = V0 + s1*(vec1) + s2*(vec2)
            Ray:   (x,y,z) = (eye_pos) + t(ray_vec)
            ==> s1*(vec1)+ s2*(vec2)- t(ray_vec) = eye- V0
            ==> Av = b
            Intersection found if 0<= s1, s2 <= 1 and t > 0
        */
        double equ_arr[3][3]; // 3*3 matrix of v1, v2, ray.
        vector<double> b; // values after rearrange the equation.
        vector<double> coef{0,0,0}; // tri vec coefficient.
        b.reserve(3);

        // b = O - V0
        b = vec_subtract(eye_pos,V0); // eye_pos - V0
        
        // A matrix
        for(int j=0;j<=2;j++){ // j-th row
        // fill matrix: right then down
        // matrix = [ vec1(↓)  vec2(↓)  -ray(↓)]
            equ_arr[j][0] = v1[j]; // i-th tri_vec, x1 -> y1 -> z1
            equ_arr[j][1] = v2[j]; // i-th tri_vec, x2 -> y2 -> z2
            equ_arr[j][2] = -ray[j];
        }
        coef = Solve_33_matrix(equ_arr,b); // return coefficients of s1, s2, t
        // cout << "(s1, s2, t) =: " << coef[0] <<", " << coef[1] <<", " << coef[2] << endl;

        // if t<0 then not intersect.
        if(coef[2]<=0){
            return res;
        }
        else{ // t>=0 and 0 <= {s1,s2} <=1 and s1+s2 <=1
            if(1>=coef[0]>=0 && 1>=coef[1]>=0){
                if(coef[0]+coef[1] <= 1 ){
                    //chk_point = vec_sum(vec_sum(V0, vec_scalar(coef[0],v1)), vec_scalar(coef[1],v2));
                    chk_point = vec_sum(eye_pos, vec_scalar(coef[2],ray));
                    res["intersection"] = chk_point;

                    // shading trace
                    vector<double> shading_ray = vec_subtract(light_pos, chk_point);
                    bool is_shade = Hit_Sphere(chk_point, shading_ray, sphere);
                    if(is_shade){
                        return res;
                    }
                    else{ // not blocked by the sphere
                        unordered_map<string, vector<double>> tri_spe_res;
                        vector<double> diffuse_color = Tri_diffuse(light_pos, chk_point, n, amb_color);
                        tri_spe_res = Tri_specular(light_pos, chk_point, n, ray, mat[6]); // color, reflection
                        color = vec_sum(vec_scalar(mat[3], amb_color),
                                        vec_scalar(mat[4], diffuse_color));
                        color = vec_sum(color, vec_scalar(mat[5],tri_spe_res["color"]));
                        res["color"] = color;
                        return res;
                    }
                }
                else{
                    return res;
                }
            }
        }
    }

    return res;
}

vector<double> Reflect_from_tri(vector<double> intersect_pos, vector<double> reflect_ray,
                                vector<double> ){
    vector<double> refle_color;
    return refle_color;
}


int main(){
    // Read the file
    string path = "hw3_input.txt";
    // spacing axes: x←, y↑ , Z↗
    // claim configuration vars
    vector<double> eye_pos; // eye position
    vector<double> light_pos; // light position
    vector<double> view_pos; // view direction & up direction
    eye_pos.reserve(3);
    light_pos.reserve(3);
    view_pos.reserve(6);

    double field_ang{0}; // field angle
    int mat_no{0}; // material number for recording which objects are applied No. texture.
    vector<int> size; // resolution, width, height
    size.reserve(2);

    vector<vector<double>> mat; // materials coefficients
    vector<vector<double>> sphere; // centre position, radius
    vector<vector<double>> tri_p; // 1 triangle defined by every 3 points

    ifstream Input_File(path);

    // import input values to variables
    if (Input_File.is_open()){
        cout << "Start reading file..." << endl;
        string line;
        char header;

        vector<double> mat_temp;
        vector<double> tri_temp;
        vector<double> sphere_temp;
        tri_temp.reserve(9);

        while (getline(Input_File, line)) {
            istringstream str2value;
            string str_value;

            // cout << line << endl;
            if (line[0] == '#' || line.empty()) continue; // empty line

            header = line[0];

            str2value.str(line.substr(1)); // put 'line' (except header) to 'istringstream' obj  
            while(str2value >> str_value){ //read every word, delimited by space
                // cout << str_value << endl;
                switch(header){
                    case 'E': // eye position
                        // cout << "String value: "<<str_value << endl;
                        // cout << "Converted value: "<<stof(str_value) << endl;
                        eye_pos.push_back(stof(str_value));
                        break;
                    case 'L': // light position
                        light_pos.push_back(stof(str_value));
                        break;
                    case 'V': // view direction
                        view_pos.push_back(stof(str_value));
                        break;
                    case 'F': // field angle
                        field_ang = stoi(str_value);
                        // cout << "Filed angle:" << field_ang << endl;
                        break;
                    case 'R': // resolution: width, height
                        size.push_back(stoi(str_value));
                        break;
                    case 'M': // mats
                        mat_temp.push_back(stof(str_value));
                        break;
                    case 'S': // sphere
                        sphere_temp.push_back(stof(str_value));
                        break;
                    case 'T': // triangles
                        tri_temp.push_back(stof(str_value));
                        break;
                }
            }
            if(header == 'M'){
                mat_no += 1;
                // cout << "Mat No." << mat_no << endl;
                mat.push_back(mat_temp);
                mat_temp.clear();
            }

            if(header == 'S'){
                sphere_temp.push_back(mat_no); // record texture No.
                sphere.push_back(sphere_temp);
                // cout << "Mat No of sphere: " << sphere.back() << endl;
                sphere_temp.clear();
            }
            if(header == 'T'){
                tri_temp.push_back(mat_no); // record texture No.
                tri_p.push_back(tri_temp);
                // cout << "Mat No of triangle: " << tri_temp.back() << endl;
                tri_temp.clear();
            }
        }
        cout << "Importing complete." << endl;

    }
    else {
        std::cerr << "Please check the input file.\n";
    }
    // cout << "Tri pos size: "<< tri_p.size() << tri_p[0].size() << tri_p[1].size() << endl;
    // cout << "Sphere size: "<< sphere.size() << sphere[0].size() << sphere[1].size() << sphere[2].size() << endl;

//============================= Parameter Setting ======================================

    int WIDTH = size[0];
    int HEIGHT = size[1];

    double FIELD_W{0}; // output screen width
    double FIELD_H{0}; // output screen height
    FIELD_W = 2*VIEW_D*tan((field_ang/2)*PI/180.0);
    FIELD_H = FIELD_W*HEIGHT/WIDTH;
    
    vector<double> view_vec(view_pos.begin(),view_pos.begin()+3); // eye-viewing direction vector
    vector<double> up_vec(view_pos.begin()+3,view_pos.end()); // prependicular with the ground
    vector<double> eye_right_vec;
    vector<double> eye_down_vec;
    eye_right_vec.reserve(3);
    eye_down_vec.reserve(3);

    vector<double> output_o; // output area origin
    output_o.reserve(3);

    // vector for output screen.
    vector<double> dw;
    vector<double> dh;
    dw.reserve(3);
    dh.reserve(3);


    // vector from light to pixel
    vector<double> ray;
    ray.reserve(3);
    
    eye_right_vec = cross_product(view_vec, up_vec);
    // assume: length of view_vec = distance from eye to output screen
    view_vec = vec_scalar(VIEW_D, unit_vec(view_vec));
    // cout << "Vector of view: " << view_vec[0]<< " " << view_vec[1] <<" "<< view_vec[2] << endl;
    
    vector<double> u_eye_right_vec(unit_vec(eye_right_vec));
    eye_right_vec = vec_scalar(FIELD_W/2, u_eye_right_vec);
    // cout << "Vector of left: " << eye_right_vec[0]<< " " << eye_right_vec[1] <<" "<< eye_right_vec[2] << endl;
    
    vector<double> view_right_sum(vec_sum(view_vec, eye_right_vec));
    // cout << "Vector of left + view: " << view_right_sum[0]<< " " << view_right_sum[1] <<" "<< view_right_sum[2] << endl;
    
    eye_down_vec = cross_product(view_right_sum, eye_right_vec);
    vector<double> u_eye_down_vec(unit_vec(eye_down_vec));
    eye_down_vec = vec_scalar(FIELD_H/2, u_eye_down_vec);
    // cout << "Vector of DOWN: " << eye_down_vec[0] <<" "<< eye_down_vec[1]<<" " << eye_down_vec[2] << endl;
    
    output_o = vec_sum(eye_pos, vec_sum(view_right_sum, eye_down_vec)); // down-right corner
    // cout << "Output origin: " << output_o[0] <<" "<< output_o[1]<<" " << output_o[2] << endl;

    dw = vec_scalar(-1.0, u_eye_right_vec); // unit vec(←)
    dh = vec_scalar(-1.0, u_eye_down_vec); // unit vec(↑)
    // cout << "Unit vector of output screen x : " << dw[0] <<" "<< dw[1] <<" "<< dw[2] << endl;
    // cout << "Unit vector of output screen y : " << dh[0] <<" "<< dh[1] <<" "<< dh[2] << endl;

    // image.resize(HEIGHT,vector<vector<double>>(WIDTH,vector<double>(3))); // init all pixel rgb = (0,0,0)
    double pw, ph; // width and height of each pixel on the output screen
    pw = FIELD_W / WIDTH;
    ph = FIELD_H / HEIGHT;
    // cout << "Pixel width and height:" << pw << ph << endl;

    // Create and open a text file
    ofstream Output_File("output.ppm");
    // Write to the file
    /* ppm rgb format:
    P3           # "P3" means this is a RGB color image in ASCII
    3 2          # "3 2" is the width and height of the image in pixels
    255          # "255" is the maximum value for each color
    */
    Output_File << "P3" << endl; // for .ppm file format
    Output_File << HEIGHT <<" "<< WIDTH << endl;
    Output_File << "255" << endl;

//================================ Ray Tracing =====================================
    cout << "Start Ray Tracing..." << endl;

    uniform_real_distribution<double> uni_f(0,1); // gen random float for distributed ray tracing
    uniform_int_distribution<int> uni_int(0,360); // random int for distributed angle.
    default_random_engine re;
    //====================== pixel check ================================
    for(int row=HEIGHT;row>0;row--){
        for(int col=WIDTH;col>0;col--){
            unordered_map<string, vector<double>> res;
            vector<double> pixel_pos; // pixel position = output_origin + dx(→) + dy(↓)
            vector<double> pixel_shift; // dx(→) + dy(↓)
            vector<double> pixel_color;
            vector<double> ray;
            vector<double> reflect_ray;
            pixel_shift.reserve(3);
            pixel_pos.reserve(3);
            pixel_color.reserve(3);
            ray.reserve(3);

            vector<double> dx(vec_scalar(col*pw, dw)); // number of pixel toward x-axis (width)
            vector<double> dy(vec_scalar(row*ph, dh)); // number of pixel toward y-axis (height)

            pixel_shift = vec_sum(dx,dy);
            pixel_pos = vec_sum(output_o, pixel_shift);  // pixle position

            ray = vec_subtract(pixel_pos, eye_pos);
            
            // distributed ray tracing, i.e. distributed eye positions in apeture range
            vector<double> r_arr{};
            vector<double> g_arr{};
            vector<double> b_arr{};
            r_arr.reserve(EYE_SAMPLE_N);
            g_arr.reserve(EYE_SAMPLE_N);
            b_arr.reserve(EYE_SAMPLE_N);
            int r{}, g{}, b{}; // init RGB
            //====================== distrubted ray tracing ===========================
            for(size_t n{0}; n < EYE_SAMPLE_N; n++){
                // init RGB for each sampled eye ray
                r=0;
                g=0;
                b=0;

                // init eye position
                eye_pos[0] = 0;
                eye_pos[1] = 0;

                // random radisu and theta for distributed.
                double rand_r{(APETURE/20)*uni_f(re)}; // notice: apeture unit is mm, convert to unit:cm here.
                int rand_theta{uni_int(re)}; // random angle unit:degree

                eye_pos[0] += rand_r*cos(rand_theta*PI/180);
                eye_pos[1] += rand_r*sin(rand_theta*PI/180);

                //======================= object rendering =============================
                // assign the pixel color
                int num_sphere = sphere.size(); // No. of the spheres
                bool is_sphere{false};
                for(size_t s{0}; s<num_sphere ;s++){ // check every spheres first.
                    if(r == 0 && g == 0 && b == 0){
                        res = Chk_Sphere(light_pos, eye_pos, ray, sphere[s], tri_p, mat, is_sphere);
                        pixel_color = res["color"];
                        // check global reflection with other objects
                        reflect_ray = res["reflection"];
                    }

                    if(!is_sphere){
                        res = Chk_Triangle(light_pos, eye_pos, ray, sphere[s], tri_p, mat);
                        pixel_color = res["color"];
                    }
                    // assign RGB value
                    r = round(pixel_color[0]);
                    g = round(pixel_color[1]);
                    b = round(pixel_color[2]);
                }

                r_arr.push_back(pixel_color[0]);
                g_arr.push_back(pixel_color[1]);
                b_arr.push_back(pixel_color[2]);
                
            }
            r = round(Array_mean(r_arr));
            g = round(Array_mean(g_arr));
            b = round(Array_mean(b_arr));
            // RGB value adjustment
            if(r>255) r=255;
            if(g>255) g=255;
            if(b>255) b=255;
            if(r<0) r=0;
            if(g<0) g=0;
            if(b<0) b=0;
            // write RGB value to the ppm file
            Output_File << r << " "
                        << g << " "
                        << b << " \n";
            pixel_color.clear();

        }
    }

    cout << "File has been written." << endl;

    // Close the file
    Output_File.close();

    return 0;
}