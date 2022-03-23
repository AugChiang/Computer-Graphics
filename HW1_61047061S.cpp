#include<string>
#include<vector>
#include<array>
#include<istream>
#include<fstream>
#include<iterator>
#include<iostream>
#include<sstream>
#include<numeric>
#include<cmath>

using namespace std;

void Write_to_pbm(vector<vector<int>> output){
    int rows = output.size();
    int cols = output[0].size();

    // Create and open a text file
    ofstream Output_File("output.pbm");
    // Write to the file
    Output_File << "P1" << endl; // for .pbm file format
    Output_File << rows <<" "<< cols << endl;

    for(int i=0;i<rows;i++){
        for(int j=0;j<cols;j++){
            Output_File << output[i][j] << " ";
        }
        Output_File << endl;
    }
    cout << "File has been written." << endl;
    // Close the file
    Output_File.close();
}

void Input_to_vec(string input_line, vector<double> &coordinates){
    // cout << "Input by line: " << input_line << endl;
    istringstream line(input_line);
    string input_pos;
    while(line>>input_pos){
        // object from the class stringstream to convert string into double
        double pos;
        stringstream x(input_pos);
        x >> pos;
        // cout << "The element that push into vec: " << pos << endl;
        // this always push a zero into coordinates because of the first char
        coordinates.push_back(pos);

    }
    coordinates.erase(coordinates.begin()); // erase the first element in the vec.
}

int Chk_Sphere(vector<double> eye, vector<double> ray, vector<double> to_centre, double radius){
    // cout << "Checking Sphere Collision..." << endl;
    double t; // parameter for ray
    double dot_product; // projection length of to_ray and to_centre vectors
    double ray_length;
    double proj_length; // length of to_centre vec project on ray_vec
    vector<double> chk_point; // destination point, (x,y,z)

    dot_product = inner_product(to_centre.begin(), to_centre.end(), ray.begin(), 0);
    ray_length = sqrt(pow(ray[0],2) + pow(ray[1],2) + pow(ray[2],2));
    proj_length = dot_product / ray_length;
    t = proj_length / ray_length;
    
    double det{0}; // det = (x-x0)^2 + (y-y0)^2 + (z-z0)^2 v.s. r**2
    for(int i=0;i<=2;i++){
        chk_point.push_back(eye[i] + t*(ray[i]/sqrt(ray_length)));
        // cout << "Destination Point: " << chk_point[i] << endl;
        det = det + pow((chk_point[i] - (eye[i] + to_centre[i])), 2);
    }
    // cout << "Result of Sphere Collision Calculation: " << det << endl;
    // cout << "Collision Result: " << det << endl;
    chk_point.clear();

    if(det <= pow(radius, 2)){
        return 0;
    }
    else{
        return 1;
    }
}


double Det33(double A3[3][3]){
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

vector<double> Solve_33_matrix(double A[3][3], double b[3]){
    vector<double> coef; // to be solved
    double A_inv[3][3];
    // A_inv = (1/det)*adj(A)
    double c = 1/Det33(A);
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
        coef.push_back(A_inv[i][0]*b[0]+A_inv[i][1]*b[1]+A_inv[i][2]*b[2]);
    }
    return coef;
}


int Chk_tri(vector<double> eye, vector<double> tri_pos, vector<double> ray){
    int tris = tri_pos.size()/9; // number of triangles
    vector<vector<double>> tri_vec; // store sets of two paired vectors of triangle.
    // generate vectors of triangles.
    for(int i=0;i<tris;i++){
        vector<double> tri_vec_temp;
        // tri_vec (x,y,z)
        double tri_x;
        double tri_y;
        double tri_z;

        for(int j=1;j<=2;j++){ // only two vectors for a triangle.
            // two vectors of the triangle, in the i-th triangle points set.
            tri_x = tri_pos[9*i+3*j]-tri_pos[9*i]; // xj - x0
            tri_y = tri_pos[9*i+3*j+1]-tri_pos[9*i+1]; // yj- y0
            tri_z = tri_pos[9*i+3*j+2]-tri_pos[9*i+2]; // zj - z0
            tri_vec_temp.push_back(tri_x);
            tri_vec_temp.push_back(tri_y);
            tri_vec_temp.push_back(tri_z);
        }
        tri_vec.push_back(tri_vec_temp);
        tri_vec_temp.clear();
    }
    // calculate intersection
    /*
    Plane: (x,y,z)= V0 + s1*(vec1)+ s2*(vec2)
    Ray:   (x,y,z)= (eye)+ t(ray_vec)
    ==> s1*(vec1)+ s2*(vec2)- t(ray_vec) = eye- V0
    ==> Av = b
    Intersection found if 0<= s1, s2 <= 1 and t > 0
    */
    double equ_arr[3][3]; // 3*3 matrix of v1, v2, ray.
    vector<double> coef{0,0,0}; // tri vec coefficient.
    for(int i=0;i<tri_vec.size();i++){
        // tri_vec structure: each row = [x1,y1,z1,x2,y2,z2]
        double b[3]; // values after rearrange the equation.

        // b = O - V0
        b[0] = eye[0]-tri_pos[9*i];
        b[1] = eye[1]-tri_pos[9*i+1];
        b[2] = eye[2]-tri_pos[9*i+2];
        
        // A matrix
        for(int j=0;j<=2;j++){
            equ_arr[j][0] = tri_vec[i][0+j]; // i-th tri_vec, x1 -> y1 -> z1
            equ_arr[j][1] = tri_vec[i][3+j]; // i-th tri_vec, x2 -> y2 -> z2
            equ_arr[j][2] = -ray[j];
        }

        coef = Solve_33_matrix(equ_arr,b); // return coefficients of s1, s2, t

        // cout << "alpha, beta, t =: " << coef[0] <<", " << coef[1] <<", " << coef[2] << endl;


        // if t<0 then not intersect.
        if(coef[2]<=0){
            return 1;
        }
        else{ // 0 <= s1,s2 <=1 and s1+s2 <=1
            if(1>=coef[0]>=0 and 1>=coef[1]>=0){
                if(coef[0]+coef[1] <= 1 ){
                    return 0;
                }
                else{
                    return 1;
                }
            }
            
        }
    }
    return 0;
}

int Ray_tracing(int width, int height, vector<vector<double>> &env){
    // cout << "Ray_tracing function activated." << endl;
    // claim the 2D output array with init value = 0.
    vector<vector<int>> output(height, vector<int>(width,1)); // init output = 1
    // eye coordinates
    vector<double> EYE = {env[0][0], env[0][1],env[0][2]};
    // screen vec of width
    double OUTPUT_w_vec = (env[1][3] - env[1][0])/width;
    // double OUTPUT_wy = (env[1][4] - env[1][1])/width;
    // double OUTPUT_wz = (env[1][5] - env[1][2]);
    // screen vec of height
    double OUTPUT_h_vec = (env[1][7] - env[1][1])/height;
    // double OUTPUT_hy = (env[1][7] - env[1][1])/height;
    // double OUTPUT_hz = (env[1][5] - env[1][2]);
    vector<double> ORIGIN = {env[1][0], env[1][1], env[1][2]};
    vector<double> to_centre = {env[2][0] - EYE[0],
                               env[2][1] - EYE[1],
                               env[2][2] - EYE[2]};
    vector<double> tri_pos = env[3]; // triangle points, 9 elements(3 points) as 1 triangle.
    vector<double> ray; // eye-to-pixel vector

    double r = env[2][3]; // sphere radius

    for(int i=0; i<height; i++){ // height
        for(int j=0; j<width; j++){ // width
            // cout << "Width Vec: "<< OUTPUT_w_vec << endl;
            // cout << "Height Vec: "<< OUTPUT_h_vec << endl;

            double dw = i*OUTPUT_w_vec;
            double dh = j*OUTPUT_h_vec;
            // cout << "Delta w&h = "<< dw << "," << dh << endl;
            double ray_x = (ORIGIN[0] + dw) - EYE[0];
            double ray_y = (ORIGIN[1] + dh) - EYE[1];
            double ray_z = ORIGIN[2] - EYE[2];

            // cout <<"Ray: "<< ray_x << "," << ray_y << "," << ray_z << endl;
            ray = {ray_x, ray_y, ray_z};
            int res; // to fill output
            res = Chk_Sphere(EYE, ray, to_centre, r);
            output[i][j] = res;
            if (output[i][j]==1){
                res = Chk_tri(EYE, tri_pos, ray);
                output[i][j] = res;
            }
            ray.clear();
            // cout << res << endl;
            // cout << "Output("<<i<<","<<j<<") = " << output_temp[j] << endl;
        }
    }
    // cout << "Size of the output matrix: "<< output.size()<<"*"<<output[0].size() << endl;
    Write_to_pbm(output);
    return 0;
}

int main(){

    cout << "Start Reading the file..." << endl;
    // Create a text string, which is used to output the text file
    string input_line;

    // resolution
    vector<int> reso;

    // vector container for the vectors below
    vector<vector<double>> env;
    // vectors for double elements
    vector<double> eye_pos; // (x,y,z)
    vector<double> output_area; // 4 points
    vector<double> sphere_pos; // centre(1 points), radius(real number)
    vector<double> square_pos; // square
    vector<double> tri_pos; // triangle points, 3 as 1 triangle

    // Read from the text file
    ifstream InputFile("hw1_input.txt"); // claim the var

    // Use a while loop together with the getline() function to read the file line by line
    while (getline (InputFile, input_line)){
        vector<double> tri_pos_temp;
        char header = input_line[0]; // distinguish the pos types.
        // istringstream line(input_line); // create iterative obj 'line', not including space.
        // string input_pos; // for storing cooridnates string
        switch(header){
            // Eye position
            case 'E':
                Input_to_vec(input_line, eye_pos);
                break;
            // Output
            case 'O':
                Input_to_vec(input_line, output_area);
                break;
            // resolution
            case 'R':
                break;
            // sphere
            case 'S':
                Input_to_vec(input_line, sphere_pos);
                break;
            case 'T':
                Input_to_vec(input_line, tri_pos_temp);
                tri_pos.insert(tri_pos.end(), tri_pos_temp.begin(),tri_pos_temp.end());
                tri_pos_temp.clear();
                break;
        }
        // Output the text from the file
        // cout << "Content:" << endl;
        // cout << input_line << endl;
    }
    env.push_back(eye_pos);
    env.push_back(output_area);
    env.push_back(sphere_pos);
    env.push_back(tri_pos);
    /*
    for(int i=0;i<tri_pos.size();i++){
        cout << tri_pos[i] << endl;
    }
    */
    // cout << "Size or vector env: " << env.size() << endl;

    // Close the file
    InputFile.close();

    cout << "Ray Tracing starts..." << endl;
    Ray_tracing(256, 256, env);
    return 0;
}

