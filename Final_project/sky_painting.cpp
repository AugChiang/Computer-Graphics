#include "algebra3.h"
#include "generator.h"
#include<iostream>
#include<math.h>
#include<vector>
#include<fstream>
#include<sstream>
#include<string>
using namespace std;
#define PI 3.1415926

int main(void){
    // color constant
    const vec3 WHITE{250};
    const vec3 SKY_BLUE{20,100,200};
    const vec3 GRAY{80,80,100};
    const vec3 SUN_ORANGE{250,200,100};
    const vec3 INV_SUN{5,5,10}; // special color of inversed lighting
                                // it's kinda strange for this arg.

    // viewing constant
    const unsigned int WIDTH{800}; // column index
    const unsigned int HEIGHT{800}; // row index
    vec3 VIEW_v{0,-1,-1}; // viewing direction, starting from the eye to object
    const double VIEW_d{1000}; // viewing distance from canvas = sky base height
    const double VIEW_ang{120}; // viewing angle [degree]
    vec3 EYE_Pos{0,2000,0};
    const vec3 ny{0,1,0};
    vec3 SUN{1,1,0}; // point to the sun direction, starting from the object
    SUN = SUN.normalize();

    // clouds scaling to control size
    const double CLOUD_X_FCT{8e-3};
    const double CLOUD_Z_FCT{4e-3};
    const unsigned int GEN_N{8}; // the high frequency parameter

    // Canvas (screen) variables
    // Canvas plane = generated noise plane
    vec3 CANVAS_O{}; // Canvas original point position, at top-left
    vec3 CANVAS_vx{}; // viewer's right direction unit vector
    vec3 CANVAS_vz{}; // viewer's downward direction unit vector

    // '^': cross product
    CANVAS_vx = (VIEW_v^ny).normalize(); 
    CANVAS_vz = (VIEW_v^CANVAS_vx).normalize();
    // cout << "Canvas vx: " << CANVAS_vx[0] <<", " << CANVAS_vx[2] << endl;
    // cout << "Canvas vz: " << CANVAS_vz[0] <<", " << CANVAS_vz[2] << endl;

    // half distance of viewing projection on the ground(cloud).
    double half_gw = tan((VIEW_ang/2)*(PI/180));
    double half_gh = half_gw*(HEIGHT/WIDTH);
    double view_t = VIEW_d / VIEW_v.length(); // coefficient
    vec3 canvas_c = EYE_Pos + view_t*VIEW_v; // Canvas center position
    CANVAS_O = {canvas_c + half_gw*(-1)*CANVAS_vx + half_gh*(-1)*CANVAS_vz}; // top-left

    // Create and open a text file
    ofstream Output_File("output.ppm");
    ofstream Output_y("output_y.txt");
    Output_File << "P3" << endl; // for .ppm file format
    Output_File << HEIGHT <<" "<< WIDTH << endl;
    Output_File << "255" << endl;
    for(size_t i{0}; i<HEIGHT; i++){
        for(size_t j{0}; j<WIDTH; j++){
            vec3 pixel = CANVAS_O +
                         i*CANVAS_vx*CLOUD_X_FCT +
                         j*CANVAS_vz*CLOUD_Z_FCT; // pixel location on the canvas
                                                  // also, adjust the scaling.
            // cout << "Pixel Pos: " << pixel[0] <<", " << pixel[2] << endl;
            vec3 ray = (pixel - EYE_Pos).normalize(); // ray vector from eye to pixel
            vec3 pixel_color{GRAY}; // pixel color, default: sky color
            vec3 normal_gen{0,1,0}; // normal vector at position x,z on ground(cloud)
            vec3 sun_reflect{}; // reflection ray of sun lighting

            double dot_n_sun{}; // normal of cloud dot sun vector
            double dot_ray_view{}; // ray vector dot sun vector
            double dot_sun_view{}; // sun vector dot viewing vector

            vector<double> y_dx_dz = Gen_detail(pixel, GEN_N);
            double y = y_dx_dz[0];
            int depth_coef = y/2048; // for ambient color
            pixel[1] = y;
            Output_y << y << " "; // write to txt file.

            normal_gen[0] = -y_dx_dz[1]; // -dx component
            normal_gen[2] = -y_dx_dz[2]; // -dz component
            normal_gen = normal_gen.normalize();

            dot_n_sun = normal_gen*SUN;
            dot_ray_view = ray*VIEW_v;
            dot_sun_view = SUN*VIEW_v;
            
            pixel_color = 0.4*(dot_n_sun*WHITE)*Shader(pixel, SUN_ORANGE) + 
                          0.3*pow(dot_ray_view, 2)*SKY_BLUE +
                          0.2*(1-dot_sun_view)*(SUN_ORANGE) +
                          0.05*Reflect2Eye(normal_gen, SUN, -VIEW_v, INV_SUN);
            
            // write RGB value to the ppm file
            int r = round(pixel_color[0]);
            int g = round(pixel_color[1]);
            int b = round(pixel_color[2]);
            if(r<20){ r = 80*depth_coef; }
            if(r>255){ r = 255; }
            if(g<20){ g = 100*depth_coef; }
            if(g>255){ g = 255; }
            if(b<20){ b = 200*depth_coef; }
            if(b>255){ b = 255; }
            Output_File << r << " " << g << " " << b << " \n";
        }
        Output_y <<"\n";
    }

    cout << "File has been written." << endl;

    // Close the file
    Output_File.close();
    Output_y.close();

    return 0;
}