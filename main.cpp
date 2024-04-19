#include<bits/stdc++.h>
#include "geometry.h"
#include "model.h"
#include "tgaimage.h"
#include "omp.h"

using namespace std;

const TGAColor white = TGAColor{255, 255, 255, 255};
const TGAColor red   = TGAColor{0,  0,   255,255};//很奇怪,不是RGBa吗？
const int width=1000;
const int height=1000;
const int depth=255;
const double eps = 1e-18;
const double PI = acos(-1.0);
vector<double>zbuffer(width*height,-numeric_limits<double>::max());
vec3 camera{1,1,5};
vec3 center{0,0,0};
vec3 up{0,1,0};
vec3 light{1,1,-3};
double now_intensity=0;
struct Vec2i{
    int x,y;
};
void line(int x0,int y0,int x1,int y1,TGAImage &image ,TGAColor color){
    bool steep = false;
    if(abs(x0-x1)<abs(y0-y1)){
        steep =true;
        swap(x0,y0);
        swap(x1,y1);
    }
    if(x0>x1){
        swap(x0,x1);
        swap(y0,y1);
    }
    int dx = x1-x0;
    int dy = abs(y1-y0);
    int y=y0;
    int pk=2*dy-dx;
    for(int x=x0;x<=x1;x++){
        if(steep){
            image.set(y,x,color);
        }else{
            image.set(x,y,color);
        }
        //calc next point;
        {
            if(pk<0){
                pk+=2*dy;
            }else{
                pk+=2*dy-2*dx;
                y += (y1>y0) ? 1 : -1;
            }
        }
    }

}
mat<4,4> translation(vec3 v) {
	mat<4,4> Tr = mat<4,4>::identity();
	Tr[0][3] = v.x;
	Tr[1][3] = v.y;
	Tr[2][3] = v.z;
	return Tr;
}
 
mat<4,4> scale(double factorX, double factorY, double factorZ)
{
	mat<4,4> Z = mat<4,4>::identity();
	Z[0][0] = factorX;
	Z[1][1] = factorY;
	Z[2][2] = factorZ;
	return Z;
}
 
mat<4,4> rotation_x(double angle)
{
	angle = angle * PI / 180;
	double sinangle = sin(angle);
	double cosangle = cos(angle);
 
	mat<4,4> R = mat<4,4>::identity();
	R[1][1] = R[2][2] = cosangle;
	R[1][2] = -sinangle;
	R[2][1] = sinangle;
	return R;
}
 
mat<4,4> rotation_y(double angle)
{
	angle = angle * PI / 180;
	double sinangle = sin(angle);
	double cosangle = cos(angle);
 
	mat<4,4> R = mat<4,4>::identity();
	R[0][0] = R[2][2] = cosangle;
	R[0][2] = sinangle;
	R[2][0] = -sinangle;
	return R;
}
 
mat<4,4> rotation_z(double angle) {
	angle = angle * PI / 180;
	double sinangle = sin(angle);
	double cosangle = cos(angle);
 
	mat<4,4> R = mat<4,4>::identity();
	R[0][0] = R[1][1] = cosangle;
	R[0][1] = -sinangle;
	R[1][0] = sinangle;
	return R;
}
mat<4,4> Projection(vec3 eye,vec3 center){
    mat<4,4> res= mat<4,4>::identity();
    res[3][2] = -1.f/(camera-center).norm();
    return res;
}
mat<4,4> Lookat(vec3 eye, vec3 center, vec3 up) {
    mat<4,4> res = mat<4,4>::identity();
    vec3 z = (eye - center).normalized();
	vec3 x = cross(up,z).normalized();
	vec3 y = cross(z,x).normalized();
	for (int i = 0; i < 3; i++) {
		res[0][i] = x[i];
		res[1][i] = y[i];
		res[2][i] = z[i];
		res[i][3] = -center[i];
	}
	return res;
}
mat<4,4> Viewport(int x, int y, int w, int h) {
    mat<4,4> m = mat<4,4>::identity();
    m[0][3] = x+w/2.f;
    m[1][3] = y+h/2.f;
    m[2][3] = depth/2.f;

    m[0][0] = w/2.f;
    m[1][1] = h/2.f;
    m[2][2] = depth/2.f;
    return m;
}
void wired_model(Model &model,TGAImage &image){
    for(int i=0;i<model.nfaces();i++){
        //vert(nface ,nidx)
        for (int j=0; j<3; j++) { 
            vec3 v0 = model.vert(i,j); 
            vec3 v1 = model.vert(i,(j+1)%3); 
            int x0 = (v0.x+1.)*width/2.; 
            int y0 = (v0.y+1.)*height/2.; 
            int x1 = (v1.x+1.)*width/2.; 
            int y1 = (v1.y+1.)*height/2.; 
            line(x0, y0, x1, y1, image, white); 
        } 
    }
}
mat<4,4> projection_matrix = Projection(camera,center);
mat<4,4> viewport_matrix = Viewport(width/8, height/8, width*3/4, height*3/4);
mat<4,4> modelview_matrix = Lookat(camera,center,up);
vec3 barycentric(vec3 *pts,vec2 p){
    vec3 ans;

    vec3 x={
        pts[1][0]-pts[0][0],
        pts[2][0]-pts[0][0],
        pts[0][0]-p[0]
    };
    vec3 y={
        pts[1][1]-pts[0][1],
        pts[2][1]-pts[0][1],
        pts[0][1]-p[1]
    };
    ans = cross(x,y);
    if(ans[2]<1)return vec3{-1,1,1};
    return vec3{
        ans[0]/ans[2],
        ans[1]/ans[2],
        1-(ans[0]+ans[1])/ans[2]
    };
}
TGAImage now_texture;
void triangle(vec3 *pts,TGAImage &image,vector<vec2> pts_uv,vector<double> light_intensity,vector<double> &zbuffer){
    int boxmin[2](image.width()-1,image.height()-1);
    int boxmax[2](0,0);
    for(int i:{0,1,2}){
        boxmin[0] = max(0,min(boxmin[0],(int)pts[i][0]));
        boxmin[1] = max(0,min(boxmin[1],(int)pts[i][1]));
        boxmax[0] = min(image.width()-1,max(boxmax[0],(int)pts[i][0]));
        boxmax[1] = min(image.height()-1,max(boxmax[1],(int)pts[i][1]));
    }
   
// #pragma omp parallel for collapse(2) schedule(guided)
    for(int x=boxmin[0];x<=boxmax[0];x++){
        for(int y=boxmin[1];y<=boxmax[1];y++){ 
            vec2 p{
                (double)x,
                (double)y
                };
            vec3 weights=barycentric(pts,p);
            if(weights[0]<0||weights[1]<0||weights[2]<0)continue;
            else{
                double z=0;
                for(int i:{0,1,2}){
                    z+=weights[i]*pts[i][2];
                }
                int idx = p[0]*width+p[1];
                vec2 uv{0,0};
                double intensity = 0.f;
                for(int i:{0,1,2}){
                    uv = uv + weights[i]*pts_uv[i];
                    intensity += weights[i]*light_intensity[i];
                }
                
                intensity = max(5e-2,intensity);
                intensity = min(1.,intensity);
                // TGAColor wwhite = TGAColor{uint8_t(255*light_intensity),
                //                             uint8_t(255*light_intensity),
                //                             uint8_t(255*light_intensity),
                //                             uint8_t(255*light_intensity)};
                if(z>zbuffer[idx]){
                    TGAColor color = now_texture.get(uv.x*now_texture.width(),uv.y*now_texture.height());
                    for(int i:{0,1,2}){
                        color[i]*=intensity;
                    }
                    image.set(p[0],p[1],color);
                    zbuffer[idx]=z;
                }
            }
        }
    }
}
int cnt=0;
void tri_model(Model &model,TGAImage &image,vec3 light){
    
    for(int i=0;i<model.nfaces();i++){
        vec3 screen_coords[3];
        vec3 world_coords[3];
        vector<vec2>pts_uv(3);
        vector<double>light_intensity(3);
        now_texture=model.diffuse();
        for(int j:{0,1,2}){
            world_coords[j] = model.vert(i,j);
            vec3 v = model.vert(i,j);
            vec4 vv;
            vec2 uv = model.uv(i,j);
            vec3 vn = model.normal(i,j);
            
            vv = embed<4,3>(v,1);
            light_intensity[(j+2)%3] = vn.normalized()*(-1*light).normalized(); //光线方向别反了
            pts_uv[(j+2)%3]=uv;  //(bug)为什么会错位？1,2,0
            //This one needs a circular permutation of barycentric coordinates in the uv-interpolation:
            //对应这个错误  
            // int x = (v.x+1.)*width/2.; 
            // int y = (v.y+1.)*height/2.; 
            // screen_coords[j].x=x;
            // screen_coords[j].y=y;
            // screen_coords[j].z=world_coords[j].z;
            if(cnt==1){
                vv = translation({0,-0.005,-0.0181})*vv;
            }
            if(cnt==2){
                vv = translation({0,-0.005,-0.0367})*vv;
            }
            vec4 tmp = viewport_matrix * 
                       projection_matrix * 
                       modelview_matrix * 
                       translation({0,0,0}) *
                        vv;
            screen_coords[j]=embed<3,4>(tmp) / tmp[3];
        }
        vec3 normLine = cross(
            world_coords[2]-world_coords[0],
            world_coords[1]-world_coords[0]
        );
        normLine=normLine.normalized();
        double camera_view = normLine * (-1*camera).normalized();
        now_intensity = normLine * light.normalized();
        now_intensity = max(1e-1,now_intensity);
        if(
            // light_intensity > eps 不应该以光线删去面
            // && 
            camera_view > eps 
        )
            triangle(screen_coords,image,pts_uv,light_intensity,zbuffer);
    }
}
int main(int argc, char** argv) {
    std::cin.sync_with_stdio(0);
    std::cin.tie(0);
    std::cout.tie(0);
    vector<Model>models;
    TGAImage image(width,height,TGAImage::RGB);
    if(argc>=2){
        for(int i=1;i<argc;i++){
            Model model = Model(argv[i]);
            models.push_back(model);
        }
    }else{
        cout<<"No model input\n";
        return 0;
    }
    projection_matrix[3][2] = -1.f/(camera-center).norm();

    for(auto &model:models){
        auto start = omp_get_wtime();
        tri_model(model,image,light);
        // wired_model(model,image);
        auto end = omp_get_wtime();
        cout<<"time: "<<end-start<<" s"<<endl;
        cnt++;
    }
    
    image.write_tga_file("output.tga");
    return 0;
}