// C++ include
#include <iostream>
#include <string>
#include <vector>
#include <tbb/tbb.h>
#include <fstream>
// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"
#include "utils.h"
#include <math.h>
#include <Eigen/Dense>//行列式矩阵包含在Dense库中
#include <sys/time.h>  
// Shortcut to avoid Eigen:: and std:: everywhere, DO NOT USE IN .h
using namespace std;
using namespace Eigen;
using namespace tbb;
MatrixXd V1=MatrixXd::Zero(3,3);
MatrixXi F1=MatrixXi::Zero(3,3);
MatrixXd V2=MatrixXd::Zero(3,3);
MatrixXi F2=MatrixXi::Zero(3,3);
const Vector3d c_1(-1,-0.15,0);
const Vector3d c_2(1,-0.15,-0.3);
const double sphere_radius = 0.3;
//new
int x_len=0;
MatrixXd C = MatrixXd::Zero(800,800); // Store the color
MatrixXd C2 = MatrixXd::Zero(800,800); // Store the color
MatrixXd C3 = MatrixXd::Zero(800,800); // Store the color
MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask
Vector3d plane_origin(-1,1,1);
Vector3d x_displacement(2.0/C.cols(),0,0);
Vector3d y_displacement(0,-2.0/C.rows(),0);
Vector3d light_position(-2,2,3);//-100,100,30
Vector3d light_position2(2,2,-3);//3,2.2,10
Vector3d e(0,1,2);
const Vector3d surface_color1(1,0.6,0.3);
const Vector3d surface_color2(1,0.3,0.8);
const Vector3d surface_color3(0.3,1,0.8);
const Vector3d surface_color4(0.8,1,0.3);
MatrixXd Rotation(3,3);
double diffuse_render(Vector3d ray_L,Vector3d ray_normal)
{
     double temp=ray_L.transpose() * ray_normal;
     double diffuse = max(temp, 0.);
     return diffuse;
}
double light_render(Vector3d ray_L, Vector3d ray_V, Vector3d ray_normal)
{
    // Compute normal at the intersection point
    Vector3d ray_H = (ray_L+ray_V).normalized();//H向量
    // Simple diffuse model
    double temp=ray_L.transpose() * ray_normal;
    double diffuse = max(temp, 0.);
    double temp2=ray_normal.transpose()*ray_H;
    double Phong = pow(max(temp2, 0.),1000);
    return diffuse+Phong;
}
double if_ray_hit_triangle(Vector3d e, Vector3d d,Vector3d Va,Vector3d Vb,Vector3d Vc)
{
    //直线与三角形的交点
    MatrixXd Mb = MatrixXd::Zero(3,3);
    MatrixXd Mr = MatrixXd::Zero(3,3);
    MatrixXd Mt = MatrixXd::Zero(3,3);
    MatrixXd Ma = MatrixXd::Zero(3,3);
    double b,r,t;
    Mb<<Va(0)-e(0), Va(0)-Vc(0), d(0),
        Va(1)-e(1), Va(1)-Vc(1), d(1),
        Va(2)-e(2), Va(2)-Vc(2), d(2);
    Mr<<Va(0)-Vb(0), Va(0)-e(0), d(0),
        Va(1)-Vb(1), Va(1)-e(1), d(1),
        Va(2)-Vb(2), Va(2)-e(2), d(2);
    Mt<<Va(0)-Vb(0), Va(0)-Vc(0), Va(0)-e(0),
        Va(1)-Vb(1), Va(1)-Vc(1), Va(1)-e(1),
        Va(2)-Vb(2), Va(2)-Vc(2), Va(2)-e(2);
    Ma<<Va(0)-Vb(0), Va(0)-Vc(0), d(0),
        Va(1)-Vb(1), Va(1)-Vc(1), d(1),
        Va(2)-Vb(2), Va(2)-Vc(2), d(2);
    b=Mb.determinant()/Ma.determinant();
    r=Mr.determinant()/Ma.determinant();
    t=Mt.determinant()/Ma.determinant();
    if (r<0||r>1)
    {
        return 1000000;
    }
    else if (b<0||b>1-r)
    {
        return 1000000;
    }
    else
    {
        return t;
    }
    
}

double if_ray_hit_sphere(Vector3d ray_intersection, Vector3d ray_L,Vector3d center, double radius)
{
    double a=ray_L.transpose()*ray_L;
    double b_2=ray_L.transpose()*(ray_intersection-center);
    double c=(ray_intersection-center).transpose()*(ray_intersection-center)-radius*radius;
    if(b_2*b_2-a*c>=0)
    {
        double t_tmp1=(-b_2+sqrt(b_2*b_2-a*c))/a;
        double t_tmp2=(-b_2-sqrt(b_2*b_2-a*c))/a;
        double t;
        if(t_tmp1>0&&t_tmp2>0)
        {
            t=min(t_tmp1,t_tmp2);
        }
        else
        {
            t=max(t_tmp1,t_tmp2);
        }
        return t;
    }
    return 1000000;
}
double ray_hit_plane(Vector3d ray_intersection, Vector3d ray_L,double y)//射线起点，方向
{
    return (y-ray_intersection(1))/ray_L(1);
}
bool if_hit_object(Vector3d ray_intersection, Vector3d ray_L)
{
    //在光源1下，绘制所有其他物体对物体2造成的阴影
    for(int k=0;k<F1.rows();k++)//物体1
    {
        int tmp1,tmp2,tmp3;
        tmp1=F1(k,0);
        tmp2=F1(k,1);
        tmp3=F1(k,2);
        Vector3d Va,Vb,Vc;
        Va<<V1(tmp1,0),V1(tmp1,1),V1(tmp1,2);
        Vb<<V1(tmp2,0),V1(tmp2,1),V1(tmp2,2);
        Vc<<V1(tmp3,0),V1(tmp3,1),V1(tmp3,2);
        double tmp=if_ray_hit_triangle(ray_intersection,ray_L,Va,Vb,Vc);
        if(tmp>0.01&&tmp<1)//tmp=1000000表示无穷大即没有撞到，在这里线段有长度，故大于1就没有阻挡，等于零表示自身
        {
            return 1;
        }
    }
    for(int k=0;k<F2.rows();k++)//物体2
    {
        int tmp1,tmp2,tmp3;
        tmp1=F2(k,0);
        tmp2=F2(k,1);
        tmp3=F2(k,2);
        Vector3d Va,Vb,Vc;
        Va<<V2(tmp1,0),V2(tmp1,1),V2(tmp1,2);
        Vb<<V2(tmp2,0),V2(tmp2,1),V2(tmp2,2);
        Vc<<V2(tmp3,0),V2(tmp3,1),V2(tmp3,2);
        double tmp=if_ray_hit_triangle(ray_intersection,ray_L,Va,Vb,Vc);
        if(tmp>0.01&&tmp<1)//tmp=1000000表示无穷大即没有撞到，在这里线段有长度，故大于1就没有阻挡，等于零表示自身
        {
            return 1;
        }
    }
    double tt=if_ray_hit_sphere(ray_intersection,ray_L,c_1,sphere_radius);
    if(tt>0.01&&tt<1)//球一
    {
        return 1;
    }
    tt=if_ray_hit_sphere(ray_intersection,ray_L,c_2,sphere_radius);
    if(tt>0.01&&tt<1)//球二
    {
        return 1;
    }
    return 0;
}
double caculate_intersection(Vector3d e, Vector3d d, Vector3d&ray_normal,int&flag)
{
    double t0=0;
    double t1=1000000;
    double t;
    Vector3d Va;
    Vector3d Vb;
    Vector3d Vc;
    //对于第一个兔子的每一个面都要计算交点，然后用t更新t1
    for(int k=0;k<F1.rows();k++)
    {
        int tmp1,tmp2,tmp3;
        tmp1=F1(k,0);
        tmp2=F1(k,1);
        tmp3=F1(k,2);
        Va<<V1(tmp1,0),V1(tmp1,1),V1(tmp1,2);
        Vb<<V1(tmp2,0),V1(tmp2,1),V1(tmp2,2);
        Vc<<V1(tmp3,0),V1(tmp3,1),V1(tmp3,2);
        t=if_ray_hit_triangle(e,d,Va,Vb,Vc);
        if(t<t0||t>=t1)
        {
            continue;
        }
        else
        {
            Vector3d ab=Vb-Va;
            Vector3d ac=Vc-Va;
            ray_normal = ab.cross(ac).normalized();//法向量
            if(ray_normal.transpose()*d>0)
            {
                ray_normal=-ray_normal;
            }
            t1=t;
            flag=0;
        }
    }
    //对于第二个物体
    for(int k=0;k<F2.rows();k++)
    {
        int tmp1,tmp2,tmp3;
        tmp1=F2(k,0);
        tmp2=F2(k,1);
        tmp3=F2(k,2);
        Va<<V2(tmp1,0),V2(tmp1,1),V2(tmp1,2);
        Vb<<V2(tmp2,0),V2(tmp2,1),V2(tmp2,2);
        Vc<<V2(tmp3,0),V2(tmp3,1),V2(tmp3,2);
        t=if_ray_hit_triangle(e,d,Va,Vb,Vc);
        if(t<t0||t>=t1)
        {
            continue;
        }
        else
        {
            Vector3d ab=Vb-Va;
            Vector3d ac=Vc-Va;
            ray_normal = ab.cross(ac).normalized();//法向量
            if(ray_normal.transpose()*d>0)
            {
                ray_normal=-ray_normal;
            }
            t1=t;
            flag=1;
        }
    }
    //碰撞小球1
    t=if_ray_hit_sphere(e,d,c_1,sphere_radius);//返回值为t的大小，不相撞则t=1000000
    if(t<t1&&t>t0)
    {
        t1=t;
        flag=2;
        Vector3d ray_intersection=e+t1*d;
        ray_normal = (ray_intersection-c_1).normalized();
    }
    //碰撞小球2
    t=if_ray_hit_sphere(e,d,c_2,sphere_radius);
    if(t<t1&&t>t0)
    {
        t1=t;
        flag=3;
        Vector3d ray_intersection=e+t1*d;
        ray_normal = (ray_intersection-c_2).normalized();
    }
    return t1;
}
double caculate_intersection2(Vector3d e, Vector3d d, Vector3d&ray_normal,int&flag)
{
    double t0=0;
    double t1=1000000;
    double t;
    Vector3d Va;
    Vector3d Vb;
    Vector3d Vc;
    //对于第一个兔子的每一个面都要计算交点，然后用t更新t1
    for(int k=0;k<F1.rows();k++)
    {
        int tmp1,tmp2,tmp3;
        tmp1=F1(k,0);
        tmp2=F1(k,1);
        tmp3=F1(k,2);
        Va<<V1(tmp1,0),V1(tmp1,1),V1(tmp1,2);
        Vb<<V1(tmp2,0),V1(tmp2,1),V1(tmp2,2);
        Vc<<V1(tmp3,0),V1(tmp3,1),V1(tmp3,2);
        t=if_ray_hit_triangle(e,d,Va,Vb,Vc);
        if(t<t0||t>=t1)
        {
            continue;
        }
        else
        {
            Vector3d ab=Vb-Va;
            Vector3d ac=Vc-Va;
            ray_normal = ab.cross(ac).normalized();//法向量
            if(ray_normal.transpose()*d>0)
            {
                ray_normal=-ray_normal;
            }
            t1=t;
            flag=0;
        }
    }
    //对于第二个物体
    for(int k=0;k<F2.rows();k++)
    {
        int tmp1,tmp2,tmp3;
        tmp1=F2(k,0);
        tmp2=F2(k,1);
        tmp3=F2(k,2);
        Va<<V2(tmp1,0),V2(tmp1,1),V2(tmp1,2);
        Vb<<V2(tmp2,0),V2(tmp2,1),V2(tmp2,2);
        Vc<<V2(tmp3,0),V2(tmp3,1),V2(tmp3,2);
        t=if_ray_hit_triangle(e,d,Va,Vb,Vc);
        if(t<t0||t>=t1)
        {
            continue;
        }
        else
        {
            Vector3d ab=Vb-Va;
            Vector3d ac=Vc-Va;
            ray_normal = ab.cross(ac).normalized();//法向量
            if(ray_normal.transpose()*d>0)
            {
                ray_normal=-ray_normal;
            }
            t1=t;
            flag=1;
        }
    }
    //碰撞小球1
    t=if_ray_hit_sphere(e,d,c_1,sphere_radius);//返回值为t的大小，不相撞则t=1000000
    if(t<t1&&t>t0)
    {
        t1=t;
        flag=2;
        Vector3d ray_intersection=e+t1*d;
        ray_normal = (ray_intersection-c_1).normalized();
    }
    //碰撞小球2
    t=if_ray_hit_sphere(e,d,c_2,sphere_radius);
    if(t<t1&&t>t0)
    {
        t1=t;
        flag=3;
        Vector3d ray_intersection=e+t1*d;
        ray_normal = (ray_intersection-c_2).normalized();
    }
    //碰撞地面
    t=ray_hit_plane(e,d,-0.45);
    if(t>1000000)
    {
        //t=999999;//太远了就设置比无穷大小一点点
    }
    if(t<t1&&t>t0)
    {
        t1=t;
        flag=4;
        ray_normal<<0,1,0;
    }
    return t1;
}


void read_mesh()
{
    int nverts = 0;//顶点个数
    int nfaces = 0;//面个数
    int nedges = 0;//边个数
    int line_count = 1;//读入行数
    MatrixXd scale1(3,3);
    scale1<<10,0,0,
            0,10,0,
            0,0,10;
    Vector3d translation1(1,-0.78,-1);
    MatrixXd scale2(3,3);
    scale2<<0.2,0,0,
            0,0.2,0,
            0,0,0.2;
    Vector3d translation2(-1,0.43,-1);
    //分别读取两个文件
    string file1_path="/home/invisible/Hubs/CS-GY-6533/Assignment_1/data/bunny.off";
    string file2_path="/home/invisible/Hubs/CS-GY-6533/Assignment_1/data/bumpy_cube.off";
    ifstream file1(file1_path.c_str());
    ifstream file2(file2_path.c_str());
	string line;
	if(!file1||!file2) // 有该文件
	{
        cout <<"no such file" << endl;
        return;
	}
    //对于第一个
	while (getline (file1, line)) // line中不包括每行的换行符
    { 
        if(line_count==2)
        {
            int pos=0;
            pos=line.find(" ",0);
            nverts=atoi(line.substr(0,pos).c_str());
            line=line.substr(pos+1,-1);
            pos=line.find(" ",0);
            nfaces=atoi(line.substr(0,pos).c_str());
            nedges=atoi(line.substr(pos+1,-1).c_str());
            line_count++;
            break;
        } 
        line_count++;
        
    }
    V1.resize(nverts,3);
    F1.resize(nfaces,3);

    while (getline (file1, line))
    {
        if(line_count<nverts+3)//505
        {
            int pos=0;
            pos=line.find(" ",0);
            V1(line_count-3,0)=atof(line.substr(0,pos).c_str());
            line=line.substr(pos+1,-1);
            pos=line.find(" ",0);
            V1(line_count-3,1)=atof(line.substr(0,pos).c_str());
            V1(line_count-3,2)=atof(line.substr(pos+1,-1).c_str());
        }
        else if(line_count<nverts+3+nfaces)//1505
        {
            int pos=0;
            pos=line.find(" ",0);
            line=line.substr(pos+1,-1);


            pos=line.find(" ",0);
            F1(line_count-3-nverts,0)=atoi(line.substr(0,pos).c_str());
            line=line.substr(pos+1,-1);
            pos=line.find(" ",0);
            F1(line_count-3-nverts,1)=atoi(line.substr(0,pos).c_str());
            F1(line_count-3-nverts,2)=atoi(line.substr(pos+1,-1).c_str());
        }
        line_count++;
    }
    V1=(scale1*(V1.transpose())).transpose();
    for(int i=0;i<V1.rows();i++)
    {
        V1(i,0)=V1(i,0)+translation1(0);
        V1(i,1)=V1(i,1)+translation1(1);
        V1(i,2)=V1(i,2)+translation1(2);
    }
    //对于第二个
    line_count=1;
	while (getline (file2, line)) // line中不包括每行的换行符
    { 
        if(line_count==2)
        {
            int pos=0;
            pos=line.find(" ",0);
            nverts=atoi(line.substr(0,pos).c_str());
            line=line.substr(pos+1,-1);
            pos=line.find(" ",0);
            nfaces=atoi(line.substr(0,pos).c_str());
            nedges=atoi(line.substr(pos+1,-1).c_str());
            line_count++;
            break;
        } 
        line_count++;
        
    }
    V2.resize(nverts,3);
    F2.resize(nfaces,3);
    while (getline (file2, line))
    {
        if(line_count<nverts+3)//505
        {
            int pos=0;
            pos=line.find(" ",0);
            V2(line_count-3,0)=atof(line.substr(0,pos).c_str());
            line=line.substr(pos+1,-1);
            pos=line.find(" ",0);
            V2(line_count-3,1)=atof(line.substr(0,pos).c_str());
            V2(line_count-3,2)=atof(line.substr(pos+1,-1).c_str());
        }
        else if(line_count<nverts+3+nfaces)//1505
        {
            int pos=0;
            pos=line.find(" ",0);
            line=line.substr(pos+1,-1);


            pos=line.find(" ",0);
            F2(line_count-3-nverts,0)=atoi(line.substr(0,pos).c_str());
            line=line.substr(pos+1,-1);
            pos=line.find(" ",0);
            F2(line_count-3-nverts,1)=atoi(line.substr(0,pos).c_str());
            F2(line_count-3-nverts,2)=atoi(line.substr(pos+1,-1).c_str());
        }
        line_count++;
    }
    V2=(scale2*(V2.transpose())).transpose();
    for(int i=0;i<V2.rows();i++)
    {
        V2(i,0)=V2(i,0)+translation2(0);
        V2(i,1)=V2(i,1)+translation2(1);
        V2(i,2)=V2(i,2)+translation2(2);
    }
}

void part1()
{
    std::cout << "Part 1.1: Simple ray tracer, mutil-sphere with orthographic projection" << std::endl;

    const std::string filename("part1.png");
    MatrixXd C = MatrixXd::Zero(800,800); // Store the color
    MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask

    Vector3d origin(0,4,1);
    Vector3d x_displacement(4.0/C.cols(),0,0);
    Vector3d y_displacement(0,-4.0/C.rows(),0);

    // Single light source
    const Vector3d light_position(-100,100,30);
    //two center
    const Vector3d sphere_1(1,2.2,0);
    const Vector3d sphere_2(3,2.2,-0.3);


    for (unsigned i=0;i<C.cols();i++)
    {
        for (unsigned j=0;j<C.rows();j++)
        {
            // Prepare the ray
            Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
            //Vector3d ray_direction = RowVector3d(0,0,-1);//not used

            // Intersect with the sphere
            // NOTE: this is a special case of a sphere centered in the origin and for orthographic rays aligned with the z axis
            Vector2d ray_on_xy(ray_origin(0),ray_origin(1));
            const double sphere_radius = 0.9;
            //hit sphere 1
            double length_square=(sphere_1(0)-ray_on_xy(0))
                        *(sphere_1(0)-ray_on_xy(0))
                        +(sphere_1(1)-ray_on_xy(1))
                        *(sphere_1(1)-ray_on_xy(1));
            //hit sphere2
            double length_square2=(sphere_2(0)-ray_on_xy(0))
                        *(sphere_2(0)-ray_on_xy(0))
                        +(sphere_2(1)-ray_on_xy(1))
                        *(sphere_2(1)-ray_on_xy(1));
            
            if(length_square < sphere_radius*sphere_radius)
            {
                // The ray hit the sphere, compute the exact intersection point
                Vector3d ray_intersection(ray_on_xy(0),ray_on_xy(1),sqrt(sphere_radius*sphere_radius - length_square));
                // Compute normal at the intersection point
                Vector3d ray_normal = (ray_intersection-sphere_1).normalized();
                // Simple diffuse model
                C(i,j) = (light_position-ray_intersection).normalized().transpose() * ray_normal;
                // Clamp to zero
                C(i,j) = 0.8*max(C(i,j),0.);

                // Disable the alpha mask for this pixel
                A(i,j) = 1;
            }
            else if (length_square2 < sphere_radius*sphere_radius)
            {
                // The ray hit the sphere, compute the exact intersection point
                Vector3d ray_intersection2(ray_on_xy(0),ray_on_xy(1),sqrt(sphere_radius*sphere_radius - length_square2));
                // Compute normal at the intersection point
                Vector3d ray_normal2 = (ray_intersection2-sphere_2).normalized();
                // Simple diffuse model
                C(i,j) = (light_position-ray_intersection2).normalized().transpose() * ray_normal2;
                // Clamp to zero
                C(i,j) = 0.8*max(C(i,j),0.);

                // Disable the alpha mask for this pixel
                A(i,j) = 1;
            }
        }
    }
    // Save to png
    write_matrix_to_png(C,C,C,A,filename);
}

void part2()
{
    std::cout << "Part 1.2: Shading" << std::endl;
    const std::string filename("part2.png");
    MatrixXd C = MatrixXd::Zero(800,800); // Store the color
    MatrixXd C2 = MatrixXd::Zero(800,800); // Store the color
    MatrixXd C3 = MatrixXd::Zero(800,800); // Store the color
    MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask
    Vector3d origin(-2,2,1);
    Vector3d x_displacement(4.0/C.cols(),0,0);
    Vector3d y_displacement(0,-4.0/C.rows(),0);
    // Single light source
    const Vector3d light_position(-100,100,30);
    const Vector3d light_position2(3,2.2,10);
    //两球圆心坐标
    const Vector3d sphere_1(-1,0.2,0);
    const Vector3d sphere_2(1,0.2,-0.3);
    const Vector3d ray_view(0,0,-1);

    for (unsigned i=0;i<C.cols();i++)
    {
        for (unsigned j=0;j<C.rows();j++)
        {
            // Prepare the ray
            Vector3d ray_origin = origin + double(i)*x_displacement + double(j)*y_displacement;
            //Vector3d ray_direction = RowVector3d(0,0,-1);//not used

            // Intersect with the sphere
            // NOTE: this is a special case of a sphere centered in the origin and for orthographic rays aligned with the z axis
            Vector2d ray_on_xy(ray_origin(0),ray_origin(1));

            const double sphere_radius = 0.9;
            //碰撞小球1
            double length_square=(sphere_1(0)-ray_on_xy(0))
                        *(sphere_1(0)-ray_on_xy(0))
                        +(sphere_1(1)-ray_on_xy(1))
                        *(sphere_1(1)-ray_on_xy(1));
            //碰撞小球2
            double length_square2=(sphere_2(0)-ray_on_xy(0))
                        *(sphere_2(0)-ray_on_xy(0))
                        +(sphere_2(1)-ray_on_xy(1))
                        *(sphere_2(1)-ray_on_xy(1));
            if(length_square < sphere_radius*sphere_radius)
            {
                // The ray hit the sphere, compute the exact intersection point
                Vector3d ray_intersection(ray_on_xy(0),ray_on_xy(1),sqrt(sphere_radius*sphere_radius - length_square));
                //光源1
                // Compute normal at the intersection point
                Vector3d ray_normal = (ray_intersection-sphere_1).normalized();//法向量
                Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
                // Simple diffuse model
                double temp=ray_L.transpose() * ray_normal;
                double diffuse = max(temp, 0.);
                //光源2
                // Compute normal at the intersection point
                Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量a
                // Simple diffuse model
                double temp3=ray_L2.transpose() * ray_normal;
                double diffuse2 = max(temp3, 0.);
                
                C(i,j)=(diffuse+diffuse2)*surface_color1(0);
                C2(i,j)=(diffuse+diffuse2)*surface_color1(1);
                C3(i,j)=(diffuse+diffuse2)*surface_color1(2);
                // Disable the alpha mask for this pixel
                A(i,j) = 1;
            }
            else if (length_square2 < sphere_radius*sphere_radius)
            {
                // The ray hit the sphere, compute the exact intersection point
                Vector3d ray_intersection2(ray_on_xy(0),ray_on_xy(1),sqrt(sphere_radius*sphere_radius - length_square2));
                //光源1
                // Compute normal at the intersection point
                Vector3d ray_normal2 = (ray_intersection2-sphere_2).normalized();//法向量
                Vector3d ray_L = (light_position-ray_intersection2).normalized();//L向量
                Vector3d ray_H = (ray_L-ray_view).normalized();//H向量/////////////////////////
                // Simple diffuse model
                double temp=ray_L.transpose() * ray_normal2;
                double diffuse = max(temp, 0.);
                double temp2=ray_normal2.transpose()*ray_H;
                double Phong = pow(max(temp2, 0.),1000);
                //光源2
                // Compute normal at the intersection point
                Vector3d ray_L2 = (light_position2-ray_intersection2).normalized();//L向量
                Vector3d ray_H2 = (ray_L2-ray_view).normalized();//H向量//////////////////////////
                // Simple diffuse model
                double temp3=ray_L2.transpose() * ray_normal2;
                double diffuse2 = max(temp3, 0.);
                double temp4=ray_normal2.transpose()*ray_H2;
                double Phong2 = pow(max(temp4, 0.),1000);
                C(i,j)=(diffuse+Phong+diffuse2+Phong2)*surface_color2(0);
                C2(i,j)=(diffuse+Phong+diffuse2+Phong2)*surface_color2(1);
                C3(i,j)=(diffuse+Phong+diffuse2+Phong2)*surface_color2(2);

                // Disable the alpha mask for this pixel
                A(i,j) = 1;
            }
        }
    }
    // Save to png
    write_matrix_to_png(C,C2,C3,A,filename);
}

void part3()
{
    std::cout << "Part 1.3: Perspective Projection" << std::endl;
    const std::string filename("part3.png");
    MatrixXd C = MatrixXd::Zero(800,800); // Store the color
    MatrixXd C2 = MatrixXd::Zero(800,800); // Store the color
    MatrixXd C3 = MatrixXd::Zero(800,800); // Store the color
    MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask
    Vector3d plane_origin(-1,1,1);
    Vector3d x_displacement(2.0/C.cols(),0,0);
    Vector3d y_displacement(0,-2.0/C.rows(),0);
    // Single light source
    const Vector3d light_position(-100,100,30);
    const Vector3d light_position2(3,2.2,10);
    //两球圆心坐标
    // const Vector3d c_1(-0.9,0,0);
    // const Vector3d c_2(0.9,0,-1);
    const Vector3d c_1(-1,0.2,0);
    const Vector3d c_2(1,0.2,-0.3);
    const Vector3d e(0,0,2);
    const double sphere_radius = 0.65;
    double t0=0;
    double t1=1000000;
    double a=0;
    double b_2=0;
    double c=0;
    double t=0;
    Vector3d ray_normal(0,0,0);
    Vector3d ray_intersection(0,0,0);
    int flag=-1;
    for (unsigned i=0;i<C.cols();i++)
    {
        for (unsigned j=0;j<C.rows();j++)
        {
            t0=0;
            t1=1000000;
            // Prepare the ray
            Vector3d s=plane_origin+double(i)*x_displacement + double(j)*y_displacement;
            Vector3d d=s-e;
            double t_tmp1=0;
            double t_tmp2=0;
            //碰撞小球1
            a=d.transpose()*d;
            b_2=d.transpose()*(e-c_1);
            c=(e-c_1).transpose()*(e-c_1)-sphere_radius*sphere_radius;
            if(b_2*b_2-a*c<0)
            {
            
            }
            else
            {
                t_tmp1=(-b_2+sqrt(b_2*b_2-a*c))/a;
                t_tmp2=(-b_2-sqrt(b_2*b_2-a*c))/a;
                if(t_tmp1>0&&t_tmp2>0)
                {
                    t=min(t_tmp1,t_tmp2);
                }
                else
                {
                    t=max(t_tmp1,t_tmp2);
                }
                if(t<t1&&t>t0)
                {
                    t1=t;
                    flag=0;
                    
                }
                
            }
            //碰撞小球2
            a=d.transpose()*d;
            b_2=d.transpose()*(e-c_2);
            c=(e-c_2).transpose()*(e-c_2)-sphere_radius*sphere_radius;
            if(b_2*b_2-a*c<0)
            {

            }
            else
            {
                t_tmp1=(-b_2+sqrt(b_2*b_2-a*c))/a;
                t_tmp2=(-b_2-sqrt(b_2*b_2-a*c))/a;
                if(t_tmp1>0&&t_tmp2>0)
                {
                    t=min(t_tmp1,t_tmp2);
                }
                else
                {
                    t=max(t_tmp1,t_tmp2);
                }
                if(t<t1&&t>t0)
                {
                    t1=t;
                    flag=1;

                }
            }
            if(t1==1000000)
            {
                //设置成灰色
                C(i,j) = 0.1;
                C2(i,j)=C(i,j);
                C3(i,j)=C(i,j);
                A(i,j) = 1;
            }
            else
            {
                // The ray hit the sphere, compute the exact intersection point
                if(flag==1)
                {
                    ray_intersection=e+t1*d;
                    ray_normal = (ray_intersection-c_2).normalized();//法向量
                    //光源1
                    // Compute normal at the intersection point
                    Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
                    Vector3d ray_H = (ray_L-d.normalized()).normalized();//H向量/////////////////////////
                    // Simple diffuse model
                    double temp=ray_L.transpose() * ray_normal;
                    double diffuse = max(temp, 0.);
                    double temp2=ray_normal.transpose()*ray_H;
                    double Phong = pow(max(temp2, 0.),1000);
                    //光源2
                    // Compute normal at the intersection point
                    Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量
                    Vector3d ray_H2 = (ray_L2-d.normalized()).normalized();//H向量//////////////////////////
                    // Simple diffuse model
                    double temp3=ray_L2.transpose() * ray_normal;
                    double diffuse2 = max(temp3, 0.);
                    double temp4=ray_normal.transpose()*ray_H2;
                    double Phong2 = pow(max(temp4, 0.),1000);
                    C(i,j)=(diffuse+Phong+diffuse2+Phong2)*surface_color2(0);
                    C2(i,j)=(diffuse+Phong+diffuse2+Phong2)*surface_color2(1);
                    C3(i,j)=(diffuse+Phong+diffuse2+Phong2)*surface_color2(2);
                    // Disable the alpha mask for this pixel
                    A(i,j) = 1;
                }
                else
                {
                    ray_intersection=e+t1*d;
                    ray_normal = (ray_intersection-c_1).normalized();//法向量
                    //光源1
                    // Compute normal at the intersection point
                    Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
                    Vector3d ray_H = (ray_L-d.normalized()).normalized();//H向量/////////////////////////
                    // Simple diffuse model
                    double temp=ray_L.transpose() * ray_normal;
                    double diffuse = max(temp, 0.);
                    double temp2=ray_normal.transpose()*ray_H;
                    double Phong = pow(max(temp2, 0.),1000);
                    //光源2
                    // Compute normal at the intersection point
                    Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量
                    Vector3d ray_H2 = (ray_L2-d.normalized()).normalized();//H向量//////////////////////////
                    // Simple diffuse model
                    double temp3=ray_L2.transpose() * ray_normal;
                    double diffuse2 = max(temp3, 0.);
                    double temp4=ray_normal.transpose()*ray_H2;
                    double Phong2 = pow(max(temp4, 0.),1000);
                    C(i,j)=(diffuse+diffuse2)*surface_color1(0);
                    C2(i,j)=(diffuse+diffuse2)*surface_color1(1);
                    C3(i,j)=(diffuse+diffuse2)*surface_color1(2);
                    // Disable the alpha mask for this pixel
                    A(i,j) = 1;
                }
            }
        }
    }
    // Save to png
    write_matrix_to_png(C,C2,C3,A,filename);
}

void part4()
{
    //read_mesh();
    cout << "Part 1.4: Ray Tracing Triangle Meshes" << std::endl;
    const std::string filename("part4.png");
    double t0=0;
    double t1=1000000;
    double t=0;
    Vector3d ray_normal(0,0,0);
    Vector3d ray_intersection(0,0,0);
    int flag=-1;
    Vector3d Va,Vb,Vc;
    for (unsigned i=0;i<C.cols();i++)
    {
        for (unsigned j=0;j<C.rows();j++)
        {
            if(j%50==0)
            {
                cout<<i<<" "<<j<<endl;
            }
            flag=-1;
            t0=0;
            t1=1000000;
            // Prepare the ray
            Vector3d s=plane_origin+double(i)*x_displacement + double(j)*y_displacement;
            Vector3d d=s-e;
            t1=caculate_intersection(e,d,ray_normal,flag);
            if(t1==1000000)
            {
                //无穷远设置成纯黑色
                C(i,j) = 0;
                C2(i,j)=C(i,j);
                C3(i,j)=C(i,j);
                A(i,j) = 1;
            }
            else
            {

                if(flag==0)
                {
                    C(i,j)=0.13;
                    C2(i,j)=0.13;
                    C3(i,j)=0.13;
                    ray_intersection=e+t1*d;
                    //光源1
                    Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
                    //光源2
                    Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量
                    C(i,j) =C(i,j)+(light_render(ray_L,-d.normalized(),ray_normal)+light_render(ray_L2,-d.normalized(),ray_normal))*surface_color3(0);
                    C2(i,j)=C2(i,j)+(light_render(ray_L,-d.normalized(),ray_normal)+light_render(ray_L2,-d.normalized(),ray_normal))*surface_color3(1);
                    C3(i,j)=C3(i,j)+(light_render(ray_L,-d.normalized(),ray_normal)+light_render(ray_L2,-d.normalized(),ray_normal))*surface_color3(2);
                    // Disable the alpha mask for this pixel
                    A(i,j) = 1;
                }
                else if(flag==1)
                {
                    C(i,j)=0.13;
                    C2(i,j)=0.13;
                    C3(i,j)=0.13;
                    ray_intersection=e+t1*d;
                    //光源1
                    Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
                    //光源2
                    Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量
                    C(i,j) =C(i,j)+(light_render(ray_L,-d.normalized(),ray_normal)+light_render(ray_L2,-d.normalized(),ray_normal))*surface_color4(0);
                    C2(i,j)=C2(i,j)+(light_render(ray_L,-d.normalized(),ray_normal)+light_render(ray_L2,-d.normalized(),ray_normal))*surface_color4(1);
                    C3(i,j)=C3(i,j)+(light_render(ray_L,-d.normalized(),ray_normal)+light_render(ray_L2,-d.normalized(),ray_normal))*surface_color4(2);
                    // Disable the alpha mask for this pixel
                    A(i,j) = 1;
                }
                else if(flag==2)
                {
                    C(i,j)=0.13;
                    C2(i,j)=0.13;
                    C3(i,j)=0.13;
                    ray_intersection=e+t1*d;
                    ray_normal = (ray_intersection-c_1).normalized();//法向量
                    //光源1
                    // Compute normal at the intersection point
                    Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
                    // Simple diffuse model
                    double temp=ray_L.transpose() * ray_normal;
                    double diffuse = max(temp, 0.);
                    //光源2
                    // Compute normal at the intersection point
                    Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量
                    // Simple diffuse model
                    double temp3=ray_L2.transpose() * ray_normal;
                    double diffuse2 = max(temp3, 0.);
                    C(i,j) =C(i,j)+(diffuse+diffuse2)*surface_color1(0);
                    C2(i,j)=C2(i,j)+(diffuse+diffuse2)*surface_color1(1);
                    C3(i,j)=C3(i,j)+(diffuse+diffuse2)*surface_color1(2);
                    // Disable the alpha mask for this pixel
                    A(i,j) = 1;
                }
                else if(flag==3)
                {
                    C(i,j)=0.13;
                    C2(i,j)=0.13;
                    C3(i,j)=0.13;
                    ray_intersection=e+t1*d;
                    ray_normal = (ray_intersection-c_2).normalized();//法向量
                    //光源1
                    // Compute normal at the intersection point
                    Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
                    //光源2
                    // Compute normal at the intersection point
                    Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量
                    // Simple diffuse model
                    C(i,j) =C(i,j)+(light_render(ray_L,-d.normalized(),ray_normal)+light_render(ray_L2,-d.normalized(),ray_normal))*surface_color2(0);
                    C2(i,j)=C2(i,j)+(light_render(ray_L,-d.normalized(),ray_normal)+light_render(ray_L2,-d.normalized(),ray_normal))*surface_color2(1);
                    C3(i,j)=C3(i,j)+(light_render(ray_L,-d.normalized(),ray_normal)+light_render(ray_L2,-d.normalized(),ray_normal))*surface_color2(2);
                    // Disable the alpha mask for this pixel
                    A(i,j) = 1;
                }

            }
        }
    }
    // Save to png
    write_matrix_to_png(C,C2,C3,A,filename);
}

void part5()
{
    //read_mesh();
    cout << "Part 1.5: Shadow" << std::endl;
    const std::string filename("part5.png");
    // MatrixXd C = MatrixXd::Zero(800,800); // Store the color
    // MatrixXd C2 = MatrixXd::Zero(800,800); // Store the color
    // MatrixXd C3 = MatrixXd::Zero(800,800); // Store the color
    // MatrixXd A = MatrixXd::Zero(800,800); // Store the alpha mask
    // Vector3d plane_origin(-1,1,1);
    // Vector3d x_displacement(2.0/C.cols(),0,0);
    // Vector3d y_displacement(0,-2.0/C.rows(),0);
    // const Vector3d light_position(-100,100,30);
    // const Vector3d light_position2(3,2.2,10);
    // const Vector3d e(0,0,2);
    double t0=0;
    double t1=1000000;
    double t=0;
    Vector3d ray_normal(0,0,0);
    Vector3d ray_intersection(0,0,0);

    int flag=-1;
    Vector3d Va,Vb,Vc;
    for (unsigned i=0;i<C.cols();i++)
    {
        for (unsigned j=0;j<C.rows();j++)
        {    
            if(j%50==0)
            {
                cout<<i<<" "<<j<<endl;
            }
            flag=-1;
            t0=0;
            t1=1000000;
            // Prepare the ray
            Vector3d s=plane_origin+double(i)*x_displacement + double(j)*y_displacement;
            Vector3d d=s-e;
            t1=caculate_intersection(e,d,ray_normal,flag);
            if(t1==1000000)
            {
                //无穷远设置成纯黑色
                C(i,j) = 0;
                C2(i,j)=C(i,j);
                C3(i,j)=C(i,j);
                A(i,j) = 1;
            }
            else
            {

                if(flag==0)
                {
                    C(i,j)=0.13;
                    C2(i,j)=0.13;
                    C3(i,j)=0.13;
                    ray_intersection=e+t1*d;
                    //光源1
                    Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
                    //光源2
                    Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量
                    if(if_hit_object(ray_intersection,ray_L)==0)
                    {
                        C2(i,j) =C2(i,j)+light_render(ray_L,-d.normalized(),ray_normal);
                    }
                    if(if_hit_object(ray_intersection,ray_L2)==0)
                    {
                        C2(i,j) =C2(i,j)+light_render(ray_L2,-d.normalized(),ray_normal);
                    }
                    C(i,j) = 0.3*C2(i,j);
                    C3(i,j)=0.8*C2(i,j);
                    // Disable the alpha mask for this pixel
                    A(i,j) = 1;
                }
                else if(flag==1)
                {
                    C(i,j)=0.13;
                    C2(i,j)=0.13;
                    C3(i,j)=0.13;
                    ray_intersection=e+t1*d;
                    //光源1
                    Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
                    //光源2
                    Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量
                    if(if_hit_object(ray_intersection,ray_L)==0)
                    {
                        C2(i,j) =C2(i,j)+light_render(ray_L,-d.normalized(),ray_normal);
                    }
                    if(if_hit_object(ray_intersection,ray_L2)==0)
                    {
                        C2(i,j) =C2(i,j)+light_render(ray_L2,-d.normalized(),ray_normal);
                    }
                    C(i,j) = 0.8*C2(i,j);
                    C3(i,j)=0.3*C(i,j);
                    // Disable the alpha mask for this pixel
                    A(i,j) = 1;
                }
                else if(flag==2)
                {
                    C(i,j)=0.13;
                    C2(i,j)=0.13;
                    C3(i,j)=0.13;
                    ray_intersection=e+t1*d;
                    ray_normal = (ray_intersection-c_1).normalized();//法向量
                    //光源1
                    // Compute normal at the intersection point
                    Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
                    //光源2
                    // Compute normal at the intersection point
                    Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量
                    if(if_hit_object(ray_intersection,ray_L)==0)
                    {
                        C(i,j) =C(i,j)+diffuse_render(ray_L,ray_normal);
                    }
                    if(if_hit_object(ray_intersection,ray_L2)==0)
                    {
                        C(i,j) =C(i,j)+diffuse_render(ray_L2,ray_normal);
                    }
                    C2(i,j)=0.6*C(i,j);
                    C3(i,j)=0.3*C(i,j);
                    // Disable the alpha mask for this pixel
                    A(i,j) = 1;
                }
                else if(flag==3)
                {
                    C(i,j)=0.13;
                    C2(i,j)=0.13;
                    C3(i,j)=0.13;
                    ray_intersection=e+t1*d;
                    ray_normal = (ray_intersection-c_2).normalized();//法向量
                    //光源1
                    // Compute normal at the intersection point
                    Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
                    //光源2
                    // Compute normal at the intersection point
                    Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量
                    if(if_hit_object(ray_intersection,ray_L)==0)
                    {
                        C(i,j) =C(i,j)+light_render(ray_L,-d.normalized(),ray_normal);
                    }
                    if(if_hit_object(ray_intersection,ray_L2)==0)
                    {
                        C(i,j) =C(i,j)+light_render(ray_L2,-d.normalized(),ray_normal);
                    }
                
                    C2(i,j)=0.3*C(i,j);
                    C3(i,j)=0.8*C(i,j);
                    // Disable the alpha mask for this pixel
                    A(i,j) = 1;
                }

            }
        }
    }
    // Save to png
    write_matrix_to_png(C,C2,C3,A,filename);
}

void part6()
{


    //read_mesh();
    cout << "Part 1.6: Floor" << std::endl;
    const std::string filename("part6.png");
    // MatrixXd C = MatrixXd::Zero(200,200); // Store the color
    // MatrixXd C2 = MatrixXd::Zero(200,200); // Store the color
    // MatrixXd C3 = MatrixXd::Zero(200,200); // Store the color
    // MatrixXd A = MatrixXd::Zero(200,200); // Store the alpha mask
    // Vector3d plane_origin(-1,1,1);
    // Vector3d x_displacement(2.0/C.cols(),0,0);
    // Vector3d y_displacement(0,-2.0/C.rows(),0);
    // const Vector3d light_position(-100,100,30);
    // const Vector3d light_position2(3,2.2,10);
    // const Vector3d light_position(1,1,1);
    // const Vector3d light_position2(3,2.2,1.3);
    //const Vector3d e(0,0,2);
    double t0=0;
    double t1=1000000;
    double t=0;
    Vector3d ray_normal(0,0,0);
    Vector3d ray_intersection(0,0,0);
    int flag=-1;
    Vector3d Va,Vb,Vc;
    struct timeval timeStart, timeEnd; 
    double runTime=0; 
    gettimeofday(&timeStart, NULL );
    for (unsigned i=0;i<C.cols();i++)
    {
        for (unsigned j=0;j<C.rows();j++)
        {    
            if(j%50==0)
            {
                cout<<i<<" "<<j<<endl;
            }
            flag=-1;
            t0=0;
            t1=1000000;
            // Prepare the ray
            Vector3d s=plane_origin+double(i)*x_displacement + double(j)*y_displacement;
            Vector3d d=s-e;
            t1=caculate_intersection2(e,d,ray_normal,flag);//算上地面
            if(t1==1000000)
            {
                //无穷远设置成纯黑色
                C(i,j) = 0;
                C2(i,j)=0;
                C3(i,j)=0;
                A(i,j) = 1;
            }
            else
            {

                if(flag==0)
                {
                    C(i,j)=0.13;
                    C2(i,j)=0.13;
                    C3(i,j)=0.13;
                    ray_intersection=e+t1*d;
                    //光源1
                    Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
                    //光源2
                    Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量
                    if(if_hit_object(ray_intersection,ray_L)==0)
                    {
                        C2(i,j) =C2(i,j)+light_render(ray_L,-d.normalized(),ray_normal);
                    }
                    if(if_hit_object(ray_intersection,ray_L2)==0)
                    {
                        C2(i,j) =C2(i,j)+light_render(ray_L2,-d.normalized(),ray_normal);
                    }
                    C(i,j) = 0.3*C2(i,j);
                    C3(i,j)=0.8*C2(i,j);
                    // Disable the alpha mask for this pixel
                    A(i,j) = 1;
                }
                else if(flag==1)
                {
                    C(i,j)=0.13;
                    C2(i,j)=0.13;
                    C3(i,j)=0.13;
                    ray_intersection=e+t1*d;
                    //光源1
                    Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
                    //光源2
                    Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量
                    if(if_hit_object(ray_intersection,ray_L)==0)
                    {
                        C2(i,j) =C2(i,j)+light_render(ray_L,-d.normalized(),ray_normal);
                    }
                    if(if_hit_object(ray_intersection,ray_L2)==0)
                    {
                        C2(i,j) =C2(i,j)+light_render(ray_L2,-d.normalized(),ray_normal);
                    }
                    C(i,j) = 0.8*C2(i,j);
                    C3(i,j)=0.3*C(i,j);
                    // Disable the alpha mask for this pixel
                    A(i,j) = 1;
                }
                else if(flag==2)
                {
                    C(i,j)=0.13;
                    C2(i,j)=0.13;
                    C3(i,j)=0.13;
                    ray_intersection=e+t1*d;
                    ray_normal = (ray_intersection-c_1).normalized();//法向量
                    //光源1
                    // Compute normal at the intersection point
                    Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
                    //光源2
                    // Compute normal at the intersection point
                    Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量
                    if(if_hit_object(ray_intersection,ray_L)==0)
                    {
                        C(i,j) =C(i,j)+diffuse_render(ray_L,ray_normal);
                    }
                    if(if_hit_object(ray_intersection,ray_L2)==0)
                    {
                        C(i,j) =C(i,j)+diffuse_render(ray_L2,ray_normal);
                    }
                    C2(i,j)=0.6*C(i,j);
                    C3(i,j)=0.3*C(i,j);
                    // Disable the alpha mask for this pixel
                    A(i,j) = 1;
                }
                else if(flag==3)
                {
                    C(i,j)=0.13;
                    C2(i,j)=0.13;
                    C3(i,j)=0.13;
                    ray_intersection=e+t1*d;
                    ray_normal = (ray_intersection-c_2).normalized();//法向量
                    //光源1
                    // Compute normal at the intersection point
                    Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
                    //光源2
                    // Compute normal at the intersection point
                    Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量
                    if(if_hit_object(ray_intersection,ray_L)==0)
                    {
                        C(i,j) =C(i,j)+light_render(ray_L,-d.normalized(),ray_normal);
                    }
                    if(if_hit_object(ray_intersection,ray_L2)==0)
                    {
                        C(i,j) =C(i,j)+light_render(ray_L2,-d.normalized(),ray_normal);
                    }
                
                    C2(i,j)=0.3*C(i,j);
                    C3(i,j)=0.8*C(i,j);
                    // Disable the alpha mask for this pixel
                    A(i,j) = 1;
                }
                else if(flag==4)
                {
                    C(i,j)=0.13;
                    C2(i,j)=0.13;
                    C3(i,j)=0.13;
                    ray_intersection=e+t1*d;
                    //光源1
                    Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
                    //光源2
                    Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量
                    if(if_hit_object(ray_intersection,ray_L)==0)
                    {
                        C(i,j) =C(i,j)+0.6*light_render(ray_L,-d.normalized(),ray_normal);
                        C2(i,j) = C(i,j);
                        C3(i,j)= C(i,j);
                    }
                    if(if_hit_object(ray_intersection,ray_L2)==0)
                    {
                        C(i,j) =C(i,j)+0.6*light_render(ray_L2,-d.normalized(),ray_normal);
                        C2(i,j) = C(i,j);
                        C3(i,j)= C(i,j);
                    }
                    //增加倒影
                    //反射光线，不论地面的某一点是否在阴影下都要计算反射光线
                    //计算反射光线与所有物体交点，取最近的那个
                    double coefficient=d.normalized().transpose()*ray_normal;
                    coefficient=2*coefficient;
                    Vector3d r=d.normalized()-coefficient*ray_normal;
                    Vector3d ray_normal2;
                    flag=-1;
                    t1=caculate_intersection(ray_intersection,r,ray_normal2,flag);
                    if(t1==1000000)//反射光没有和任何物体相撞
                    {
                        C(i,j) = C(i,j)+0;
                        C2(i,j)=C2(i,j)+0;
                        C3(i,j)=C3(i,j)+0;
                        A(i,j) = 1;
                    }
                    else
                    {
                        if(flag==0)
                        {
                            Vector3d ray_intersection2=ray_intersection+t1*r;
                            //光源1
                            Vector3d ray_L = (light_position-ray_intersection2).normalized();//L向量
                            //光源2
                            Vector3d ray_L2 = (light_position2-ray_intersection2).normalized();//L向量
                            if(if_hit_object(ray_intersection2,ray_L)==0)
                            {
                                C(i,j) =C(i,j)+0.6*0.3*light_render(ray_L,-r.normalized(),ray_normal2);
                                C2(i,j)=C2(i,j)+0.6*light_render(ray_L,-r.normalized(),ray_normal2);
                                C3(i,j)=C3(i,j)+0.6*0.8*light_render(ray_L,-r.normalized(),ray_normal2);
                            }
                            if(if_hit_object(ray_intersection2,ray_L2)==0)
                            {
                                C(i,j) =C(i,j)+0.6*0.3*light_render(ray_L2,-r.normalized(),ray_normal2);
                                C2(i,j)=C2(i,j)+0.6*light_render(ray_L2,-r.normalized(),ray_normal2);
                                C3(i,j)=C3(i,j)+0.6*0.8*light_render(ray_L2,-r.normalized(),ray_normal2);
                            }
                            
                        }
                        else if (flag==1)
                        {
                            Vector3d ray_intersection2=ray_intersection+t1*r;
                            //光源1
                            Vector3d ray_L = (light_position-ray_intersection2).normalized();//L向量
                            //光源2
                            Vector3d ray_L2 = (light_position2-ray_intersection2).normalized();//L向量
                            if(if_hit_object(ray_intersection2,ray_L)==0)
                            {
                                C(i,j) =C(i,j)+0.6*0.8*light_render(ray_L,-r.normalized(),ray_normal2);
                                C2(i,j)=C2(i,j)+0.6*light_render(ray_L,-r.normalized(),ray_normal2);
                                C3(i,j)=C3(i,j)+0.6*0.3*light_render(ray_L,-r.normalized(),ray_normal2);
                            }
                            if(if_hit_object(ray_intersection2,ray_L2)==0)
                            {
                                C(i,j) =C(i,j)+0.6*0.8*light_render(ray_L2,-r.normalized(),ray_normal2);
                                C2(i,j)=C2(i,j)+0.6*light_render(ray_L2,-r.normalized(),ray_normal2);
                                C3(i,j)=C3(i,j)+0.6*0.3*light_render(ray_L2,-r.normalized(),ray_normal2);
                            }
                        }
                        else if (flag==2)
                        {
                            Vector3d ray_intersection2=ray_intersection+t1*r;
                            //光源1
                            Vector3d ray_L = (light_position-ray_intersection2).normalized();//L向量
                            //光源2
                            Vector3d ray_L2 = (light_position2-ray_intersection2).normalized();//L向量
                            if(if_hit_object(ray_intersection2,ray_L)==0)
                            {
                                C(i,j) =C(i,j)+0.6*diffuse_render(ray_L,ray_normal2);
                                C2(i,j)=C2(i,j)+0.6*0.8*diffuse_render(ray_L,ray_normal2);
                                C3(i,j)=C3(i,j)+0.6*0.3*diffuse_render(ray_L,ray_normal2);
                            }
                            if(if_hit_object(ray_intersection2,ray_L2)==0)
                            {
                                C(i,j) =C(i,j)+0.6*diffuse_render(ray_L2,ray_normal2);
                                C2(i,j)=C2(i,j)+0.6*0.6*diffuse_render(ray_L2,ray_normal2);
                                C3(i,j)=C3(i,j)+0.6*0.3*diffuse_render(ray_L2,ray_normal2);
                            }
                        }
                        else if (flag==3)
                        {
                            Vector3d ray_intersection2=ray_intersection+t1*r;
                            //光源1
                            Vector3d ray_L = (light_position-ray_intersection2).normalized();//L向量
                            //光源2
                            Vector3d ray_L2 = (light_position2-ray_intersection2).normalized();//L向量
                            if(if_hit_object(ray_intersection2,ray_L)==0)
                            {
                                C(i,j) =C(i,j)+0.6*light_render(ray_L,-r.normalized(),ray_normal2);
                                C2(i,j)=C2(i,j)+0.6*0.3*light_render(ray_L,-r.normalized(),ray_normal2);
                                C3(i,j)=C3(i,j)+0.6*0.8*light_render(ray_L,-r.normalized(),ray_normal2);
                            }
                            if(if_hit_object(ray_intersection2,ray_L2)==0)
                            {
                                C(i,j) =C(i,j)+0.6*light_render(ray_L2,-r.normalized(),ray_normal2);
                                C2(i,j)=C2(i,j)+0.6*0.3*light_render(ray_L2,-r.normalized(),ray_normal2);
                                C3(i,j)=C3(i,j)+0.6*0.8*light_render(ray_L2,-r.normalized(),ray_normal2);
                            }
                        }
                    }
                    // Disable the alpha mask for this pixel
                    A(i,j) = 1;
                }

            }
        }
    }
    gettimeofday( &timeEnd, NULL ); 
    runTime = (timeEnd.tv_sec - timeStart.tv_sec ) + (double)(timeEnd.tv_usec -timeStart.tv_usec)/1000000;  
    printf("Part6 runTime is %lf\n", runTime); 
    // Save to png
    write_matrix_to_png(C,C2,C3,A,filename);
}
void draw_each_coloum(int j)
{
    if(j%50==0)
    {
        cout<<x_len<<" "<<j<<endl;
    }
    int flag=-1;
    double t0=0;
    double t1=1000000;
    double t=0;
    Vector3d ray_normal(0,0,0);
    Vector3d ray_intersection(0,0,0);
    Vector3d Va,Vb,Vc;
    // Prepare the ray
    Vector3d s=plane_origin+double(x_len)*x_displacement+double(j)*y_displacement;
    Vector3d d=s-e;
    t1=caculate_intersection2(e,d,ray_normal,flag);//算上地面
    //cout<<x_len<<"  "<<j<<"   "<<flag<<endl;
    if(t1==1000000)
    {
        //无穷远设置成纯黑色
        C(x_len,j) = 0;
        C2(x_len,j)=0;
        C3(x_len,j)=0;
        A(x_len,j) = 1;
    }
    else
    {
        
        if(flag==0)
        {
            C(x_len,j)=0.13;
            C2(x_len,j)=0.13;
            C3(x_len,j)=0.13;
            ray_intersection=e+t1*d;
            //光源1
            Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
            //光源2
            Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量
            if(if_hit_object(ray_intersection,ray_L)==0)
            {
                C2(x_len,j) =C2(x_len,j)+light_render(ray_L,-d.normalized(),ray_normal);
            }
            if(if_hit_object(ray_intersection,ray_L2)==0)
            {
                C2(x_len,j) =C2(x_len,j)+light_render(ray_L2,-d.normalized(),ray_normal);
            }
            C(x_len,j) = 0.3*C2(x_len,j);
            C3(x_len,j)=0.8*C2(x_len,j);
            // Disable the alpha mask for this pixel
            A(x_len,j) = 1;
        }
        else if(flag==1)
        {
            C(x_len,j)=0.13;
            C2(x_len,j)=0.13;
            C3(x_len,j)=0.13;
            ray_intersection=e+t1*d;
            //光源1
            Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
            //光源2
            Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量
            if(if_hit_object(ray_intersection,ray_L)==0)
            {
                C2(x_len,j) =C2(x_len,j)+light_render(ray_L,-d.normalized(),ray_normal);
            }
            if(if_hit_object(ray_intersection,ray_L2)==0)
            {
                C2(x_len,j) =C2(x_len,j)+light_render(ray_L2,-d.normalized(),ray_normal);
            }
            C(x_len,j) = 0.8*C2(x_len,j);
            C3(x_len,j)=0.3*C(x_len,j);
            // Disable the alpha mask for this pixel
            A(x_len,j) = 1;
        }
        else if(flag==2)
        {
            C(x_len,j)=0.13;
            C2(x_len,j)=0.13;
            C3(x_len,j)=0.13;
            ray_intersection=e+t1*d;
            ray_normal = (ray_intersection-c_1).normalized();//法向量
            //光源1
            // Compute normal at the intersection point
            Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
            //光源2
            // Compute normal at the intersection point
            Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量
            if(if_hit_object(ray_intersection,ray_L)==0)
            {
                C(x_len,j) =C(x_len,j)+diffuse_render(ray_L,ray_normal);
            }
            if(if_hit_object(ray_intersection,ray_L2)==0)
            {
                C(x_len,j) =C(x_len,j)+diffuse_render(ray_L2,ray_normal);
            }
            C2(x_len,j)=0.6*C(x_len,j);
            C3(x_len,j)=0.3*C(x_len,j);
            // Disable the alpha mask for this pixel
            A(x_len,j) = 1;
        }
        else if(flag==3)
        {
            C(x_len,j)=0.13;
            C2(x_len,j)=0.13;
            C3(x_len,j)=0.13;
            ray_intersection=e+t1*d;
            ray_normal = (ray_intersection-c_2).normalized();//法向量
            //光源1
            // Compute normal at the intersection point
            Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
            //光源2
            // Compute normal at the intersection point
            Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量
            if(if_hit_object(ray_intersection,ray_L)==0)
            {
                C(x_len,j) =C(x_len,j)+light_render(ray_L,-d.normalized(),ray_normal);
            }
            if(if_hit_object(ray_intersection,ray_L2)==0)
            {
                C(x_len,j) =C(x_len,j)+light_render(ray_L2,-d.normalized(),ray_normal);
            }
        
            C2(x_len,j)=0.3*C(x_len,j);
            C3(x_len,j)=0.8*C(x_len,j);
            // Disable the alpha mask for this pixel
            A(x_len,j) = 1;
        }
        else if(flag==4)
        {
            C(x_len,j)=0.13;
            C2(x_len,j)=0.13;
            C3(x_len,j)=0.13;
            ray_intersection=e+t1*d;
            //光源1
            Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
            //光源2
            Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量
            if(if_hit_object(ray_intersection,ray_L)==0)
            {
                C(x_len,j) =C(x_len,j)+0.6*light_render(ray_L,-d.normalized(),ray_normal);
                C2(x_len,j) = C(x_len,j);
                C3(x_len,j)= C(x_len,j);
            }
            if(if_hit_object(ray_intersection,ray_L2)==0)
            {
                C(x_len,j) =C(x_len,j)+0.6*light_render(ray_L2,-d.normalized(),ray_normal);
                C2(x_len,j) = C(x_len,j);
                C3(x_len,j)= C(x_len,j);
            }
            //增加倒影
            //反射光线，不论地面的某一点是否在阴影下都要计算反射光线
            //计算反射光线与所有物体交点，取最近的那个
            double coefficient=d.normalized().transpose()*ray_normal;
            coefficient=2*coefficient;
            Vector3d r=d.normalized()-coefficient*ray_normal;
            Vector3d ray_normal2;
            flag=-1;
            t1=caculate_intersection(ray_intersection,r,ray_normal2,flag);
            if(t1==1000000)//反射光没有和任何物体相撞
            {
                C(x_len,j) = C(x_len,j)+0;
                C2(x_len,j)=C2(x_len,j)+0;
                C3(x_len,j)=C3(x_len,j)+0;
                A(x_len,j) = 1;
            }
            else
            {
                if(flag==0)
                {
                    Vector3d ray_intersection2=ray_intersection+t1*r;
                    //光源1
                    Vector3d ray_L = (light_position-ray_intersection2).normalized();//L向量
                    //光源2
                    Vector3d ray_L2 = (light_position2-ray_intersection2).normalized();//L向量
                    if(if_hit_object(ray_intersection2,ray_L)==0)
                    {
                        C(x_len,j) =C(x_len,j)+0.6*0.3*light_render(ray_L,-r.normalized(),ray_normal2);
                        C2(x_len,j)=C2(x_len,j)+0.6*light_render(ray_L,-r.normalized(),ray_normal2);
                        C3(x_len,j)=C3(x_len,j)+0.6*0.8*light_render(ray_L,-r.normalized(),ray_normal2);
                    }
                    if(if_hit_object(ray_intersection2,ray_L2)==0)
                    {
                        C(x_len,j) =C(x_len,j)+0.6*0.3*light_render(ray_L2,-r.normalized(),ray_normal2);
                        C2(x_len,j)=C2(x_len,j)+0.6*light_render(ray_L2,-r.normalized(),ray_normal2);
                        C3(x_len,j)=C3(x_len,j)+0.6*0.8*light_render(ray_L2,-r.normalized(),ray_normal2);
                    }
                    
                }
                else if (flag==1)
                {
                    Vector3d ray_intersection2=ray_intersection+t1*r;
                    //光源1
                    Vector3d ray_L = (light_position-ray_intersection2).normalized();//L向量
                    //光源2
                    Vector3d ray_L2 = (light_position2-ray_intersection2).normalized();//L向量
                    if(if_hit_object(ray_intersection2,ray_L)==0)
                    {
                        C(x_len,j) =C(x_len,j)+0.6*0.8*light_render(ray_L,-r.normalized(),ray_normal2);
                        C2(x_len,j)=C2(x_len,j)+0.6*light_render(ray_L,-r.normalized(),ray_normal2);
                        C3(x_len,j)=C3(x_len,j)+0.6*0.3*light_render(ray_L,-r.normalized(),ray_normal2);
                    }
                    if(if_hit_object(ray_intersection2,ray_L2)==0)
                    {
                        C(x_len,j) =C(x_len,j)+0.6*0.8*light_render(ray_L2,-r.normalized(),ray_normal2);
                        C2(x_len,j)=C2(x_len,j)+0.6*light_render(ray_L2,-r.normalized(),ray_normal2);
                        C3(x_len,j)=C3(x_len,j)+0.6*0.3*light_render(ray_L2,-r.normalized(),ray_normal2);
                    }
                }
                else if (flag==2)
                {
                    Vector3d ray_intersection2=ray_intersection+t1*r;
                    //光源1
                    Vector3d ray_L = (light_position-ray_intersection2).normalized();//L向量
                    //光源2
                    Vector3d ray_L2 = (light_position2-ray_intersection2).normalized();//L向量
                    if(if_hit_object(ray_intersection2,ray_L)==0)
                    {
                        C(x_len,j) =C(x_len,j)+0.6*diffuse_render(ray_L,ray_normal2);
                        C2(x_len,j)=C2(x_len,j)+0.6*0.8*diffuse_render(ray_L,ray_normal2);
                        C3(x_len,j)=C3(x_len,j)+0.6*0.3*diffuse_render(ray_L,ray_normal2);
                    }
                    if(if_hit_object(ray_intersection2,ray_L2)==0)
                    {
                        C(x_len,j) =C(x_len,j)+0.6*diffuse_render(ray_L2,ray_normal2);
                        C2(x_len,j)=C2(x_len,j)+0.6*0.6*diffuse_render(ray_L2,ray_normal2);
                        C3(x_len,j)=C3(x_len,j)+0.6*0.3*diffuse_render(ray_L2,ray_normal2);
                    }
                }
                else if (flag==3)
                {
                    Vector3d ray_intersection2=ray_intersection+t1*r;
                    //光源1
                    Vector3d ray_L = (light_position-ray_intersection2).normalized();//L向量
                    //光源2
                    Vector3d ray_L2 = (light_position2-ray_intersection2).normalized();//L向量
                    if(if_hit_object(ray_intersection2,ray_L)==0)
                    {
                        C(x_len,j) =C(x_len,j)+0.6*light_render(ray_L,-r.normalized(),ray_normal2);
                        C2(x_len,j)=C2(x_len,j)+0.6*0.3*light_render(ray_L,-r.normalized(),ray_normal2);
                        C3(x_len,j)=C3(x_len,j)+0.6*0.8*light_render(ray_L,-r.normalized(),ray_normal2);
                    }
                    if(if_hit_object(ray_intersection2,ray_L2)==0)
                    {
                        C(x_len,j) =C(x_len,j)+0.6*light_render(ray_L2,-r.normalized(),ray_normal2);
                        C2(x_len,j)=C2(x_len,j)+0.6*0.3*light_render(ray_L2,-r.normalized(),ray_normal2);
                        C3(x_len,j)=C3(x_len,j)+0.6*0.8*light_render(ray_L2,-r.normalized(),ray_normal2);
                    }
                }
            }
            // Disable the alpha mask for this pixel
            A(x_len,j) = 1;
        }

    }
}
void draw_each_coloum2(int j)
{
    if(j%50==0)
    {
        cout<<x_len<<" "<<j<<endl;
    }
    int flag=-1;
    double t0=0;
    double t1=1000000;
    double t=0;
    Vector3d ray_normal(0,0,0);
    Vector3d ray_intersection(0,0,0);
    Vector3d Va,Vb,Vc;
    // Prepare the ray
    Vector3d s=plane_origin+double(x_len)*x_displacement+double(j)*y_displacement;
    s=Rotation*s;   
    Vector3d d=s-e;
    t1=caculate_intersection2(e,d,ray_normal,flag);//算上地面
    //cout<<x_len<<"  "<<j<<"   "<<flag<<endl;
    if(t1==1000000)
    {
        //无穷远设置成纯黑色
        C(x_len,j) = 0;
        C2(x_len,j)=0;
        C3(x_len,j)=0;
        A(x_len,j) = 1;
    }
    else
    {
        
        if(flag==0)
        {
            C(x_len,j)=0.13;
            C2(x_len,j)=0.13;
            C3(x_len,j)=0.13;
            ray_intersection=e+t1*d;
            //光源1
            Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
            //光源2
            Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量
            if(if_hit_object(ray_intersection,ray_L)==0)
            {
                C2(x_len,j) =C2(x_len,j)+light_render(ray_L,-d.normalized(),ray_normal);
            }
            if(if_hit_object(ray_intersection,ray_L2)==0)
            {
                C2(x_len,j) =C2(x_len,j)+light_render(ray_L2,-d.normalized(),ray_normal);
            }
            C(x_len,j) = 0.3*C2(x_len,j);
            C3(x_len,j)=0.8*C2(x_len,j);
            // Disable the alpha mask for this pixel
            A(x_len,j) = 1;
        }
        else if(flag==1)
        {
            C(x_len,j)=0.13;
            C2(x_len,j)=0.13;
            C3(x_len,j)=0.13;
            ray_intersection=e+t1*d;
            //光源1
            Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
            //光源2
            Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量
            if(if_hit_object(ray_intersection,ray_L)==0)
            {
                C2(x_len,j) =C2(x_len,j)+light_render(ray_L,-d.normalized(),ray_normal);
            }
            if(if_hit_object(ray_intersection,ray_L2)==0)
            {
                C2(x_len,j) =C2(x_len,j)+light_render(ray_L2,-d.normalized(),ray_normal);
            }
            C(x_len,j) = 0.8*C2(x_len,j);
            C3(x_len,j)=0.3*C(x_len,j);
            // Disable the alpha mask for this pixel
            A(x_len,j) = 1;
        }
        else if(flag==2)
        {
            C(x_len,j)=0.13;
            C2(x_len,j)=0.13;
            C3(x_len,j)=0.13;
            ray_intersection=e+t1*d;
            ray_normal = (ray_intersection-c_1).normalized();//法向量
            //光源1
            // Compute normal at the intersection point
            Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
            //光源2
            // Compute normal at the intersection point
            Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量
            if(if_hit_object(ray_intersection,ray_L)==0)
            {
                C(x_len,j) =C(x_len,j)+diffuse_render(ray_L,ray_normal);
            }
            if(if_hit_object(ray_intersection,ray_L2)==0)
            {
                C(x_len,j) =C(x_len,j)+diffuse_render(ray_L2,ray_normal);
            }
            C2(x_len,j)=0.6*C(x_len,j);
            C3(x_len,j)=0.3*C(x_len,j);
            // Disable the alpha mask for this pixel
            A(x_len,j) = 1;
        }
        else if(flag==3)
        {
            C(x_len,j)=0.13;
            C2(x_len,j)=0.13;
            C3(x_len,j)=0.13;
            ray_intersection=e+t1*d;
            ray_normal = (ray_intersection-c_2).normalized();//法向量
            //光源1
            // Compute normal at the intersection point
            Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
            //光源2
            // Compute normal at the intersection point
            Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量
            if(if_hit_object(ray_intersection,ray_L)==0)
            {
                C(x_len,j) =C(x_len,j)+light_render(ray_L,-d.normalized(),ray_normal);
            }
            if(if_hit_object(ray_intersection,ray_L2)==0)
            {
                C(x_len,j) =C(x_len,j)+light_render(ray_L2,-d.normalized(),ray_normal);
            }
        
            C2(x_len,j)=0.3*C(x_len,j);
            C3(x_len,j)=0.8*C(x_len,j);
            // Disable the alpha mask for this pixel
            A(x_len,j) = 1;
        }
        else if(flag==4)
        {
            C(x_len,j)=0.13;
            C2(x_len,j)=0.13;
            C3(x_len,j)=0.13;
            ray_intersection=e+t1*d;
            //光源1
            Vector3d ray_L = (light_position-ray_intersection).normalized();//L向量
            //光源2
            Vector3d ray_L2 = (light_position2-ray_intersection).normalized();//L向量
            if(if_hit_object(ray_intersection,ray_L)==0)
            {
                C(x_len,j) =C(x_len,j)+0.6*light_render(ray_L,-d.normalized(),ray_normal);
                C2(x_len,j) = C(x_len,j);
                C3(x_len,j)= C(x_len,j);
            }
            if(if_hit_object(ray_intersection,ray_L2)==0)
            {
                C(x_len,j) =C(x_len,j)+0.6*light_render(ray_L2,-d.normalized(),ray_normal);
                C2(x_len,j) = C(x_len,j);
                C3(x_len,j)= C(x_len,j);
            }
            //增加倒影
            //反射光线，不论地面的某一点是否在阴影下都要计算反射光线
            //计算反射光线与所有物体交点，取最近的那个
            double coefficient=d.normalized().transpose()*ray_normal;
            coefficient=2*coefficient;
            Vector3d r=d.normalized()-coefficient*ray_normal;
            Vector3d ray_normal2;
            flag=-1;
            t1=caculate_intersection(ray_intersection,r,ray_normal2,flag);
            if(t1==1000000)//反射光没有和任何物体相撞
            {
                C(x_len,j) = C(x_len,j)+0;
                C2(x_len,j)=C2(x_len,j)+0;
                C3(x_len,j)=C3(x_len,j)+0;
                A(x_len,j) = 1;
            }
            else
            {
                if(flag==0)
                {
                    Vector3d ray_intersection2=ray_intersection+t1*r;
                    //光源1
                    Vector3d ray_L = (light_position-ray_intersection2).normalized();//L向量
                    //光源2
                    Vector3d ray_L2 = (light_position2-ray_intersection2).normalized();//L向量
                    if(if_hit_object(ray_intersection2,ray_L)==0)
                    {
                        C(x_len,j) =C(x_len,j)+0.6*0.3*light_render(ray_L,-r.normalized(),ray_normal2);
                        C2(x_len,j)=C2(x_len,j)+0.6*light_render(ray_L,-r.normalized(),ray_normal2);
                        C3(x_len,j)=C3(x_len,j)+0.6*0.8*light_render(ray_L,-r.normalized(),ray_normal2);
                    }
                    if(if_hit_object(ray_intersection2,ray_L2)==0)
                    {
                        C(x_len,j) =C(x_len,j)+0.6*0.3*light_render(ray_L2,-r.normalized(),ray_normal2);
                        C2(x_len,j)=C2(x_len,j)+0.6*light_render(ray_L2,-r.normalized(),ray_normal2);
                        C3(x_len,j)=C3(x_len,j)+0.6*0.8*light_render(ray_L2,-r.normalized(),ray_normal2);
                    }
                    
                }
                else if (flag==1)
                {
                    Vector3d ray_intersection2=ray_intersection+t1*r;
                    //光源1
                    Vector3d ray_L = (light_position-ray_intersection2).normalized();//L向量
                    //光源2
                    Vector3d ray_L2 = (light_position2-ray_intersection2).normalized();//L向量
                    if(if_hit_object(ray_intersection2,ray_L)==0)
                    {
                        C(x_len,j) =C(x_len,j)+0.6*0.8*light_render(ray_L,-r.normalized(),ray_normal2);
                        C2(x_len,j)=C2(x_len,j)+0.6*light_render(ray_L,-r.normalized(),ray_normal2);
                        C3(x_len,j)=C3(x_len,j)+0.6*0.3*light_render(ray_L,-r.normalized(),ray_normal2);
                    }
                    if(if_hit_object(ray_intersection2,ray_L2)==0)
                    {
                        C(x_len,j) =C(x_len,j)+0.6*0.8*light_render(ray_L2,-r.normalized(),ray_normal2);
                        C2(x_len,j)=C2(x_len,j)+0.6*light_render(ray_L2,-r.normalized(),ray_normal2);
                        C3(x_len,j)=C3(x_len,j)+0.6*0.3*light_render(ray_L2,-r.normalized(),ray_normal2);
                    }
                }
                else if (flag==2)
                {
                    Vector3d ray_intersection2=ray_intersection+t1*r;
                    //光源1
                    Vector3d ray_L = (light_position-ray_intersection2).normalized();//L向量
                    //光源2
                    Vector3d ray_L2 = (light_position2-ray_intersection2).normalized();//L向量
                    if(if_hit_object(ray_intersection2,ray_L)==0)
                    {
                        C(x_len,j) =C(x_len,j)+0.6*diffuse_render(ray_L,ray_normal2);
                        C2(x_len,j)=C2(x_len,j)+0.6*0.8*diffuse_render(ray_L,ray_normal2);
                        C3(x_len,j)=C3(x_len,j)+0.6*0.3*diffuse_render(ray_L,ray_normal2);
                    }
                    if(if_hit_object(ray_intersection2,ray_L2)==0)
                    {
                        C(x_len,j) =C(x_len,j)+0.6*diffuse_render(ray_L2,ray_normal2);
                        C2(x_len,j)=C2(x_len,j)+0.6*0.6*diffuse_render(ray_L2,ray_normal2);
                        C3(x_len,j)=C3(x_len,j)+0.6*0.3*diffuse_render(ray_L2,ray_normal2);
                    }
                }
                else if (flag==3)
                {
                    Vector3d ray_intersection2=ray_intersection+t1*r;
                    //光源1
                    Vector3d ray_L = (light_position-ray_intersection2).normalized();//L向量
                    //光源2
                    Vector3d ray_L2 = (light_position2-ray_intersection2).normalized();//L向量
                    if(if_hit_object(ray_intersection2,ray_L)==0)
                    {
                        C(x_len,j) =C(x_len,j)+0.6*light_render(ray_L,-r.normalized(),ray_normal2);
                        C2(x_len,j)=C2(x_len,j)+0.6*0.3*light_render(ray_L,-r.normalized(),ray_normal2);
                        C3(x_len,j)=C3(x_len,j)+0.6*0.8*light_render(ray_L,-r.normalized(),ray_normal2);
                    }
                    if(if_hit_object(ray_intersection2,ray_L2)==0)
                    {
                        C(x_len,j) =C(x_len,j)+0.6*light_render(ray_L2,-r.normalized(),ray_normal2);
                        C2(x_len,j)=C2(x_len,j)+0.6*0.3*light_render(ray_L2,-r.normalized(),ray_normal2);
                        C3(x_len,j)=C3(x_len,j)+0.6*0.8*light_render(ray_L2,-r.normalized(),ray_normal2);
                    }
                }
            }
            // Disable the alpha mask for this pixel
            A(x_len,j) = 1;
        }

    }
}
void part7()
{
    //read_mesh();
    cout << "Part 1.7: TBB" << std::endl;
    const std::string filename("part7.png");
    struct timeval timeStart, timeEnd; 
    double runTime=0; 
    gettimeofday(&timeStart, NULL );
    for (x_len=0;x_len<C.cols();x_len++)
    {
        parallel_for(0,int(C.rows()),draw_each_coloum);
    }
    gettimeofday( &timeEnd, NULL ); 
    runTime = (timeEnd.tv_sec - timeStart.tv_sec ) + (double)(timeEnd.tv_usec -timeStart.tv_usec)/1000000;  
    printf("Part7 with tbb runTime is %lf\n", runTime); 
    // Save to png
    write_matrix_to_png(C,C2,C3,A,filename);
}
void part8()
{
    //read_mesh();
    cout << "Part 1.8: Muti-frame" << std::endl;
    //const std::string filename("part8_1.png");
    //e(1)=1-1.3;
    light_position(2)=light_position(2)-13*0.5;
    light_position2(2)=light_position2(2)+13*0.5;
    double angle=0;//y axis
    for(int p=0;p<23;p++)
    {
        string p1="part8_";
        string p2=to_string(p);
        string p3=".png";
        string filename=p1+p2+p3;
        if(p>=0&&p<=12)
        {
            e(1)=e(1)-0.1;
            light_position(2)=light_position(2)-0.5;
            light_position2(2)=light_position2(2)+0.5;
            for (x_len=0;x_len<C.cols();x_len++)
            {
                parallel_for(0,int(C.rows()),draw_each_coloum);
            }
            write_matrix_to_png(C,C2,C3,A,filename);
        }
        else
        {
            angle=angle+10.0;
            double angle_value=angle*M_PI/180.0;
            light_position(2)=light_position(2)-0.5;
            light_position2(2)=light_position2(2)+0.5;

            Rotation<<cos(angle_value),0,sin(angle_value),
                        0,1,0,
                        -sin(angle_value),0,cos(angle_value);

            e<<0,-0.3,2;
            e=Rotation*e; 
            for (x_len=0;x_len<C.cols();x_len++)
            {
                parallel_for(0,int(C.rows()),draw_each_coloum2);
            }
            write_matrix_to_png(C,C2,C3,A,filename);
        }

    }
}
int main()
{
    read_mesh();
    part1();
    part2();
    part3();
    part4();
    part5();
    part6();
    //In my computer tbb is installed to the directory /usr/lib, and in the cmake file I add a command to link the execution file to the lbb in /usr/lib
    //if your computer cannot find tbb, just comment the part7 and part8 which used tbb.
    part7();
    //part8();
    return 0;
}
