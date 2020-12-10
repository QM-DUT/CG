// This example is heavily based on the tutorial at https://open.gl

// OpenGL Helpers to reduce the clutter
#include "Helpers.h"
//将头文件变成cpp文件
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
// GLFW is necessary to handle the OpenGL context
#include <GLFW/glfw3.h>

// Linear Algebra Library
#include <Eigen/Core>
#include <Eigen/Dense>
// Timer
#include <chrono>
#include<iostream>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <fstream>
using namespace std;
using namespace glm;
using namespace Eigen;

glm::vec3 cameraPos   = glm::vec3(0.0f, 0.0f,  3.0f);
//glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
glm::vec3 cameraUp    = glm::vec3(0.0f, 1.0f,  0.0f);
glm::vec3 lightPos(1.1f, 1.2f, 2.1f);

float deltaTime = 0.0f; // 当前帧与上一帧的时间差
float lastFrame = 0.0f; // 上一帧的时间

// VertexBufferObject wrapper

VertexBufferObject VBO_cube;
VertexBufferObject VBO_N_cube;
VertexBufferObject VBO_avg_N_cube;

VertexBufferObject VBO_bunny;
VertexBufferObject VBO_N_bunny;
VertexBufferObject VBO_avg_N_bunny;


VertexBufferObject VBO_bumpy_cube;
VertexBufferObject VBO_N_bumpy_cube;
VertexBufferObject VBO_avg_N_bumpy_cube;

MatrixXf V_cube(3,3);
MatrixXf N_cube(3,3);
MatrixXf avg_N_cube(3,3);

MatrixXf V_bunny(3,3);
MatrixXf N_bunny(3,3);
MatrixXf avg_N_bunny(3,3);
MatrixXf V_bumpy_cube(3,3);
MatrixXf N_bumpy_cube(3,3);
MatrixXf avg_N_bumpy_cube(3,3);
glm::mat4 view          = glm::mat4(1.0f);
glm::mat4 projection    = glm::mat4(1.0f);
glm::mat4 model         =glm::mat4(1.0f);
glm::mat4 camera         =glm::mat4(1.0f);

typedef struct object_copy
{
    unsigned int type;
    vec3 object_color=vec3(1.0f);
    mat4 model=mat4(1.0f);
    mat4 test=mat4(1.0f);
    // mat4 translation= mat4(1.0f);
    // mat4 rotation= mat4(1.0f);
    // mat4 scale= mat4(1.0f);
    unsigned int draw_type;
}my_copy;
vector<my_copy>copys;
vector<my_copy>::iterator iter_global;
vector<my_copy>::iterator iter_global_old;
bool if_select_object=false;
bool if_mouse_click=false;
bool if_first_mouse_move=true;
double if_ray_hit_triangle(Vector3f e, Vector3f d,Vector3f Va,Vector3f Vb,Vector3f Vc)
{
    //直线与三角形的交点
    MatrixXf Mb = MatrixXf::Zero(3,3);
    MatrixXf Mr = MatrixXf::Zero(3,3);
    MatrixXf Mt = MatrixXf::Zero(3,3);
    MatrixXf Ma = MatrixXf::Zero(3,3);
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
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    // Get the position of the mouse in the window
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);

    // Get the size of the window
    int width, height;
    glfwGetWindowSize(window, &width, &height);

    // Convert screen position to world coordinates
    double x_nds = ((xpos/double(width))*2)-1;
    double y_nds = (((height-1-ypos)/double(height))*2)-1; // NOTE: y axis is flipped in glfw
    vec3 ray_nds=vec3(x_nds, y_nds, 1.0f);
    vec4 ray_clip = vec4(ray_nds.x, ray_nds.y,-1.0, 1.0);
    vec4 ray_eye = inverse(projection) * ray_clip;
    ray_eye = vec4(ray_eye.x,ray_eye.y, -1.0, 0.0);//////////////////////这里为什么是-1而不是-0.1
    vec3 ray_wor = vec3((inverse(view) *inverse(camera)* ray_eye).x,(inverse(view)*inverse(camera) * ray_eye).y,(inverse(view) *inverse(camera)* ray_eye).z);
    //vec4 ray_wor_pre=inverse(camera)*inverse(view) * ray_eye;



    ray_wor = normalize(ray_wor);//(0,0,3)+td
    // Update the position of the first vertex if the left button is pressed
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    {
        //开始判断射线与哪个个物体相交
        cout<<"press"<<endl;
        //遍历整个copys,得到最近的那个物体
        double t0=0;
        double t1=1000000;
        double t=0;;
        if_select_object=false;
        if_mouse_click=true;
        for(vector<my_copy>::iterator iter=copys.begin();iter!=copys.end();iter++)
        {
            model=mat4(1.0f);
            Vector3f e(cameraPos[0],cameraPos[1],cameraPos[2]);
            Vector3f d(ray_wor[0],ray_wor[1],ray_wor[2]);
            iter->object_color=vec3(1.0,1.0,1.0);
            //先判断是哪种类型的object，从对应的内存中读出来
            switch (iter->type)
            {
            case 1:
                for(int i=0;i<V_cube.cols();i=i+3)
                {
                    Vector3f Va=V_cube.col(i+0);
                    Vector3f Vb=V_cube.col(i+1);
                    Vector3f Vc=V_cube.col(i+2);
                    MatrixXf temp(4,4);
                    mat4 temp2=(iter->model)*(iter->test);
                    temp<<temp2[0][0],temp2[1][0],temp2[2][0],temp2[3][0],
                        temp2[0][1],temp2[1][1],temp2[2][1],temp2[3][1],
                        temp2[0][2],temp2[1][2],temp2[2][2],temp2[3][2],
                        temp2[0][3],temp2[1][3],temp2[2][3],temp2[3][3];
                    Vector4f Va4(Va(0),Va(1),Va(2),1.0);
                    Vector4f Vb4(Vb(0),Vb(1),Vb(2),1.0);
                    Vector4f Vc4(Vc(0),Vc(1),Vc(2),1.0);
                    Vector4f Va4temp,Vb4temp,Vc4temp;
                    Va4temp=temp*Va4;
                    Vb4temp=temp*Vb4;
                    Vc4temp=temp*Vc4;
                    Va<<Va4temp(0),Va4temp(1),Va4temp(2);
                    Vb<<Vb4temp(0),Vb4temp(1),Vb4temp(2);
                    Vc<<Vc4temp(0),Vc4temp(1),Vc4temp(2);

                    t=if_ray_hit_triangle(e,d,Va,Vb,Vc);
                    if(t>t0&&t<t1)
                    {
                        t1=t;
                        //iter_global_old=iter_global;
                        iter_global=iter;
                        if_select_object=true;
                    }
                }
                break;
            case 2:
                for(int i=0;i<V_bunny.cols();i=i+3)
                {
                    Vector3f Va=V_bunny.col(i+0);
                    Vector3f Vb=V_bunny.col(i+1);
                    Vector3f Vc=V_bunny.col(i+2);
                    MatrixXf temp(4,4);
                    mat4 temp2=(iter->model)*(iter->test);
                    temp<<temp2[0][0],temp2[1][0],temp2[2][0],temp2[3][0],
                        temp2[0][1],temp2[1][1],temp2[2][1],temp2[3][1],
                        temp2[0][2],temp2[1][2],temp2[2][2],temp2[3][2],
                        temp2[0][3],temp2[1][3],temp2[2][3],temp2[3][3];
                    Vector4f Va4(Va(0),Va(1),Va(2),1.0);
                    Vector4f Vb4(Vb(0),Vb(1),Vb(2),1.0);
                    Vector4f Vc4(Vc(0),Vc(1),Vc(2),1.0);
                    Vector4f Va4temp,Vb4temp,Vc4temp;
                    Va4temp=temp*Va4;
                    Vb4temp=temp*Vb4;
                    Vc4temp=temp*Vc4;
                    Va<<Va4temp(0),Va4temp(1),Va4temp(2);
                    Vb<<Vb4temp(0),Vb4temp(1),Vb4temp(2);
                    Vc<<Vc4temp(0),Vc4temp(1),Vc4temp(2);
                    t=if_ray_hit_triangle(e,d,Va,Vb,Vc);
                    if(t>t0&&t<t1)
                    {
                        t1=t;
                        //iter_global_old=iter_global;
                        iter_global=iter;
                        if_select_object=true;
                    }
                }
                break;
            case 3:
                for(int i=0;i<V_bumpy_cube.cols();i=i+3)
                {
                    Vector3f Va=V_bumpy_cube.col(i+0);
                    Vector3f Vb=V_bumpy_cube.col(i+1);
                    Vector3f Vc=V_bumpy_cube.col(i+2);
                    MatrixXf temp(4,4);
                    mat4 temp2=(iter->model)*(iter->test);
                    temp<<temp2[0][0],temp2[1][0],temp2[2][0],temp2[3][0],
                        temp2[0][1],temp2[1][1],temp2[2][1],temp2[3][1],
                        temp2[0][2],temp2[1][2],temp2[2][2],temp2[3][2],
                        temp2[0][3],temp2[1][3],temp2[2][3],temp2[3][3];
                    Vector4f Va4(Va(0),Va(1),Va(2),1.0);
                    Vector4f Vb4(Vb(0),Vb(1),Vb(2),1.0);
                    Vector4f Vc4(Vc(0),Vc(1),Vc(2),1.0);
                    Vector4f Va4temp,Vb4temp,Vc4temp;
                    Va4temp=temp*Va4;
                    Vb4temp=temp*Vb4;
                    Vc4temp=temp*Vc4;
                    Va<<Va4temp(0),Va4temp(1),Va4temp(2);
                    Vb<<Vb4temp(0),Vb4temp(1),Vb4temp(2);
                    Vc<<Vc4temp(0),Vc4temp(1),Vc4temp(2);
                    t=if_ray_hit_triangle(e,d,Va,Vb,Vc);
                    if(t>t0&&t<t1)
                    {
                        t1=t;
                        //iter_global_old=iter_global;
                        iter_global=iter;
                        if_select_object=true;
                    }
                }
                break;
            default:
                break;
            }
            

        }
        if(if_select_object==true)
        {
            //iter_global_old->object_color=vec3(1.0,1.0,1.0);
            iter_global->object_color=vec3(0.3,0.2,0.9);
        }
    }

    
}
void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
    // Get the size of the window
    int width, height;
    glfwGetWindowSize(window, &width, &height);
    // Convert screen position to world coordinates
    double x_nds = ((xpos/double(width))*2)-1;
    double y_nds = (((height-1-ypos)/double(height))*2)-1; // NOTE: y axis is flipped in glfw
    vec3 ray_nds=vec3(x_nds, y_nds, 1.0f);
    vec4 ray_clip = vec4(ray_nds.x, ray_nds.y,-1.0, 1.0);
    vec4 ray_eye = inverse(projection) * ray_clip;
    ray_eye = vec4(ray_eye.x,ray_eye.y, -1.0, 0.0);//////////////////////这里为什么是-1而不是-0.1
    vec3 ray_wor = vec3((inverse(view) * ray_eye).x,(inverse(view) * ray_eye).y,(inverse(view) * ray_eye).z);
    //ray_wor = normalize(ray_wor);//(0,0,3)+td
    // if(if_select_object==true)
    // {
    //     iter_global->translation[3][0]=ray_wor[0];
    //     iter_global->translation[3][1]=ray_wor[1];
    //     // mat4 temp=iter_global->model;
    //     // iter_global->model=(iter_global->translation)*temp;
    // }
    

}
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    // Update the position of the first vertex if the keys 1,2, or 3 are pressed
    float cameraSpeed = 300 * deltaTime;
    if(action == GLFW_PRESS)
    {
        switch (key)
        {
            case  GLFW_KEY_1:
            {
                my_copy copy_temp;
                copy_temp.type=1;//cube
                copy_temp.draw_type=2;
                copy_temp.model[0][0]=0.5;
                copy_temp.model[1][1]=0.5;
                copy_temp.model[2][2]=0.5;
                //copy_temp.model[3][1]=-0.5;
                //copy_temp.model = glm::translate(copy_temp.model, glm::vec3(0.0, -0.5, 0.0));
                //copy_temp.test= glm::translate(copy_temp.test, glm::vec3(0.0, -0.5, 0.0));
                copys.push_back(copy_temp);
                break;
            } 
            case  GLFW_KEY_2:
            {
                my_copy copy_temp;
                copy_temp.type=2;//bunny
                copy_temp.draw_type=2;
                copy_temp.model[0][0]=3;
                copy_temp.model[1][1]=3;
                copy_temp.model[2][2]=3;
                copy_temp.test= glm::translate(copy_temp.test, glm::vec3(0.0, -0.4, 0.0));
                copy_temp.test[3][1]=-0.1;
                // mat4 temp=mat4(1.0);
                // temp[3][1]=-0.1;
                // copy_temp.model=copy_temp.model*temp;

                copys.push_back(copy_temp);
                break;
            } 
            case  GLFW_KEY_3:
            {
                my_copy copy_temp;
                copy_temp.type=3;//bumty_cube
                copy_temp.draw_type=2;
                copy_temp.model[0][0]=0.1;
                copy_temp.model[1][1]=0.1;
                copy_temp.model[2][2]=0.1;
                copys.push_back(copy_temp);
                break;
            } 
            case  GLFW_KEY_W:
                cameraPos += cameraSpeed * glm::normalize(-cameraPos);//摄像机前移
                break;  
            case  GLFW_KEY_S:
                cameraPos -= cameraSpeed * glm::normalize(-cameraPos);//摄像机后移
                break;
            case  GLFW_KEY_Z://摄像机左移
                //cameraPos -= glm::normalize(glm::cross(cameraPos, cameraUp)) * cameraSpeed;
                cameraPos+=vec3(-1,0,0);
                break;  
            case  GLFW_KEY_C://摄像机右移
                //cameraPos += glm::normalize(glm::cross(cameraPos, cameraUp)) * cameraSpeed;
                cameraPos+=vec3(1,0,0);
                break;
            case  GLFW_KEY_X://摄像机上移
                //cameraPos -= glm::normalize(glm::cross(cameraPos, cameraUp)) * cameraSpeed;
                cameraPos+=vec3(0,1,0);
                break;  
            case  GLFW_KEY_V://摄像机下移
                //cameraPos += glm::normalize(glm::cross(cameraPos, cameraUp)) * cameraSpeed;
                cameraPos+=vec3(0,-1,0);
                break;
            case  GLFW_KEY_A:
            {
                //cameraPos -= glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;//摄像机左移
                float left_angle=20.0;
                mat4 temp=mat4(1.0);
                vec3 left=glm::cross(-cameraPos, cameraUp);
                vec3 axi=glm::normalize(glm::cross(left, cameraPos));
                temp = glm::rotate(temp, glm::radians(left_angle),axi);
                vec4 cameraPos_v4=temp*vec4(cameraPos,1.0);
                cameraPos=vec3(cameraPos_v4.x,cameraPos_v4.y,cameraPos_v4.z);
                break;
            }
            case  GLFW_KEY_D:
            {
                //cameraPos -= glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;//摄像机左移
                float left_angle=-20.0;
                mat4 temp=mat4(1.0);
                vec3 left=glm::cross(-cameraPos, cameraUp);
                vec3 axi=glm::normalize(glm::cross(left, cameraPos));
                temp = glm::rotate(temp, glm::radians(left_angle),axi);
                vec4 cameraPos_v4=temp*vec4(cameraPos,1.0);
                cameraPos=vec3(cameraPos_v4.x,cameraPos_v4.y,cameraPos_v4.z);
                break;
            }
            case  GLFW_KEY_Q://上旋
            {
                float up_angle=-40.0;
                mat4 temp=mat4(1.0);
                vec3 axi=glm::normalize(glm::cross(-cameraPos, cameraUp)); 
                if(axi[0]<0)
                {
                    axi=-axi;
                }
                temp = glm::rotate(temp, glm::radians(up_angle),axi);
                vec4 cameraPos_v4=temp*vec4(cameraPos,1.0);
                cameraPos=vec3(cameraPos_v4.x,cameraPos_v4.y,cameraPos_v4.z);
                if(cameraPos[2]<=0)
                {
                    cameraUp=vec3(0,-1,0);
                }
                else
                {
                    cameraUp=vec3(0,1,0);
                }
                break;
            }
            case  GLFW_KEY_E:
            {
                float up_angle=40.0;
                mat4 temp=mat4(1.0);
                vec3 axi=glm::normalize(glm::cross(-cameraPos, cameraUp));
                if(axi[0]<0)
                {
                    axi=-axi;
                }
                temp = glm::rotate(temp, glm::radians(up_angle),axi);
                vec4 cameraPos_v4=temp*vec4(cameraPos,1.0);
                cameraPos=vec3(cameraPos_v4.x,cameraPos_v4.y,cameraPos_v4.z);
                if(cameraPos[2]<=0)
                {
                    cameraUp=vec3(0,-1,0);
                }
                else
                {
                    cameraUp=vec3(0,1,0);
                }
                break;
            }


            case  GLFW_KEY_U:
            {
                if(copys.size()!=0&&if_select_object==true)
                {
                    iter_global->draw_type=1;
                }
                break;
            } 
            case  GLFW_KEY_I:
            {
                if(copys.size()!=0&&if_select_object==true)
                {
                    iter_global->draw_type=2;
                }
                break;
            }       
            case  GLFW_KEY_O:
            {
                if(copys.size()!=0&&if_select_object==true)
                {
                    iter_global->draw_type=3;
                }
                break;
            }    
            case  GLFW_KEY_P://deletion
            {
                if(copys.size()!=0&&if_select_object==true)
                {
                    //swap
                    // if(copys.size()==1)
                    // {
                    //     copys.pop_back();
                    // }
                    // else if(copys.size()>1)
                    // {
                    //     std::swap(*iter_global,copys[-1]);
                    //     copys.pop_back();
                    // }
                    // cout<<44444444444444<<endl;
                    copys.erase(iter_global);
                    if_select_object=false;
                    //cout<<copys.size()<<endl;
                }
                break;
            }
            case  GLFW_KEY_KP_4://left
            {
                if(if_select_object==true)
                {
                    mat4 temp=mat4(1.0f);
                    temp[3][0]=-0.5;
                    iter_global->model=temp*(iter_global->model);
                }
                break;
            }
            case  GLFW_KEY_KP_6://right
            {
                if(if_select_object==true)
                {
                    mat4 temp=mat4(1.0f);
                    temp[3][0]=0.5;
                    iter_global->model=temp*(iter_global->model);
                }
                break;
            }
            case  GLFW_KEY_KP_8://up
            {
                if(if_select_object==true)
                {
                    mat4 temp=mat4(1.0f);
                    temp[3][1]=0.5;
                    iter_global->model=temp*(iter_global->model);
                }
                break;
            }
            case  GLFW_KEY_KP_2://down
            {
                if(if_select_object==true)
                {
                    mat4 temp=mat4(1.0f);
                    temp[3][1]=-0.5;
                    iter_global->model=temp*(iter_global->model);
                }
                break;
            }
            case  GLFW_KEY_KP_1://forward
            {
                if(if_select_object==true)
                {
                    mat4 temp=mat4(1.0f);
                    temp[3][2]=-0.5;
                    iter_global->model=temp*(iter_global->model);
                }
                break;
            }
            case  GLFW_KEY_KP_3://back
            {
                if(if_select_object==true)
                {
                    mat4 temp=mat4(1.0f);
                    temp[3][2]=0.5;
                    iter_global->model=temp*(iter_global->model);
                }
                break;
            }
            case  GLFW_KEY_KP_7://left rotate
            {
                if(if_select_object==true)
                {
                   // mat4 temp=mat4(1.0f);
                    //temp[3][1]=-0.5;
                    //iter_global->model=temp*(iter_global->model);
                    float left_angle=20.0;
                    iter_global->model = glm::rotate(iter_global->model, glm::radians(left_angle), glm::vec3(0.0f, 0.0f, 1.0f));
                }
                break;
            }
            case  GLFW_KEY_KP_9://right rotate
            {
                if(if_select_object==true)
                {
                    //mat4 temp=mat4(1.0f);
                    // temp[3][1]=-0.5;
                    // iter_global->model=temp*(iter_global->model);
                    float right_angle=-20.0;
                    iter_global->model = glm::rotate(iter_global->model, glm::radians(right_angle), glm::vec3(0.0f, 0.0f, 1.0f));
                }
                break;
            }
            case  GLFW_KEY_KP_ADD://lagre
            {
                if(if_select_object==true)
                {
                    mat4 temp=mat4(1.0f);
                    temp[0][0]=2;
                    temp[1][1]=2;
                    temp[2][2]=2;
                    iter_global->model=(iter_global->model)*temp;
                    //iter_global->model = glm::scale(iter_global->model, glm::vec3(2.0, 2.0, 2.0)); 

                }
                break;
            }
            case  GLFW_KEY_KP_SUBTRACT://small
            {
                if(if_select_object==true)
                {
                    // mat4 temp=mat4(1.0f);
                    // temp[0][0]=0.5;
                    // temp[1][1]=0.5;
                    // temp[2][2]=0.5;
                    //iter_global->model=temp*(iter_global->model);
                    iter_global->model = glm::scale(iter_global->model, glm::vec3(0.5, 0.5, 0.5)); 
                }
                break;
            }
            default:
                break;
        }
    }
}
void read_mesh(string&file_path,MatrixXf&V,MatrixXi&F)
{
    int nverts = 0;//顶点个数
    int nfaces = 0;//面个数
    int nedges = 0;//边个数
    int line_count = 1;//读入行数
    
    ifstream file(file_path.c_str());
	string line;
	if(!file) // 有该文件
	{
        cout <<"no such file" << endl;
        return;
	}
	while (getline (file, line)) // line中不包括每行的换行符
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
    V.resize(nverts,3);
    F.resize(nfaces,3);

    while (getline (file, line))
    {
        if(line_count<nverts+3)//505
        {
            int pos=0;
            pos=line.find(" ",0);
            V(line_count-3,0)=atof(line.substr(0,pos).c_str());
            line=line.substr(pos+1,-1);
            pos=line.find(" ",0);
            V(line_count-3,1)=atof(line.substr(0,pos).c_str());
            V(line_count-3,2)=atof(line.substr(pos+1,-1).c_str());
        }
        else if(line_count<nverts+3+nfaces)//1505
        {
            int pos=0;
            pos=line.find(" ",0);
            line=line.substr(pos+1,-1);


            pos=line.find(" ",0);
            F(line_count-3-nverts,0)=atoi(line.substr(0,pos).c_str());
            line=line.substr(pos+1,-1);
            pos=line.find(" ",0);
            F(line_count-3-nverts,1)=atoi(line.substr(0,pos).c_str());
            F(line_count-3-nverts,2)=atoi(line.substr(pos+1,-1).c_str());
        }
        line_count++;
    }
    V.transposeInPlace();
    F.transposeInPlace();
}
void calculate_face_normal(MatrixXf &V_temp,MatrixXi &E_temp,MatrixXf &V_cube,MatrixXf &N_cube)
{
    V_cube.resize(3,E_temp.cols()*3);
    N_cube.resize(3,E_temp.cols()*3);
    for(int i=0;i<E_temp.cols();i++)
    {
        V_cube.col(3*i+0)=V_temp.col(E_temp(0,i));
        V_cube.col(3*i+1)=V_temp.col(E_temp(1,i));
        V_cube.col(3*i+2)=V_temp.col(E_temp(2,i));
        Vector3f Normal;
        Vector3f ab=V_cube.col(3*i+1)-V_cube.col(3*i+0);
        Vector3f ac=V_cube.col(3*i+2)-V_cube.col(3*i+0);
        Normal = ab.cross(ac).normalized();//法向量

        N_cube.col(3*i+0)=Normal;
        N_cube.col(3*i+1)=Normal;
        N_cube.col(3*i+2)=Normal;
    }
}
void calculate_vertex_normal(MatrixXf &V_temp,MatrixXi &E_temp,MatrixXf &N_cube,MatrixXf &avg_N_cube)
{
    avg_N_cube.resize(3,E_temp.cols()*3);
    Vector3f normal(0.0,0.0,0.0);
    //对于每个顶点，通过V_temp知道有多少个点，通过E_temp算出所有与他相关的面
    for(int i=0;i<V_temp.cols();i++)
    {
        for(int j=0;j<E_temp.cols();j++)
        {
            if(E_temp(0,j)==i||E_temp(1,j)==i||E_temp(2,j)==i)//找到了这个面j
            {
                normal+=N_cube.col(3*j);
            }
        }
        //标准化
        normal.normalize();
        for(int k=0;k<E_temp.cols();k++)
        {
            if(E_temp(0,k)==i)//找到了这个面k
            {
                avg_N_cube.col(3*k+0)=normal;
            }
            else if(E_temp(1,k)==i)//找到了这个面k
            {
                avg_N_cube.col(3*k+1)=normal;
            }
            else if(E_temp(2,k)==i)//找到了这个面k
            {
                avg_N_cube.col(3*k+2)=normal;
            }
        }
    }
}
int main(void)
{


    GLFWwindow* window;
    // Initialize the library
    if (!glfwInit())
        return -1;

    // Activate supersampling
    glfwWindowHint(GLFW_SAMPLES, 8);

    // Ensure that we get at least a 3.2 context
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);

    // On apple we have to load a core profile with forward compatibility
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // Create a windowed mode window and its OpenGL context
    window = glfwCreateWindow(640, 480, "Hello World", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    // Make the window's context current
    glfwMakeContextCurrent(window);

    #ifndef __APPLE__
      glewExperimental = true;
      GLenum err = glewInit();
      if(GLEW_OK != err)
      {
        /* Problem: glewInit failed, something is seriously wrong. */
       fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
      }
      glGetError(); // pull and savely ignonre unhandled errors like GL_INVALID_ENUM
      fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));
    #endif

    int major, minor, rev;
    major = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MAJOR);
    minor = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MINOR);
    rev = glfwGetWindowAttrib(window, GLFW_CONTEXT_REVISION);
    printf("OpenGL version recieved: %d.%d.%d\n", major, minor, rev);
    printf("Supported OpenGL is %s\n", (const char*)glGetString(GL_VERSION));
    printf("Supported GLSL is %s\n", (const char*)glGetString(GL_SHADING_LANGUAGE_VERSION));

    




    MatrixXf V_temp(3,3);
    MatrixXi E_temp(3,3);
    VertexArrayObject VAO_cube;
    VAO_cube.init();
    VertexArrayObject VAO_avg_cube;
    VAO_avg_cube.init();
    string cube_file_path="../data/cube.off";
    read_mesh(cube_file_path,V_temp,E_temp);
    calculate_face_normal(V_temp,E_temp,V_cube,N_cube);
    cout<<V_cube<<endl;
    cout<<N_cube<<endl;
    calculate_vertex_normal(V_temp,E_temp,N_cube,avg_N_cube);
    cout<<avg_N_cube<<endl;
    VBO_cube.init();
    VBO_cube.update(V_cube);
    VBO_N_cube.init();
    VBO_N_cube.update(N_cube);
    VBO_avg_N_cube.init();
    VBO_avg_N_cube.update(avg_N_cube);


    VertexArrayObject VAO_bunny;
    VAO_bunny.init();
    VertexArrayObject VAO_avg_bunny;
    VAO_avg_bunny.init();
    string bunny_file_path="../data/bunny.off";
    read_mesh(bunny_file_path,V_temp,E_temp);
    calculate_face_normal(V_temp,E_temp,V_bunny,N_bunny);
    calculate_vertex_normal(V_temp,E_temp,N_bunny,avg_N_bunny);
    VBO_bunny.init();
    VBO_bunny.update(V_bunny);
    VBO_N_bunny.init();
    VBO_N_bunny.update(N_bunny);
    VBO_avg_N_bunny.init();
    VBO_avg_N_bunny.update(avg_N_bunny);
    
    VertexArrayObject VAO_bumpy_cube;
    VAO_bumpy_cube.init();
    VertexArrayObject VAO_avg_bumpy_cube;
    VAO_avg_bumpy_cube.init();
    string bumpy_cube_file_path="../data/bumpy_cube.off";
    read_mesh(bumpy_cube_file_path,V_temp,E_temp);
    calculate_face_normal(V_temp,E_temp,V_bumpy_cube,N_bumpy_cube);
    calculate_vertex_normal(V_temp,E_temp,N_bumpy_cube,avg_N_bumpy_cube);
    VBO_bumpy_cube.init();
    VBO_bumpy_cube.update(V_bumpy_cube);
    VBO_N_bumpy_cube.init();
    VBO_N_bumpy_cube.update(N_bumpy_cube);
    VBO_avg_N_bumpy_cube.init();
    VBO_avg_N_bumpy_cube.update(avg_N_bumpy_cube);


    Program program;
    const GLchar* vertex_shader =
            "#version 330 core\n"//这里必须敲\n
                    "layout (location = 0) in vec3 position;"
                    "layout (location = 1) in vec3 aNormal;"
                    "out vec3 Normal;"
                    "out vec3 WorldPos;"
                    "uniform mat4 projection;"
                    "uniform mat4 view;"
                    "uniform mat4 model;"
                    "uniform mat4 camera;"
                    "uniform mat4 test;"
                    "void main()"
                    "{"
                    "    WorldPos = vec3(model *test* vec4(position, 1.0));"
                    "    gl_Position = projection * camera* view *vec4(WorldPos, 1.0);"
                    //"    Normal = aNormal;"
                    "    Normal = mat3(transpose(inverse(model))) * aNormal;"
                    "}";

    const GLchar* fragment_shader =
            "#version 330 core\n"
                    "out vec4 outColor;"
                    "in vec3 Normal;"
                    "in vec3 WorldPos;" 
                    "uniform vec3 lightPos;"
                    "uniform vec3 lightColor;"
                    "uniform vec3 objectColor;"
                    "uniform vec3 viewPos;"
                    "void main()"
                    "{"
                    "    float ambientCoefficient = 0.1;"
                    "    vec3 ambient = ambientCoefficient * lightColor;"
                    "    vec3 norm = normalize(Normal);"
                    "    vec3 ray_L = normalize(lightPos - WorldPos);"
                    "    vec3 ray_V = normalize(viewPos - WorldPos);"
                    "    if(dot(norm, ray_V)<0)"
                    "    {"
                    "        norm=-norm;"
                    "    }"
                    "    vec3 ray_H =normalize(ray_L+ray_V);"
                    "    float temp=dot(ray_L,norm);"
                    "    float diffuse = max(temp, 0.);"
                    "    float temp2=dot(norm,ray_H);"
                    "    float Phong = pow(max(temp2, 0.),1000);"
                    "    vec3 result = (ambient + diffuse+ Phong) * objectColor;"
                    "    outColor = vec4(result, 1.0);"
                    "}";


    program.init(vertex_shader,fragment_shader,"outColor");
    program.bind();


    VAO_cube.bind();
    program.bindVertexAttribArray("position",VBO_cube);
    program.bindVertexAttribArray("aNormal",VBO_N_cube);

    VAO_avg_cube.bind();
    program.bindVertexAttribArray("position",VBO_cube);
    program.bindVertexAttribArray("aNormal",VBO_avg_N_cube);

    VAO_bunny.bind();
    program.bindVertexAttribArray("position",VBO_bunny);
    program.bindVertexAttribArray("aNormal",VBO_N_bunny);

    VAO_avg_bunny.bind();
    program.bindVertexAttribArray("position",VBO_bunny);
    program.bindVertexAttribArray("aNormal",VBO_avg_N_bunny);

    VAO_bumpy_cube.bind();
    program.bindVertexAttribArray("position",VBO_bumpy_cube);
    program.bindVertexAttribArray("aNormal",VBO_N_bumpy_cube);

    VAO_avg_bumpy_cube.bind();
    program.bindVertexAttribArray("position",VBO_bumpy_cube);
    program.bindVertexAttribArray("aNormal",VBO_avg_N_bumpy_cube);

    // Save the current time --- it will be used to dynamically change the triangle color
    auto t_start = std::chrono::high_resolution_clock::now();

    // Register the keyboard callback
    glfwSetKeyCallback(window, key_callback);

    // Register the mouse callback
    glfwSetMouseButtonCallback(window, mouse_button_callback);

    // Update viewport
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    
    //Register the Cursor callback
    glfwSetCursorPosCallback(window, cursor_position_callback);
    
    // Loop until the user closes the window
    float aspect_ratio_pre=1.0/640;
    float aspect_ratio2_pre=1.0/480;
    while (!glfwWindowShouldClose(window))
    {

        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;
        // Bind your program
        program.bind();
        glEnable(GL_DEPTH_TEST);
        // Set the uniform value depending on the time difference
        auto t_now = std::chrono::high_resolution_clock::now();
        float time = std::chrono::duration_cast<std::chrono::duration<float>>(t_now - t_start).count();
        glUniform3f(program.uniform("triangleColor"), (float)(sin(time * 4.0f) + 1.0f) / 2.0f, 0.0f, 0.0f);

        // Clear the framebuffer
        glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        



        

        // Get size of the window
        int width, height;
        glfwGetWindowSize(window, &width, &height);
        float aspect_ratio = 1.0/float(width); // corresponds to the necessary width scaling
        float aspect_ratio2=1.0/float(height);
        camera[0][0]=camera[0][0]*(1.0/aspect_ratio_pre)*aspect_ratio;
        camera[1][1]=camera[1][1]*(1.0/aspect_ratio2_pre)*aspect_ratio2;
        aspect_ratio_pre = 1.0/float(width); // corresponds to the necessary width scaling
        aspect_ratio2_pre=1.0/float(height);
        program.setMat4("camera", camera);

        //view = glm::lookAt(cameraPos, vec3(0.0,0.0,0.0), cameraUp);
        view = glm::lookAt(cameraPos, vec3(0.0,0.0,0.0), cameraUp);
        projection = glm::perspective(glm::radians(45.0f), (float)640 / (float)480, 0.1f, 100.0f);
        unsigned int viewLoc  = program.uniform("view");
        unsigned int projectionLoc  = program.uniform("projection");
        glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view));
        glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, glm::value_ptr(projection));
        program.setMat4("model", model);


        //draw cube
        for(vector<my_copy>::iterator iter=copys.begin();iter!=copys.end();iter++)
        {
            model=mat4(1.0f);
            program.setMat4("model", iter->model);
            program.setMat4("test", iter->test);
            //这些东西是公用的，不随copy变化而变化
            program.setVec3("lightPos", lightPos);
            program.setVec3("objectColor", iter->object_color);
            program.setVec3("lightColor", 1.0f, 1.0f, 1.0f);
            program.setVec3("viewPos", cameraPos); 
            //根据物体种类
            switch (iter->type)
            {
            case 1:
            {
                switch (iter->draw_type)
                {
                case 1:
                    VAO_cube.bind();
                    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                    break;
                case 2:
                    VAO_cube.bind();
                    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                    break;
                case 3:
                    VAO_avg_cube.bind();
                    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                    break;
                default:
                    break;
                }
                //VAO_cube.bind();
                glDrawArrays(GL_TRIANGLES, 0, 36);
                break;
            }
            case 2:
            {
                switch (iter->draw_type)
                {
                case 1:
                    //VBO_N_bunny.update(N_bunny);
                    VAO_bunny.bind();
                    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                    break;
                case 2:
                    //VBO_N_bunny.update(N_bunny);
                    VAO_bunny.bind();
                    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                    break;
                case 3:
                    VAO_avg_bunny.bind();
                    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                    //VBO_N_bunny.update(avg_N_bunny);
                    break;
                default:
                    break;
                }
                //VAO_bunny.bind();
                glDrawArrays(GL_TRIANGLES, 0, 3000);
                break;
            }
            case 3:
            {
                switch (iter->draw_type)
                {
                case 1:
                    //VBO_N_bumpy_cube.update(N_bumpy_cube);
                    VAO_bumpy_cube.bind();
                    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
                    break;
                case 2:
                    //VBO_N_bumpy_cube.update(N_bumpy_cube);
                    VAO_bumpy_cube.bind();
                    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                    break;
                case 3:
                    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
                    VAO_avg_bumpy_cube.bind();
                    //VBO_N_bumpy_cube.update(avg_N_bumpy_cube);
                    break;
                default:
                    break;
                }
                //VAO_bumpy_cube.bind();
                glDrawArrays(GL_TRIANGLES, 0, 3000);
                break;
            }  
            default:
                break;
            }
            
            
        }

        // Swap front and back buffers
        glfwSwapBuffers(window);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        // Poll for and process events
        glfwPollEvents();
    }

    // Deallocate opengl memory
    program.free();

    // Deallocate glfw internals
    glfwTerminate();
    return 0;
}
