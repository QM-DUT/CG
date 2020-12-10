// This example is heavily based on the tutorial at https://open.gl
// OpenGL Helpers to reduce the clutter
#include "Helpers.h"
//将头文件变成cpp文件
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION	// include之前必须定义
#include "stb_image_write.h"
// GLFW is necessary to handle the OpenGL context
#include <GLFW/glfw3.h>

// Linear Algebra Library
#include <Eigen/Core>
#include <iostream>
#include <opencv2/opencv.hpp>
using namespace cv;
using namespace std;
using namespace glm;
using namespace Eigen;

int filter_num=8;
int active_index=0;
// VertexBufferObject wrapper
VertexBufferObject VBO_vertex;
VertexBufferObject VBO_texture;
VertexArrayObject*VAO=new VertexArrayObject[filter_num];
ElementBufferObject EBO;
// Contains the vertex positions
Eigen::MatrixXf V_vertex(2,3);
Eigen::MatrixXf V_texture(2,3);
Eigen::MatrixXi V_index(2,3);
VideoCapture cap;
Mat frame;
long currentFrame = 0;
enum mode_set{video_mode,image_mode} my_mode;
bool if_can_read=1;
glm::mat4 view          = glm::mat4(1.0f);
glm::mat4 model         =glm::mat4(1.0f);
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{

}
void set_texture()
{
    unsigned int texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    // set the texture wrapping/filtering options (on the currently bound texture object)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);	
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    // load and generate the texture
    int w = frame.cols;
	int h = frame.rows;
	if (frame.channels() == 3)
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_BGR_EXT, GL_UNSIGNED_BYTE, frame.data);
	else if (frame.channels() == 4)
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_BGRA_EXT, GL_UNSIGNED_BYTE, frame.data);
	else if (frame.channels() == 1)
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, frame.data);
    glGenerateMipmap(GL_TEXTURE_2D);
}
void grab(GLFWwindow* &window)
{
    // Get the size of the window
    int width, height;
    glfwGetWindowSize(window, &width, &height);
    cout<<width<<endl;
    cout<<height<<endl;
    stbi_flip_vertically_on_write(true);
    GLubyte* pPixelData;
    GLint    PixelDataLength;
    GLint    i, j;
    PixelDataLength = i * height;
    i = width * 3;   
    while( i%4 != 0 )      
        ++i;       
     
    pPixelData = (GLubyte*)malloc(PixelDataLength);
    if( pPixelData == 0 )
        exit(0);
    glReadPixels(0, 0, width, height,
        GL_RGB, GL_UNSIGNED_BYTE, pPixelData);
    //stbi_write_png("../processed2.png", width, height, 3, pPixelData, 0);
    //stbi_image_free(pPixelData);
    free(pPixelData);

}
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if(action == GLFW_PRESS)
    {
        switch (key)
        {
            case  GLFW_KEY_1:
                active_index=0;
                break;
            case  GLFW_KEY_2:
                active_index=1;
                break;
            case  GLFW_KEY_3:
                active_index=2;
                break;
            case  GLFW_KEY_4:
                active_index=3;
                break;
            case  GLFW_KEY_5:
                active_index=4;
                break;
            case  GLFW_KEY_6:
                active_index=5;
                break; 
            case  GLFW_KEY_7:
                active_index=6;
                break; 
            case  GLFW_KEY_8:
                active_index=7;
                break; 
            case  GLFW_KEY_SPACE:
            {
                if(if_can_read)
                    if_can_read=0;
                else
                    if_can_read=1;
                break; 
            }
            case GLFW_KEY_LEFT:
            {
                currentFrame=currentFrame-30;
                if(if_can_read==0)
                {
                    cap.set(CV_CAP_PROP_POS_FRAMES, currentFrame);
                    if(cap.read(frame))
                    {
                        set_texture();
                        currentFrame++;
                    }
                }
                break;
            }
            case GLFW_KEY_RIGHT:
            {
                currentFrame=currentFrame+30;
                if(if_can_read==0)
                {
                    cap.set(CV_CAP_PROP_POS_FRAMES, currentFrame);
                    if(cap.read(frame))
                    {
                        set_texture();
                        currentFrame++;
                    }
                }
                break;
            }
            case  GLFW_KEY_KP_4://left
            {
                mat4 temp=mat4(1.0f);
                temp[3][0]=+0.5;
                model=temp*model;
                break;
            }
            case  GLFW_KEY_KP_6://right
            {
                mat4 temp=mat4(1.0f);
                temp[3][0]=-0.5;
                model=temp*model;
                break;
            }
            case  GLFW_KEY_KP_8://up
            {
                mat4 temp=mat4(1.0f);
                temp[3][1]=-0.5;
                model=temp*model;
                break;
            }
            case  GLFW_KEY_KP_2://down
            {
                mat4 temp=mat4(1.0f);
                temp[3][1]=+0.5;
                model=temp*model;
                break;
            }
            case  GLFW_KEY_KP_ADD://lagre
            {
                mat4 temp=mat4(1.0f);
                temp[0][0]=2;
                temp[1][1]=2;
                temp[2][2]=2;
                model=temp*model;
                break;
            }
            case  GLFW_KEY_KP_SUBTRACT://small
            {
                mat4 temp=mat4(1.0f);
                temp[0][0]=0.5;
                temp[1][1]=0.5;
                temp[2][2]=0.5;
                model=temp*model;
                break;
            }
            case  GLFW_KEY_P://截图
            {
                grab(window);
                

            }
            default:
                break;
        }
    }
}

void read_one_frame()
{
    
    if(if_can_read)
    {
        cap.set(CV_CAP_PROP_POS_FRAMES, currentFrame);
        if(cap.read(frame))
        {
            set_texture();
            currentFrame++;
        }
    }
    
}

int main(int argc, char** argv )
{
    if ( argc != 3)
    {
        printf("usage: -video <video_add> or -image <image_add>>\n");
        return -1;
    }
    if(strcmp(argv[1],"-video")==0)
    {
        my_mode=video_mode;
        
    }
    else if(strcmp(argv[1],"-image")==0)
    {
        my_mode=image_mode;
    }
    else
    {
        printf("please input correct argument\n");
        return -1;
    }
    
    
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
    window = glfwCreateWindow(640, 480, "Image Processing in OpenGL", NULL, NULL);
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

    unsigned int texture;
    if(my_mode==video_mode)
    {
        cap.open(argv[2]);
        if(!cap.isOpened())
        {
            printf("can not open video\n");
            return -1;
        }
    }
    else
    {
        glGenTextures(1, &texture);
        glBindTexture(GL_TEXTURE_2D, texture);
        // set the texture wrapping/filtering options (on the currently bound texture object)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);	
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        // load and generate the texture
        int width, height, nrChannels;
        unsigned char *data = stbi_load(argv[2], &width, &height, &nrChannels, 0);
        if (data)
        {
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
            glGenerateMipmap(GL_TEXTURE_2D);
        }
        else
        {
            std::cout << "Failed to load texture" << std::endl;
        }
        stbi_image_free(data);
    }
    for(int i=0;i<filter_num;i++)
    {
        VAO[i].init();
    }
    
    
    VBO_texture.init();
    VBO_vertex.init();
    EBO.init();
    V_vertex.resize(2,4);
    V_vertex << -1.0,1.0,1.0,-1.0,
                -1.0,-1.0,1.0,1.0;
    V_texture.resize(2,4);
    V_texture << 0.0,1.0,1.0,0.0,
                 0.0,0.0,1.0,1.0;
    V_index.resize(3,2);
    V_index<<0,0,
            1,2,
            2,3;
    VBO_vertex.update(V_vertex);
    VBO_texture.update(V_texture);
    EBO.update(V_index);
    // Initialize the OpenGL Program
    // A program controls the OpenGL pipeline and it must contains
    // at least a vertex shader and a fragment shader to be valid

    Program*program=new Program[filter_num];
    string*titles=new string[filter_num];
    // Compile the two shaders and upload the binary to the GPU
    // Note that we have to explicitly specify that the output "slot" called outColor
    // is the one that we want in the fragment buffer (and thus on screen)
    program[0].init_by_path("../src/vertex.shader","../src/fragment.shader","outColor");
    program[1].init_by_path("../src/vertex.shader","../src/fragment_gaussian.shader","outColor");
    program[2].init_by_path("../src/vertex.shader","../src/fragment_gaussian5.shader","outColor");
    program[3].init_by_path("../src/vertex.shader","../src/fragment_gaussian7.shader","outColor");
    program[4].init_by_path("../src/vertex.shader","../src/fragment_embossing.shader","outColor");
    program[5].init_by_path("../src/vertex.shader","../src/fragment_sobel.shader","outColor");
    program[6].init_by_path("../src/vertex.shader","../src/fragment_dilation.shader","outColor");
    program[7].init_by_path("../src/vertex.shader","../src/fragment_erosion.shader","outColor");
    titles[0]="origin";
    titles[1]="gaussian";
    titles[2]="gaussian5*5";
    titles[3]="gaussian7*7";
    titles[4]="embossing";
    titles[5]="sobel";
    titles[6]="dilation";
    titles[7]="erosion";
    for(int i=0;i<filter_num+1;i++)
    {
        EBO.bind();
        VAO[i].bind();
        program[i].bindVertexAttribArray("position",VBO_vertex);
        program[i].bindVertexAttribArray("aTexCoord",VBO_texture);
    };
    // Register the keyboard callback
    glfwSetKeyCallback(window, key_callback);

    // Register the mouse callback
    glfwSetMouseButtonCallback(window, mouse_button_callback);

    // Update viewport
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    
    // Loop until the user closes the window
    while (!glfwWindowShouldClose(window))
    {
        // Bind your program
        program[active_index].bind();
        program[active_index].setMat4("model", model);
        program[active_index].setMat4("view", view);
        glfwSetWindowTitle(window,titles[active_index].c_str());
        // Clear the framebuffer
        glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);


        if(my_mode==video_mode)
        {
            read_one_frame();
        }
        else
        {
            // bind Texture
            glBindTexture(GL_TEXTURE_2D, texture);
        }
        VAO[active_index].bind();

        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);



        // Swap front and back buffers
        glfwSwapBuffers(window);

        // Poll for and process events
        glfwPollEvents();
        waitKey(30);
    };

    // Deallocate opengl memory
    
    VBO_vertex.free();
    VBO_texture.free();
    for(int i=0;i<filter_num;i++)
    {
        program[i].free();
        VAO[i].free();
    }
    
    EBO.free();
    cap.release();
    // Deallocate glfw internals
    glfwTerminate();
    return 0;
}