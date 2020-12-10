#version 330 core
in vec2 position;
in vec2 aTexCoord;
out vec2 TexCoord;
uniform mat4 view;
uniform mat4 model;
void main()
{
    gl_Position = view*model*vec4(position.x,-position.y,0.0,1.0);
    TexCoord = vec2(aTexCoord.x, aTexCoord.y);
}