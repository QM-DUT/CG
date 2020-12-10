#version 330

in vec2 TexCoord;

uniform sampler2D textureMap;

out vec4 outColor;

void main()
{
    ivec2 texSize = textureSize(textureMap, 0);
    float s = float(texSize.s);
    float t = float(texSize.t);
    vec2 up=vec2(0.0,1.0/t);
    vec2 down=vec2(0.0,-1.0/t);
    vec2 left=vec2(-1.0/s,0.0);
    vec2 right=vec2(1.0/s,0.0);
    vec3 Colors[]=vec3[](
    texture(textureMap, TexCoord+left+up).rgb,
    texture(textureMap, TexCoord+up).rgb,
    texture(textureMap, TexCoord+right+up).rgb,
    texture(textureMap, TexCoord+left).rgb,
    texture(textureMap, TexCoord).rgb,
    texture(textureMap, TexCoord+right).rgb,
    texture(textureMap, TexCoord+left+down).rgb,
    texture(textureMap, TexCoord+down).rgb,
    texture(textureMap, TexCoord+right+down).rgb
    );
    vec3 color = vec3(1000., 1000., 1000.);
    for(int i=0;i<9;i++)
    {
        color=min(Colors[i],color);
    }
    outColor = vec4(color, 1.0);
}