#version 330

in vec2 TexCoord;

uniform sampler2D textureMap;

out vec4 outColor;

void main()
{

    float core[]=float[](2.,  0.,  0.,
                         0.,  -1.,  0.,
                          0.,  0.,  -1.);
    float factor = 1.0;
    float gray=0.0;
    ivec2 texSize = textureSize(textureMap, 0);
    float s = float(texSize.s);
    float t = float(texSize.t);
    
    vec2 up=vec2(0.0,1.0/t);
    vec2 down=vec2(0.0,-1.0/t);
    vec2 left=vec2(-1.0/s,0.0);
    vec2 right=vec2(1.0/s,0.0);   
    vec3 color = vec3(0., 0., 0.);

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
    for(int i=0;i<9;i++)
    {
        color +=Colors[i] * core[i];
    }
    color /= factor;
    gray=max(color.r+0.5,color.g+0.5);
    gray=max(color.b+0.5,gray);

    outColor = vec4(gray ,gray,gray,1.0);
}
