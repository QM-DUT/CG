
#version 330
in vec2 TexCoord;
uniform sampler2D textureMap;
out vec4 outColor;
void main()
{
    float core[] = float[]( 0.0144,0.0281,0.0351,0.0281,0.0144,
                            0.0281,0.0547,0.0683,0.0547,0.0281,
                            0.0351,0.0683,0.0853,0.0683,0.0351,
                            0.0281,0.0547,0.0683,0.0547,0.0281,
                            0.0144,0.0281,0.0351,0.0281,0.0144);
    float factor = 1.0;

    ivec2 texSize = textureSize(textureMap, 0);
    float s = float(texSize.s);
    float t = float(texSize.t);

    vec3 color = vec3(0., 0., 0.);

    vec3 Colors[]=vec3[](
    texture(textureMap, TexCoord+vec2(0.0,-2.0/t)+vec2(-2.0/s,0.0)).rgb,
    texture(textureMap, TexCoord+vec2(0.0,-2.0/t)+vec2(-1.0/s,0.0)).rgb,
    texture(textureMap, TexCoord+vec2(0.0,-2.0/t)+vec2(0.0/s,0.0)).rgb,
    texture(textureMap, TexCoord+vec2(0.0,-2.0/t)+vec2(1.0/s,0.0)).rgb,
    texture(textureMap, TexCoord+vec2(0.0,-2.0/t)+vec2(2.0/s,0.0)).rgb,

    texture(textureMap, TexCoord+vec2(0.0,-1.0/t)+vec2(-2.0/s,0.0)).rgb,
    texture(textureMap, TexCoord+vec2(0.0,-1.0/t)+vec2(-1.0/s,0.0)).rgb,
    texture(textureMap, TexCoord+vec2(0.0,-1.0/t)+vec2(0.0/s,0.0)).rgb,
    texture(textureMap, TexCoord+vec2(0.0,-1.0/t)+vec2(1.0/s,0.0)).rgb,
    texture(textureMap, TexCoord+vec2(0.0,-1.0/t)+vec2(2.0/s,0.0)).rgb,

    texture(textureMap, TexCoord+vec2(0.0,0.0/t)+vec2(-2.0/s,0.0)).rgb,
    texture(textureMap, TexCoord+vec2(0.0,0.0/t)+vec2(-1.0/s,0.0)).rgb,
    texture(textureMap, TexCoord+vec2(0.0,0.0/t)+vec2(0.0/s,0.0)).rgb,
    texture(textureMap, TexCoord+vec2(0.0,0.0/t)+vec2(1.0/s,0.0)).rgb,
    texture(textureMap, TexCoord+vec2(0.0,0.0/t)+vec2(2.0/s,0.0)).rgb,

    texture(textureMap, TexCoord+vec2(0.0,1.0/t)+vec2(-2.0/s,0.0)).rgb,
    texture(textureMap, TexCoord+vec2(0.0,1.0/t)+vec2(-1.0/s,0.0)).rgb,
    texture(textureMap, TexCoord+vec2(0.0,1.0/t)+vec2(0.0/s,0.0)).rgb,
    texture(textureMap, TexCoord+vec2(0.0,1.0/t)+vec2(1.0/s,0.0)).rgb,
    texture(textureMap, TexCoord+vec2(0.0,1.0/t)+vec2(2.0/s,0.0)).rgb,

    texture(textureMap, TexCoord+vec2(0.0,2.0/t)+vec2(-2.0/s,0.0)).rgb,
    texture(textureMap, TexCoord+vec2(0.0,2.0/t)+vec2(-1.0/s,0.0)).rgb,
    texture(textureMap, TexCoord+vec2(0.0,2.0/t)+vec2(0.0/s,0.0)).rgb,
    texture(textureMap, TexCoord+vec2(0.0,2.0/t)+vec2(1.0/s,0.0)).rgb,
    texture(textureMap, TexCoord+vec2(0.0,2.0/t)+vec2(2.0/s,0.0)).rgb
    );
    for(int i=0;i<25;i++)
    {
        color +=Colors[i] * core[i];
    }
    color /= factor;
    outColor = vec4(color, 1.0);
}