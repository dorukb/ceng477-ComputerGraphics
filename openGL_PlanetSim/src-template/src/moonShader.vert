#version 410

layout (location = 0) in vec3 VertexPosition;
layout (location = 1) in vec3 VertexNormal;
layout (location = 2) in vec2 VertexTex;

uniform vec3 lightPosition;
uniform vec3 cameraPosition;

uniform mat4 ProjectionMatrix;
uniform mat4 ViewMatrix;
uniform mat4 NormalMatrix;
uniform mat4 MVP;
uniform mat4 moonMVP;

uniform sampler2D TexColor;
uniform sampler2D TexGrey;
uniform float textureOffset;
uniform float orbitDegree;

uniform float heightFactor;
uniform float imageWidth;
uniform float imageHeight;


out Data
{
    vec3 Position;
    vec3 Normal;
    vec2 TexCoord;
} data;


out vec3 LightVector;// Vector from Vertex to Light;
out vec3 CameraVector;// Vector from Vertex to Camera;


void main()
{
    vec2 textureCoordinate = VertexTex;
    vec3 calculated_pos = VertexPosition;

    CameraVector = normalize(cameraPosition - calculated_pos);
    LightVector = normalize(lightPosition - calculated_pos);
    gl_Position = moonMVP * vec4(calculated_pos.xyz, 1.0f);


    // Calculate texture coordinate based on data.TexCoord
    //vec2 textureCoordinate = vec2(0, 0);
   // vec4 texColor = texture(TexGrey, textureCoordinate);
    data.TexCoord = VertexTex;
    data.Position = gl_Position.xyz;
    data.Normal = VertexNormal;

    // get texture value, compute height
    // compute normal vector


   // set gl_Position variable correctly to give the transformed vertex position
    //gl_Position = vec4(VertexPosition.x, VertexPosition.y, VertexPosition.z, 1.0);

}