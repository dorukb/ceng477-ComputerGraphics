#version 410

in Data
{
    vec3 Position;
    vec3 Normal;
    vec2 TexCoord;
} data;
in vec3 LightVector;
in vec3 CameraVector;

uniform vec3 lightPosition;
uniform sampler2D TexColor;
uniform sampler2D MoonTexColor;
uniform sampler2D TexGrey;
uniform float textureOffset;

out vec4 FragColor;

vec3 ambientReflectenceCoefficient = vec3(0.5f,0.5f,0.5f);
vec3 ambientLightColor = vec3(0.6f, 0.6f, 0.6f);
vec3 specularReflectenceCoefficient= vec3(1.0f, 1.0f, 1.0f);
vec3 specularLightColor = vec3(1.0f,1.0f, 1.0f);
float SpecularExponent = 10;
vec3 diffuseReflectenceCoefficient; //
vec3 diffuseLightColor = vec3(1.0f,1.0f, 1.0f);

void main()
{
    // Calculate texture coordinate based on data.TexCoord

    
    vec4 textureColor = texture(MoonTexColor, data.TexCoord);
    
    diffuseReflectenceCoefficient = vec3(textureColor.xyz);
    vec3 ambient = (ambientReflectenceCoefficient * ambientLightColor).xyz;
    float theta = max(dot(data.Normal, LightVector), 0.0f);
    vec3 diffuse = theta * (diffuseLightColor * diffuseReflectenceCoefficient).xyz;
    vec3 reflected = reflect(-normalize(LightVector - 0.5), data.Normal);
    float alpha = max(dot(reflected, normalize(CameraVector + 5)), 0.0f);
    vec3 spec = pow(alpha,SpecularExponent) * (specularLightColor * specularReflectenceCoefficient).xyz;

    FragColor = vec4(ambient+diffuse+spec, 1.0f);
    //FragColor = vec4(clamp(FragColor.xyz, 0.0,1.0), 1.0);
    FragColor = vec4(clamp(textureColor.xyz * FragColor.xyz, 0.0,1.0), 1.0);


    //FragColor = vec4(1.0f, 1.0f, 1.0f, 1.0f);

}
