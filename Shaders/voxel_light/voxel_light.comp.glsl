#version 450

layout (local_size_x = 4, local_size_y = 4, local_size_z = 4) in;

#include "compiled.inc"
#include "std/math.glsl"
#include "std/gbuffer.glsl"
#include "std/shadows.glsl"
#include "std/imageatomic.glsl"

uniform vec3 lightPos;
uniform vec3 lightColor;
uniform int lightType;
uniform vec3 lightDir;
uniform vec2 spotData;
uniform int clipmapLevel;
#ifdef _ShadowMap
uniform int lightShadow;
uniform vec2 lightProj;
uniform float shadowsBias;
uniform mat4 LVP;
#endif


uniform layout(binding = 0, rgba8) readonly image3D voxelsOpac;
// uniform layout(binding = 1, r32ui) readonly uimage3D voxelsNor;
// uniform layout(binding = 2, rgba8) writeonly image3D voxels;
uniform layout(binding = 1, rgba8) writeonly image3D voxels;
#ifdef _ShadowMap
uniform layout(binding = 2) sampler2D shadowMap;
uniform layout(binding = 3) samplerCube shadowMapCube;
#endif

void main() {
    // Adjust resolution based on clipmapLevel
    ivec3 adjustedID = ivec3(gl_GlobalInvocationID.xyz) / int(pow(2.0, float(clipmapLevel)));
    vec4 col = imageLoad(voxelsOpac, adjustedID);

    const vec3 hres = voxelgiResolution / 2;
    vec3 wposition = ((adjustedID - hres) / hres) * voxelgiResolution;

    float visibility;
    vec3 lp = lightPos - wposition;
    vec3 l;
    if (lightType == 0) { l = lightDir; visibility = 1.0; }
    else { l = normalize(lp); visibility = attenuate(distance(wposition, lightPos)); }

#ifdef _ShadowMap
    if (lightShadow == 1) {
        vec4 lightPosition = LVP * vec4(wposition, 1.0);
        vec3 lPos = lightPosition.xyz / lightPosition.w;
        if (texture(shadowMap, lPos.xy).r < lPos.z - shadowsBias) visibility = 0.0;
    }
    else if (lightShadow == 2) visibility *= float(texture(shadowMapCube, -l).r + shadowsBias > lpToDepth(lp, lightProj));
#endif

    if (lightType == 2) {
        float spotEffect = dot(lightDir, l);
        if (spotEffect < spotData.x) {
            visibility *= smoothstep(spotData.y, spotData.x, spotEffect);
        }
    }

    col.rgb += visibility * lightColor;
    col = clamp(col, vec4(0.0), vec4(1.0));

    imageStore(voxels, ivec3(gl_GlobalInvocationID.xyz), vec4(1.0));
}
