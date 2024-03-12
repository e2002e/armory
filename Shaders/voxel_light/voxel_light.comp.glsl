#version 450

layout (local_size_x = 8, local_size_y = 8, local_size_z = 8) in;

#include "compiled.inc"
#include "std/math.glsl"
#include "std/gbuffer.glsl"
#include "std/shadows.glsl"
#include "std/imageatomic.glsl"
#include "std/conetrace.glsl"

uniform vec3 lightPos;
uniform vec3 lightColor;
uniform int lightType;
uniform vec3 lightDir;
uniform vec2 spotData;
#ifdef _ShadowMap
uniform int lightShadow;
uniform vec2 lightProj;
uniform float shadowsBias;
uniform mat4 LVP;
#endif

uniform vec3 clipmap_center;
uniform vec3 eye;
uniform int clipmapLevel;

uniform layout(r32ui) uimage3D voxels;
uniform layout(r32ui) uimage3D voxelsNor;
uniform sampler3D voxelsSampler;
uniform layout(r32ui) uimage3D voxelsEmission;
uniform layout(r32ui) uimage3D voxelsLight;
#ifdef _ShadowMap
uniform sampler2DShadow shadowMap;
uniform sampler2DShadow shadowMapSpot;
uniform samplerCubeShadow shadowMapPoint;
#endif

void main() {
	int res = voxelgiResolution.x;

	for (int i = 0; i < 6; i++)
	{
		ivec3 src = ivec3(gl_GlobalInvocationID.xyz);
		src.x += i * res;
		ivec3 dst = src;

		vec4 radiance = vec4(0.0);

		vec3 wposition = (gl_GlobalInvocationID.xyz + 0.5) / voxelgiResolution.x;
		wposition = wposition * 2.0 - 1.0;
		wposition *= voxelgiVoxelSize * pow(2.0, clipmapLevel);
		wposition *= voxelgiResolution.x;
		wposition += clipmap_center;

		radiance = convRGBA8ToVec4(imageLoad(voxels, src).r);
		vec4 emission = convRGBA8ToVec4(imageLoad(voxelsEmission, src).r);
		vec3 normal = decNor(imageLoad(voxelsNor, src).r);

		float visibility;
		vec3 lp = lightPos - wposition;
		vec3 l;
		if (lightType == 0) { l = lightDir; visibility = 1.0; }
		else { l = normalize(lp); visibility = attenuate(distance(wposition, lightPos)); }

		// float dotNL = max(dot(wnormal, l), 0.0);
		// if (dotNL == 0.0) return;

	#ifdef _ShadowMap
		if (lightShadow == 1) {
			vec4 lightPosition = LVP * vec4(wposition, 1.0);
			vec3 lPos = lightPosition.xyz / lightPosition.w;
			visibility = texture(shadowMap, vec3(lPos.xy, lPos.z - shadowsBias)).r;
		}
		else if (lightShadow == 2) {
			vec4 lightPosition = LVP * vec4(wposition, 1.0);
			vec3 lPos = lightPosition.xyz / lightPosition.w;
			visibility *= texture(shadowMapSpot, vec3(lPos.xy, lPos.z - shadowsBias)).r;
		}
		else if (lightShadow == 3) {
			visibility *= texture(shadowMapPoint, vec4(-l, lpToDepth(lp, lightProj) - shadowsBias)).r;
		}
	#endif

		if (lightType == 2) {
			float spotEffect = dot(lightDir, l);
			if (spotEffect < spotData.x) {
				visibility *= smoothstep(spotData.y, spotData.x, spotEffect);
			}
		}

		vec4 indirect = traceDiffuse(wposition, normal, voxelsSampler, eye);
		radiance.rgb *= (visibility * lightColor) / 3.1415 + indirect.rgb;
		radiance += emission;
		radiance = clamp(radiance, vec4(0.0), vec4(1.0));

		imageAtomicMax(voxelsLight, dst, convVec4ToRGBA8(radiance));
	}
}
