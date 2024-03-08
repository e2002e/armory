#version 450

// layout (local_size_x = 4, local_size_y = 4, local_size_z = 4) in;
layout (local_size_x = 8, local_size_y = 8, local_size_z = 8) in;

#include "compiled.inc"
#include "std/math.glsl"
#include "std/gbuffer.glsl"
#include "std/shadows.glsl"
#include "std/imageatomic.glsl"
#include "std/voxels_constants.h"

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
uniform vec3 clipmap_center_last;
uniform int clipmapLevel;

uniform layout(r32ui) uimage3D voxels;
uniform layout(r32ui) uimage3D voxelsB;
uniform layout(r32ui) uimage3D voxelsEmission;
uniform layout(r32ui) uimage3D voxelsLight;
#ifdef _ShadowMap
uniform sampler2DShadow shadowMap;
uniform sampler2DShadow shadowMapSpot;
uniform samplerCubeShadow shadowMapPoint;
#endif

void main() {
	int res = voxelgiResolution.x;
	vec4 aniso_colors[6];

	for (int i = 0; i < 6 + DIFFUSE_CONE_COUNT; i++)
	{
		ivec3 src = ivec3(gl_GlobalInvocationID.xyz);
		src.x += i * res;
		ivec3 dst = src;
		dst.y += clipmapLevel * res;

		vec4 radiance = vec4(0.0);

		if (i < 6)
		{
			vec3 wposition = (gl_GlobalInvocationID.xyz + 0.5) / voxelgiResolution.x;
			wposition = wposition * 2.0 - 1.0;
			wposition *= voxelgiVoxelSize * pow(2.0, clipmapLevel);
			wposition *= voxelgiResolution.x;
			wposition += clipmap_center;

			radiance = convRGBA8ToVec4(imageLoad(voxels, src).r);
			vec4 emission = convRGBA8ToVec4(imageLoad(voxelsEmission, src).r);

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

			radiance.rgb *= visibility * lightColor;// * dotNL;

			if (radiance.a > 0.0)
			{
				if (any(notEqual(clipmap_center_last, vec3(0.0))))
				{
					ivec3 coords = ivec3(dst - clipmap_center_last);
					int aniso_face_start_x = i * res;
					int aniso_face_end_x = aniso_face_start_x + res;
					int clipmap_face_start_y = clipmapLevel * res;
					int clipmap_face_end_y = clipmap_face_start_y + res;
					if (
						coords.x >= aniso_face_start_x && coords.x < aniso_face_end_x &&
						coords.y >= clipmap_face_start_y && coords.y < clipmap_face_end_y &&
						coords.z >= 0 && coords.z < res
					)
						radiance = mix(convRGBA8ToVec4(imageLoad(voxelsB, dst).r), radiance, 0.5);
				}
				else
					radiance = mix(convRGBA8ToVec4(imageLoad(voxelsB, dst).r), radiance, 0.5);
			}
			else
				radiance = vec4(0.0);

			radiance = clamp(radiance + emission, vec4(0.0), vec4(1.0));
			aniso_colors[i] = radiance;
		}
		else {
			// precompute cone sampling:
			vec3 coneDirection = DIFFUSE_CONE_DIRECTIONS[i - 6];
			vec3 aniso_direction = -coneDirection;
			uvec3 face_offsets = uvec3(
				aniso_direction.x > 0 ? 0 : 1,
				aniso_direction.y > 0 ? 2 : 3,
				aniso_direction.z > 0 ? 4 : 5
			);
			vec3 direction_weights = abs(coneDirection);
			vec4 sam =
				aniso_colors[face_offsets.x] * direction_weights.x +
				aniso_colors[face_offsets.y] * direction_weights.y +
				aniso_colors[face_offsets.z] * direction_weights.z
				;
			radiance = sam;
		}

		imageAtomicAdd(voxelsLight, dst, convVec4ToRGBA8(radiance));
	}
}
