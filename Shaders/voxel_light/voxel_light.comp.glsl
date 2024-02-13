#version 450 core

#include "compiled.inc"
#include "std/math.glsl"
#include "std/gbuffer.glsl"
#include "std/shadows.glsl"
#include "std/imageatomic.glsl"

layout (local_size_x = 8, local_size_y = 8, local_size_z = 8) in;

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

uniform layout(rgba8) image3D voxelsOpac;
uniform layout(rgba8) image3D voxelsNor;
uniform layout(rgba8) image3D voxels;
#ifdef _EmissionShaded
uniform sampler2D gbuffer_emission;
#endif

#ifdef _ShadowMap
uniform sampler2DShadow shadowMap;
uniform sampler2DShadow shadowMapSpot;
uniform samplerCubeShadow shadowMapPoint;
#endif

void main() {
	const ivec3 src = ivec3(gl_GlobalInvocationID.xyz);

	for (int i = 0; i < 6; i++) {
		ivec3 dst = src;
		vec3 wposition = vec3(dst);
		dst.y += clipmapLevel * voxelgiResolution.x;
		dst.x += i * voxelgiResolution.x;

		vec4 col = imageLoad(voxelsOpac, dst);
		if (col.a == 0.0) return;

		//uint unor = imageLoad(voxelsNor, adjustedID).r;
		//vec3 wnormal = normalize(decNor(unor));

		//wposition -= wnormal * 0.01; // Offset

		float visibility;
		vec3 ld = lightPos - wposition;
		vec3 l;
		if (lightType == 0) { l = lightDir; visibility = 1.0; }
		else { l = normalize(ld); visibility = attenuate(distance(wposition, lightPos)); }

		//float dotNL = max(dot(wnormal, l), 0.0);
		//if (dotNL == 0.0) return;

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
			visibility *= texture(shadowMapPoint, vec4(-l, lpToDepth(ld, lightProj) - shadowsBias)).r;
		}
#endif

		if (lightType == 2) {
			float spotEffect = dot(lightDir, l);
			if (spotEffect < spotData.x) {
				visibility *= smoothstep(spotData.y, spotData.x, spotEffect);
			}
		}

		col.rgb *= visibility * lightColor;// * dotNL;
		col = clamp(col, vec4(0.0), vec4(1.0));

		#ifdef _EmissionShaded
		col += texture(gbuffer_emission, wposition.xy);
		#endif

		imageStore(voxels, dst, col);
	}
}
