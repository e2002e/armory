/*
Copyright (c) 2024 Turánszki János

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 */
#ifndef _CONETRACE_GLSL_
#define _CONETRACE_GLSL_

#include "std/voxels_constants.h"

// References
// https://github.com/Friduric/voxel-cone-tracing
// https://github.com/Cigg/Voxel-Cone-Tracing
// https://github.com/GreatBlambo/voxel_cone_tracing/
// http://simonstechblog.blogspot.com/2013/01/implementing-voxel-cone-tracing.html
// http://leifnode.com/2015/05/voxel-cone-traced-global-illumination/
// http://www.seas.upenn.edu/%7Epcozzi/OpenGLInsights/OpenGLInsights-SparseVoxelization.pdf
// https://research.nvidia.com/sites/default/files/publications/GIVoxels-pg2011-authors.pdf

const float MAX_DISTANCE = 1.73205080757 * voxelgiRange;

// uniform sampler3D voxels;
// uniform sampler3D voxelsLast;

// vec3 orthogonal(const vec3 u) {
// 	// Pass normalized u
// 	const vec3 v = vec3(0.99146, 0.11664, 0.05832); // Pick any normalized vector
// 	return abs(dot(u, v)) > 0.99999 ? cross(u, vec3(0.0, 1.0, 0.0)) : cross(u, v);
// }

vec3 tangent(const vec3 n) {
	vec3 t1 = cross(n, vec3(0, 0, 1));
	vec3 t2 = cross(n, vec3(0, 1, 0));
	if (length(t1) > length(t2)) return normalize(t1);
	else return normalize(t2);
}

vec3 faceIndices(const vec3 dir) {
 	vec3 ret;
 	ret.x = (dir.x < 0.0) ? 0 : 1;
 	ret.y = (dir.y < 0.0) ? 2 : 3;
 	ret.z = (dir.z < 0.0) ? 4 : 5;
 	return ret / (6 + DIFFUSE_CONE_COUNT);
}

#ifdef _VoxelGI
vec4 sampleVoxel(vec3 P, sampler3D voxels, vec3 dir, const vec3 indices, const int precomputed_direction, const vec3 clipmap_center, const float clipmap_index, const float lod) {
 	dir = abs(dir);
 	vec4 col = vec4(0.0);
	float voxelSize = pow(2.0, clipmap_index) * voxelgiVoxelSize;
 	vec3 tc = (P - clipmap_center) / (voxelSize * voxelgiResolution.x);
	tc = tc * 0.5 + 0.5;
	float half_texel = 0.5 / voxelgiResolution.x;
	tc = clamp(tc, half_texel, 1.0 - half_texel);
	tc.x = (tc.x + precomputed_direction) / (6 + DIFFUSE_CONE_COUNT);
	tc.y = (tc.y + clipmap_index) / voxelgiClipmapCount;

	if (precomputed_direction == 0) {
		col += dir.x * textureLod(voxels, vec3(tc.x + indices.x, tc.y, tc.z), lod)
			+ dir.y * textureLod(voxels, vec3(tc.x + indices.y, tc.y, tc.z), lod)
			+ dir.z * textureLod(voxels, vec3(tc.x + indices.z, tc.y, tc.z), lod);
	}
	else
		col = textureLod(voxels, tc, lod);

	return col;
}
#else
float sampleVoxel(vec3 P, sampler3D voxels, vec3 dir, const vec3 indices, const int precomputed_direction, const vec3 clipmap_center, const float clipmap_index, const float lod) {
 	dir = abs(dir);
 	float opac = 0.0;
	float voxelSize = pow(2.0, clipmap_index) * voxelgiVoxelSize;
 	vec3 tc = (P - clipmap_center) / (voxelSize * voxelgiResolution.x);
	tc = tc * 0.5 + 0.5;
	float half_texel = 0.5 / voxelgiResolution.x;
	tc = clamp(tc, half_texel, 1.0 - half_texel);
	tc.x = (tc.x + precomputed_direction) / (6 + DIFFUSE_CONE_COUNT);
	tc.y = (tc.y + clipmap_index) / voxelgiClipmapCount;

	if (precomputed_direction == 0) {
		opac += dir.x * textureLod(voxels, vec3(tc.x + indices.x, tc.y, tc.z), lod).r
			+ dir.y * textureLod(voxels, vec3(tc.x + indices.y, tc.y, tc.z), lod).r
			+ dir.z * textureLod(voxels, vec3(tc.x + indices.z, tc.y, tc.z), lod).r;
	}
	else
		opac = textureLod(voxels, tc, lod).r;

	return opac;
}
#endif


#ifdef _VoxelGI
vec4 traceCone(sampler3D voxels, vec3 origin, vec3 n, vec3 dir, const int precomputed_direction, const float aperture, const float maxDist, const vec3 clipmap_center) {
    vec4 sampleCol = vec4(0.0);
	vec3 indices = faceIndices(dir);
	float voxelSize0 = voxelgiVoxelSize * 2.0 * voxelgiOffset;
	float dist = voxelSize0;
	float step_dist = dist;
	vec3 samplePos;
	vec3 start_pos = origin + n * voxelSize0;
	float coneCoefficient = 2.0 * tan(aperture * 0.5);
	int clipmap_index0 = 0;

    while (sampleCol.a < 1.0 && dist < maxDist && clipmap_index0 < voxelgiClipmapCount) {
		vec4 mipSample = vec4(0.0);
		float diam = max(voxelSize0, dist * coneCoefficient);
        float lod = clamp(log2(diam / voxelSize0), clipmap_index0, voxelgiClipmapCount - 1);

        float clipmap_index = floor(lod);
		float clipmap_blend = fract(lod);

		float voxelSize = pow(2.0, clipmap_index) * voxelgiVoxelSize;
		vec3 p0 = start_pos + dir * dist;

        samplePos = (p0 - clipmap_center) / (voxelSize * voxelgiResolution.x);
		samplePos = samplePos * 0.5 + 0.5;

		if (any(notEqual(samplePos, clamp(samplePos, 0.0, 1.0)))) {
			clipmap_index0++;
			continue;
		}

		mipSample = sampleVoxel(p0, voxels, dir, indices, precomputed_direction, clipmap_center, clipmap_index, 0.0);

		if(clipmap_blend > 0.0 && clipmap_index < voxelgiClipmapCount - 1) {
			vec4 mipSampleNext = sampleVoxel(p0, voxels, dir, indices, precomputed_direction, clipmap_center, clipmap_index + 1.0, 0.0);
			mipSample = mix(mipSample, mipSampleNext, clipmap_blend);
		}

		mipSample *= step_dist / voxelSize;
		sampleCol += (1.0 - sampleCol) * mipSample;
		step_dist = diam * voxelgiStep;
		dist += step_dist;
	}
    return sampleCol;
}


vec4 traceDiffuse(const vec3 origin, const vec3 normal, sampler3D voxels, const vec3 clipmap_center) {
	vec4 sampleCol = vec4(0.0);
	float sum = 0.0;
	vec4 amount = vec4(0.0);
	for (int i = 0; i < DIFFUSE_CONE_COUNT; i++) {
		vec3 coneDir = DIFFUSE_CONE_DIRECTIONS[i];
		int precomputed_direction = 6 + i;
		const float cosTheta = dot(normal, coneDir);
		if (cosTheta <= 0)
			continue;
		amount += traceCone(voxels, origin, normal, coneDir, precomputed_direction, DIFFUSE_CONE_APERTURE, MAX_DISTANCE, clipmap_center);
		sum += cosTheta;
	}
	amount /= sum;
	sampleCol = max(vec4(0.0), amount);
	return sampleCol * voxelgiOcc;;
}


vec4 traceSpecular(const vec3 origin, const vec3 normal, sampler3D voxels, const vec3 viewDir, const float roughness, const vec3 clipmap_center) {
	float specularAperture = clamp(tan((3.14159265 / 2) * roughness), 0.0174533 * 3.0, 3.14159265);
	vec3 specularDir = reflect(viewDir, normal);
	return (traceCone(voxels, origin, normal, specularDir, 0, specularAperture, MAX_DISTANCE, clipmap_center) + traceCone(voxels, origin, normal, specularDir, 0, specularAperture / 3.14159265, MAX_DISTANCE, clipmap_center) / 2.0) * voxelgiOcc;
}


vec3 traceRefraction(const vec3 origin, const vec3 normal, sampler3D voxels, const vec3 viewDir, const float ior, const float roughness, const vec3 clipmap_center) {
 	const float transmittance = 1.0;
 	vec3 refractionDir = refract(viewDir, normal, 1.0 / ior);
 	float refractiveAperture = clamp(tan((3.14159265 / 2) * roughness), 0.0174533 * 3.0, 3.14159265);
 	return transmittance * traceCone(voxels, origin, normal, refractionDir, 0, refractiveAperture, MAX_DISTANCE, clipmap_center).xyz * voxelgiOcc;
}
#endif


#ifdef _VoxelAOvar
float traceConeAO(sampler3D voxels, vec3 origin, vec3 n, vec3 dir, const int precomputed_direction, const float aperture, const float maxDist, const vec3 clipmap_center) {
    float sampleCol = 0.0;;
	vec3 indices = faceIndices(dir);
	float voxelSize0 = voxelgiVoxelSize * 2.0 * voxelgiOffset;
	float dist = voxelSize0;
	float step_dist = dist;
	vec3 samplePos;
	vec3 start_pos = origin + n * voxelSize0;
	float coneCoefficient = 2.0 * tan(aperture * 0.5);
	int clipmap_index0 = 0;

    while (sampleCol < 1.0 && dist < maxDist && clipmap_index0 < voxelgiClipmapCount) {
		float mipSample = 0.0;
		float diam = max(voxelSize0, dist * coneCoefficient);
        float lod = clamp(log2(diam / voxelSize0), clipmap_index0, voxelgiClipmapCount - 1);
		float clipmap_index = floor(lod);
		float clipmap_blend = fract(lod);
		vec3 p0 = start_pos + dir * dist;
		float voxelSize = pow(2.0, clipmap_index) * voxelgiVoxelSize;

        samplePos = (p0 - clipmap_center) / (voxelSize * voxelgiResolution.x);
		samplePos = samplePos * 0.5 + 0.5;

		if ((any(notEqual(clamp(samplePos, 0.0, 1.0), samplePos)))) {
			clipmap_index0++;
			continue;
		}

		mipSample = sampleVoxel(p0, voxels, dir, indices, precomputed_direction, clipmap_center, clipmap_index, 0.0);

		if(clipmap_blend > 0.0 && clipmap_index < voxelgiClipmapCount - 1) {
			float mipSampleNext = sampleVoxel(p0, voxels, dir, indices, precomputed_direction, clipmap_center, clipmap_index + 1.0, 0.0);
			mipSample = mix(mipSample, mipSampleNext, clipmap_blend);
		}

		mipSample *= step_dist / voxelSize;
		sampleCol += (1.0 - sampleCol) * mipSample;
		step_dist = diam * voxelgiStep;
		dist += step_dist;
	}
    return sampleCol;
}


float traceAO(const vec3 origin, const vec3 normal, sampler3D voxels, const vec3 clipmap_center) {
	float sampleCol = 0.0;
	float sum = 0.0;
	float amount = 0.0;
	for (int i = 0; i < DIFFUSE_CONE_COUNT; i++) {
		vec3 coneDir = DIFFUSE_CONE_DIRECTIONS[i];
		int precomputed_direction = 6 + i;
		const float cosTheta = dot(normal, coneDir);
		if (cosTheta <= 0)
			continue;
		amount += traceConeAO(voxels, origin, normal, coneDir, precomputed_direction, DIFFUSE_CONE_APERTURE, MAX_DISTANCE, clipmap_center);
		sum += cosTheta;
	}
	amount /= sum;
	sampleCol = max(0.0, amount);
	return sampleCol * voxelgiOcc;
}
#endif


#ifdef _VoxelShadow
float traceConeShadow(sampler3D voxels, const vec3 origin, vec3 n, vec3 dir, const float aperture, const float maxDist, const vec3 clipmap_center) {
    dir = normalize(dir);
    uvec3 indices = faceIndices(dir);
    float sampleCol = 0.0;
	float voxelSize0 = voxelgiVoxelSize * 2.0 * voxelgiOffset;
	float dist = voxelSize0;
	float step_dist = dist;
	vec3 samplePos;
	vec3 start_pos = origin + n * voxelSize0;
	float coneCoefficient = 2.0 * tan(aperture * 0.5);
	int clipmap_index0 = 0;

    while (sampleCol < 1.0 && dist < maxDist && clipmap_index0 < voxelgiClipmapCount) {
		float mipSample = 0.0;
		float diam = max(voxelSize0, dist * coneCoefficient);
        float lod = clamp(log2(diam / voxelSize0), clipmap_index0, voxelgiClipmapCount - 1);
		float clipmap_index = floor(lod);
		float clipmap_blend = fract(lod);
		vec3 p0 = start_pos + dir * dist;
		float voxelSize = pow(2.0, clipmap_index) * voxelgiVoxelSize;

        samplePos = ((start_pos + dir * dist) - clipmap_center) / (voxelSize * voxelgiResolution.x);
		samplePos = samplePos * 0.5 + 0.5;
		if ((any(notEqual(samplePos, clamp(samplePos, 0.0, 1.0))))) {
			clipmap_index0++;
			continue;
		}
		samplePos.y = (samplePos.y + clipmap_index) / voxelgiClipmapCount;

		#ifdef _VoxelAOvar
		mipSample = sampleVoxel(p0, voxels, samplePos, dir, indices, clipmap_center, clipmap_index, 0.0);
		#else
		mipSample = sampleVoxel(p0, voxels, samplePos, dir, indices, clipmap_center, clipmap_index, 0.0).a;
		#endif

		if(clipmap_blend > 0.0 && clipmap_index < voxelgiClipmapCount - 1) {
			#ifdef _VoxelAOvar
			float mipSampleNext = sampleVoxel(p0, voxels, samplePos, dir, indices, clipmap_center, clipmap_index + 1.0, 0.0);
			#else
			float mipSampleNext = sampleVoxel(p0, voxels, samplePos, dir, indices, clipmap_center, clipmap_index + 1.0, 0.0).a;
			#endif
			mipSample = mix(mipSample, mipSampleNext, clipmap_blend);
		}

		mipSample *= step_dist / voxelSize;
		sampleCol += (1.0 - sampleCol) * mipSample;
		step_dist = diam * voxelgiStep;
		dist += step_dist;
	}
	return sampleCol;
}


float traceShadow(const vec3 origin, const vec3 normal, sampler3D voxels, const vec3 dir, const vec3 clipmap_center) {
	return traceConeShadow(voxels, origin, normal, dir, voxelgiAperture, MAX_DISTANCE, clipmap_center) * voxelgiOcc;
}
#endif
#endif // _CONETRACE_GLSL_
