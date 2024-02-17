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
#version 450

#include "compiled.inc"
#include "std/math.glsl"
#include "std/gbuffer.glsl"
#include "std/imageatomic.glsl"
#include "std/voxels_constants.h"

#ifdef _VoxelGI
uniform layout(r32ui) uimage3D voxels;
uniform layout(r32ui) uimage3D voxelsOut;
uniform layout(r32ui) uimage3D voxelsOutB;
#else
uniform layout(r8) image3D voxels;
uniform layout(r8) image3D voxelsOut;
uniform layout(r8) image3D voxelsOutB;
#endif

uniform vec3 clipmap_center_last;
uniform int clipmapLevel;
uniform float voxelBlend;

layout (local_size_x = 8, local_size_y = 8, local_size_z = 8) in;

void main() {
	int res = voxelgiResolution.x;
	ivec3 src = ivec3(gl_GlobalInvocationID.xyz);
	#ifdef _VoxelGI
	vec4 col;
	vec4 aniso_colors[6];
	const float voxelSize = 2.0 * voxelgiVoxelSize;
	#else
	float opac;
	float aniso_colors[6];
	#endif

	for (int i = 0; i < 6 + 12; i++)
	{
		ivec3 dst = src;
		dst.x += i * res;
		dst.y += clipmapLevel * res;
		#ifdef _VoxelGI
		col = vec4(0.0);
		#else
		opac = 0.0;
		#endif

		if (i < 6) {
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
					#ifdef _VoxelGI
					col = mix(convRGBA8ToVec4(imageLoad(voxelsOutB, dst).r), convRGBA8ToVec4(imageLoad(voxelsOut, dst).r), voxelBlend);
					#else
					opac = mix(imageLoad(voxelsOutB, dst).r, imageLoad(voxelsOut, dst).r, voxelBlend);
					#endif
			}
			else
				#ifdef _VoxelGI
				col = mix(convRGBA8ToVec4(imageLoad(voxelsOutB, dst).r), convRGBA8ToVec4(imageLoad(voxelsOut, dst).r), voxelBlend);
				#else
				opac = mix(imageLoad(voxelsOutB, dst).r, imageLoad(voxelsOut, dst).r, voxelBlend);
				#endif
			#ifdef _VoxelGI
			aniso_colors[i] = col;
			#else
			aniso_colors[i] = opac;
			#endif
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
			#ifdef _VoxelGI
			vec4 sam =
				aniso_colors[face_offsets.x] * direction_weights.x +
				aniso_colors[face_offsets.y] * direction_weights.y +
				aniso_colors[face_offsets.z] * direction_weights.z
				;
			col = sam;
			#else
			float sam =
				aniso_colors[face_offsets.x] * direction_weights.x +
				aniso_colors[face_offsets.y] * direction_weights.y +
				aniso_colors[face_offsets.z] * direction_weights.z
				;
			opac = sam;
			#endif
		}
		#ifdef _VoxelGI
		imageAtomicAdd(voxels, dst, convVec4ToRGBA8(col));
		#else
		imageStore(voxels, dst, vec4(opac));
		#endif
	}
}
