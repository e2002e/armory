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

uniform writeonly image3D voxels;
#ifdef _VoxelGI
uniform layout(rgba8) image3D voxelsOut;
uniform sampler2D gbuffer1;
#else
uniform layout(r8) image3D voxelsOut;
#endif

uniform vec3 clipmap_center_last;
uniform int clipmapLevel;
uniform float voxelBlend;

layout (local_size_x = 8, local_size_y = 8, local_size_z = 8) in;

void main() {
	int res = voxelgiResolution.x;
	ivec3 src = ivec3(gl_GlobalInvocationID.xyz);
	ivec3 dst = src;

	#ifdef _VoxelGI
	vec4 col;
	const float voxelSize = 2.0 * voxelgiVoxelSize;
	vec3 wposition = vec3(dst);
	vec3 basecol = textureLod(gbuffer1, wposition.xy, 0.0).rgb;
	#else
	float opac;
	#endif

	dst.y += clipmapLevel * res;

	for (int i = 0; i < 6; i++)
	{
		#ifdef _VoxelGI
		col = vec4(0.0);
		#else
		opac = 0.0;
		#endif

		if (clipmap_center_last.x != 0 || clipmap_center_last.y != 0 || clipmap_center_last.z != 0)
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
				col = mix(imageLoad(voxelsOut, dst), vec4(basecol, 1.0), voxelBlend);
				#else
				opac = mix(imageLoad(voxelsOut, dst).r, 1.0, voxelBlend);
				#endif
		}
		else
			#ifdef _VoxelGI
			col = mix(imageLoad(voxelsOut, dst), vec4(basecol, 1.0), voxelBlend);
			#else
			opac = mix(imageLoad(voxelsOut, dst).r, 1.0, voxelBlend);
			#endif

		#ifdef _VoxelGI
		imageStore(voxels, dst, col);
		#else
		imageStore(voxels, dst, vec4(opac));
		#endif
	}
}
