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

uniform writeonly image3D voxelsOut;
#ifdef _VoxelGI
uniform layout(rgba8) image3D voxels;
#else
uniform layout(r8) image3D voxels;
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
	#else
	float opac;
	#endif

	dst.y += clipmapLevel * res;

	for (int i = 0; i < 6; i++)
	{
		dst.x += i * res;
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
				col = imageLoad(voxels, coords);
				#else
				opac = imageLoad(voxels, coords).r;
				#endif
			else
				#ifdef _VoxelGI
				col = vec4(0.0);
				#else
				opac = 0.0;
				#endif
		}
		else
			#ifdef _VoxelGI
			col = imageLoad(voxels, dst);
			#else
			opac = imageLoad(voxels, dst).r;
			#endif

		#ifdef _VoxelGI
		imageStore(voxelsOut, dst, col);
		#else
		imageStore(voxelsOut, dst, vec4(opac));
		#endif
	}
}
