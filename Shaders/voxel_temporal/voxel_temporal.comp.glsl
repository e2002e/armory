#version 450

#include "compiled.inc"
#include "std/math.glsl"
#include "std/gbuffer.glsl"
#include "std/imageatomic.glsl"

uniform layout(r32ui) uimage3D voxelsOut;
uniform layout(r32ui) uimage3D voxels;
uniform layout(r32ui) uimage3D voxelsB;

uniform vec3 clipmap_center_last;
uniform int clipmapLevel;
uniform float voxelBlend;

layout (local_size_x = 8, local_size_y = 8, local_size_z = 8) in;

void main() {
	int res = voxelgiResolution.x / 8;
	ivec3 src = ivec3(gl_GlobalInvocationID.xyz);
	ivec3 dest = src;
	dest.y += clipmapLevel * res;

	vec4 opac = vec4(0.0);
	//for (int i = 0; i < 6; i++)
	{
		ivec3 coords = ivec3(dest - clipmap_center_last);
		//int aniso_face_start_x = i * res;
		//int aniso_face_end_x = aniso_face_start_x + res;
		int clipmap_face_start_y = clipmapLevel * res;
		int clipmap_face_end_y = clipmap_face_start_y + res;
		if (all(notEqual(clipmap_center_last, vec3(0.0)))) {
			if (
				coords.x >= 0 && coords.x < res &&
				coords.y >= clipmap_face_start_y && coords.y < clipmap_face_end_y &&
				coords.z >= 0 && coords.z < res
			)
				opac = mix(convRGBA8ToVec4(imageLoad(voxelsB, coords).r), convRGBA8ToVec4(imageLoad(voxels, coords).r),  voxelBlend);
		}
		else
			opac = mix(convRGBA8ToVec4(imageLoad(voxelsB, coords).r), convRGBA8ToVec4(imageLoad(voxels, coords).r),  voxelBlend);

		imageAtomicAdd(voxelsOut, dest, convVec4ToRGBA8(vec4(opac)));
	}
}
