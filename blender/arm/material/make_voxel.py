import bpy

import arm.utils
import arm.assets as assets
import arm.material.cycles as cycles
import arm.material.mat_state as mat_state
import arm.material.mat_utils as mat_utils
import arm.material.make_particle as make_particle
import arm.make_state as state

import arm.utils
import arm.assets as assets
import arm.material.mat_state as mat_state

if arm.is_reload(__name__):
    arm.utils = arm.reload_module(arm.utils)
    assets = arm.reload_module(assets)
    mat_state = arm.reload_module(mat_state)
else:
    arm.enable_reload(__name__)

def make(context_id):
    rpdat = arm.utils.get_rp()
    if rpdat.rp_voxels == 'Voxel GI':
        con = make_gi(context_id)
    else:
        con = make_ao(context_id)

    assets.vs_equal(con, assets.shader_cons['voxel_vert'])
    assets.fs_equal(con, assets.shader_cons['voxel_frag'])
    assets.gs_equal(con, assets.shader_cons['voxel_geom'])

    return con

def make_gi(context_id):
    con_voxel = mat_state.data.add_context({ 'name': context_id, 'depth_write': False, 'compare_mode': 'always', 'cull_mode': 'none', 'color_write_red': False, 'color_write_green': False, 'color_write_blue': False, 'color_write_alpha': False, 'conservative_raster': True })
    wrd = bpy.data.worlds['Arm']

    vert = con_voxel.make_vert()
    frag = con_voxel.make_frag()
    geom = con_voxel.make_geom()
    tesc = None
    tese = None
    geom.ins = vert.outs
    frag.ins = geom.outs

    vert.add_include('compiled.inc')
    geom.add_include('compiled.inc')
    frag.add_include('compiled.inc')
    frag.add_include('std/math.glsl')
    frag.add_include('std/imageatomic.glsl')
    frag.add_include('std/gbuffer.glsl')
    frag.add_include('std/brdf.glsl')

    rpdat = arm.utils.get_rp()
    frag.add_uniform('writeonly layout(rgba8) image3D voxels')
    frag.add_uniform('writeonly layout(rgba8) image3D voxelsNor')

    frag.write('vec3 wposition;')#dummy
    frag.write('vec3 basecol;')
    frag.write('float roughness;') #
    frag.write('float metallic;') #
    frag.write('float occlusion;') #
    frag.write('float specular;') #
    frag.write('vec3 emissionCol = vec3(0.0);')
    parse_opacity = rpdat.arm_voxelgi_refraction
    if parse_opacity:
        frag.write('float opacity;')
        frag.write('float ior;')

    frag.write('float dotNV = 0.0;')
    cycles.parse(mat_state.nodes, con_voxel, vert, frag, geom, tesc, tese, parse_opacity=False, parse_displacement=False, basecol_only=True)

    # Voxelized particles
    particle = mat_state.material.arm_particle_flag
    if particle and rpdat.arm_particles == 'On':
        # make_particle.write(vert, particle_info=cycles.particle_info)
        frag.write_pre = True
        frag.write('const float p_index = 0;')
        frag.write('const float p_age = 0;')
        frag.write('const float p_lifetime = 0;')
        frag.write('const vec3 p_location = vec3(0);')
        frag.write('const float p_size = 0;')
        frag.write('const vec3 p_velocity = vec3(0);')
        frag.write('const vec3 p_angular_velocity = vec3(0);')
        frag.write_pre = False

    if not frag.contains('vec3 n ='):
        frag.write_pre = True
        frag.write('vec3 n;')
        frag.write_pre = False

    export_mpos = frag.contains('mposition') and not frag.contains('vec3 mposition')
    if export_mpos:
        vert.add_out('vec3 mpositionGeom')
        vert.write_pre = True
        vert.write('mpositionGeom = pos;')
        vert.write_pre = False

    export_bpos = frag.contains('bposition') and not frag.contains('vec3 bposition')
    if export_bpos:
        vert.add_out('vec3 bpositionGeom')
        vert.add_uniform('vec3 dim', link='_dim')
        vert.add_uniform('vec3 hdim', link='_halfDim')
        vert.write_pre = True
        vert.write('bpositionGeom = (pos.xyz + hdim) / dim;')
        vert.write_pre = False

    vert.add_uniform('mat4 W', '_worldMatrix')
    vert.add_uniform('mat3 N', '_normalMatrix')
    vert.add_out('vec3 voxpositionGeom')
    vert.add_out('vec3 voxnormalGeom')

    if con_voxel.is_elem('col'):
        vert.add_out('vec3 vcolorGeom')
        vert.write('vcolorGeom = col.rgb;')

    if con_voxel.is_elem('tex'):
        vert.add_out('vec2 texCoordGeom')
        vert.write('texCoordGeom = tex;')

    vert.add_uniform('vec3 clipmap_center', '_voxelSize')
    vert.add_uniform('int clipmapLevel', '_clipmapLevel')

    vert.write('vec3 P = vec3(W * vec4(pos.xyz, 1.0));')
    vert.write('float voxelSize = voxelgiVoxelSize * pow(2.0, clipmapLevel);')
    vert.write('float texelSize = 2.0 * voxelSize;')
    vert.write('vec3 clipmap_center = floor(eye / texelSize) * texelSize;')
    vert.write('voxpositionGeom = (P - clipmap_center) / voxelSize * 1.0 / voxelgiResolution.x;')
    vert.write('voxnormalGeom = normalize(N * vec3(nor.xy, pos.w));')

    geom.add_out('vec3 voxposition')
    geom.add_out('vec3 voxnormal')

    if con_voxel.is_elem('col'):
        geom.add_out('vec3 vcolor')
    if con_voxel.is_elem('tex'):
        geom.add_out('vec2 texCoord')
    if export_mpos:
        geom.add_out('vec3 mposition')
    if export_bpos:
        geom.add_out('vec3 bposition')

    geom.write('vec3 p1 = voxpositionGeom[1] - voxpositionGeom[0];')
    geom.write('vec3 p2 = voxpositionGeom[2] - voxpositionGeom[0];')

    geom.write('vec3 p = abs(cross(p1, p2));')
    geom.write('for (uint i = 0; i < 3; ++i) {')
    geom.write('    voxposition = voxpositionGeom[i];')
    geom.write('    voxnormal = voxnormalGeom[i];')
    if con_voxel.is_elem('col'):
        geom.write('    vcolor = vcolorGeom[i];')
    if con_voxel.is_elem('tex'):
        geom.write('    texCoord = texCoordGeom[i];')
    if export_mpos:
        geom.write('    mposition = mpositionGeom[i];')
    if export_bpos:
        geom.write('    bposition = bpositionGeom[i];')
    geom.write('    if (p.z > p.x && p.z > p.y) {')
    geom.write('        gl_Position = vec4(voxposition.x, voxposition.y, 0.0, 1.0);')
    geom.write('    }')
    geom.write('    else if (p.x > p.y && p.x > p.z) {')
    geom.write('        gl_Position = vec4(voxposition.y, voxposition.z, 0.0, 1.0);')
    geom.write('    }')
    geom.write('    else {')
    geom.write('        gl_Position = vec4(voxposition.x, voxposition.z, 0.0, 1.0);')
    geom.write('    }')
    geom.write('    EmitVertex();')
    geom.write('}')
    geom.write('EndPrimitive();')

    frag.write('if (abs(voxposition.z) > ' + rpdat.rp_voxelgi_resolution_z + ' || abs(voxposition.x) > 1 || abs(voxposition.y) > 1) return;')

    frag.add_uniform('int clipmapLevel', '_clipmapLevel')
    frag.write('vec3 uvw = (voxposition * 0.5 + 0.5) * voxelgiResolution.x;')
    frag.write('uvw.y += clipmapLevel * voxelgiResolution.x;')
    frag.write('vec3 face_offsets = vec3(')
    frag.write('	voxnormal.x > 0 ? 0 : 1,')
    frag.write('	voxnormal.y > 0 ? 2 : 3,')
    frag.write('	voxnormal.z > 0 ? 4 : 5')
    frag.write('	) * voxelgiResolution.x;')
    frag.write('vec3 direction_weights = abs(voxnormal);')

    frag.write('if (direction_weights.x > 0.0) {')
    frag.write('    vec3 basecol_direction = basecol * direction_weights.x;')
    frag.write('    uvw.x += face_offsets.x;')
    frag.write('    imageStore(voxels, ivec3(uvw), vec4(min(basecol_direction, vec3(1.0)), 1.0));')
    frag.write('}')

    frag.write('if (direction_weights.y > 0.0) {')
    frag.write('    vec3 basecol_direction = basecol * direction_weights.y;')
    frag.write('    uvw.x += face_offsets.y;')
    frag.write('    imageStore(voxels, ivec3(uvw), vec4(min(basecol_direction, vec3(1.0)), 1.0));')
    frag.write('}')

    frag.write('if (direction_weights.z > 0.0) {')
    frag.write('    vec3 basecol_direction = basecol * direction_weights.z;')
    frag.write('    uvw.x += face_offsets.z;')
    frag.write('    imageStore(voxels, ivec3(uvw), vec4(min(basecol_direction, vec3(1.0)), 1.0));')
    frag.write('}')

    return con_voxel


def make_ao(context_id):
    con_voxel = mat_state.data.add_context({ 'name': context_id, 'depth_write': False, 'compare_mode': 'always', 'cull_mode': 'none', 'color_writes_red': [False], 'color_writes_green': [False], 'color_writes_blue': [False], 'color_writes_alpha': [False], 'conservative_raster': False })
    wrd = bpy.data.worlds['Arm']
    rpdat = arm.utils.get_rp()

    vert = con_voxel.make_vert()
    frag = con_voxel.make_frag()
    geom = con_voxel.make_geom()
    tesc = None
    tese = None

    geom.ins = vert.outs
    frag.ins = geom.outs

    frag.add_include('compiled.inc')
    geom.add_include('compiled.inc')
    frag.add_include('std/math.glsl')
    frag.add_include('std/imageatomic.glsl')
    frag.write_header('#extension GL_ARB_shader_image_load_store : enable')
    frag.add_uniform('layout(r8) writeonly image3D voxels')

    vert.add_include('compiled.inc')
    vert.add_uniform('mat4 W', '_worldMatrix')
    vert.add_uniform('mat3 N', '_normalMatrix')
    vert.add_out('vec3 voxpositionGeom')
    vert.add_out('vec3 voxnormalGeom')

    vert.add_uniform('vec3 eye', '_cameraPosition')
    vert.add_uniform('int clipmapLevel', '_clipmapLevel')

    vert.write('vec3 P = vec3(W * vec4(pos.xyz, 1.0));')
    vert.write('float voxelSize = voxelgiVoxelSize * pow(2.0, clipmapLevel);')
    vert.write('float texelSize = 2.0 * voxelSize;')
    vert.write('vec3 clipmap_center = floor(eye / texelSize) * texelSize;')
    vert.write('voxpositionGeom = (P - clipmap_center) / voxelSize * 1.0 / voxelgiResolution.x;')
    vert.write('voxnormalGeom = normalize(N * vec3(nor.xy, pos.w));')

    geom.add_out('vec3 voxposition')
    geom.add_out('vec3 voxnormal')

    geom.write('vec3 p1 = voxpositionGeom[1] - voxpositionGeom[0];')
    geom.write('vec3 p2 = voxpositionGeom[2] - voxpositionGeom[0];')
    geom.write('vec3 p = abs(cross(p1, p2));')
    geom.write('for (uint i = 0; i < 3; ++i) {')
    geom.write('    voxposition = voxpositionGeom[i];')
    geom.write('    voxnormal = voxnormalGeom[i];')
    geom.write('    if (p.z > p.x && p.z > p.y) {')
    geom.write('        gl_Position = vec4(voxposition.x, voxposition.y, 0.0, 1.0);')
    geom.write('    }')
    geom.write('    else if (p.x > p.y && p.x > p.z) {')
    geom.write('        gl_Position = vec4(voxposition.y, voxposition.z, 0.0, 1.0);')
    geom.write('    }')
    geom.write('    else {')
    geom.write('        gl_Position = vec4(voxposition.x, voxposition.z, 0.0, 1.0);')
    geom.write('    }')
    geom.write('    EmitVertex();')
    geom.write('}')
    geom.write('EndPrimitive();')

    frag.write('if (abs(voxposition.z) > ' + rpdat.rp_voxelgi_resolution_z + ' || abs(voxposition.x) > 1 || abs(voxposition.y) > 1) return;')

    frag.add_uniform('int clipmapLevel', '_clipmapLevel')
    frag.write('vec3 uvw = (voxposition * 0.5 + 0.5) * voxelgiResolution.x;')
    frag.write('uvw.y += clipmapLevel * voxelgiResolution.x;')
    frag.write('vec3 face_offsets = vec3(')
    frag.write('	voxnormal.x > 0 ? 0 : 1,')
    frag.write('	voxnormal.y > 0 ? 2 : 3,')
    frag.write('	voxnormal.z > 0 ? 4 : 5')
    frag.write('	) * voxelgiResolution.x;')
    frag.write('vec3 direction_weights = abs(voxnormal);')

    frag.write('if (direction_weights.x > 0.0) {')
    frag.write('    float opac_direction = 1.0 * direction_weights.x;')
    frag.write('    uvw.x += face_offsets.x;')
    frag.write('    imageStore(voxels, ivec3(uvw), vec4(opac_direction));')
    frag.write('}')

    frag.write('if (direction_weights.y > 0.0) {')
    frag.write('    float opac_direction = 1.0 * direction_weights.y;')
    frag.write('    uvw.x += face_offsets.y;')
    frag.write('    imageStore(voxels, ivec3(uvw), vec4(opac_direction));')
    frag.write('}')

    frag.write('if (direction_weights.z > 0.0) {')
    frag.write('    float opac_direction = 1.0 * direction_weights.z;')
    frag.write('    uvw.x += face_offsets.z;')
    frag.write('    imageStore(voxels, ivec3(uvw), vec4(opac_direction));')
    frag.write('}')

    return con_voxel
