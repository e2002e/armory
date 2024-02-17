uint convVec4ToRGBA8(vec4 val) {
    vec4 col = val * 255.0;
    return (uint(col.w) & 0x000000FF) << 24U
         | (uint(col.z) & 0x000000FF) << 16U
         | (uint(col.y) & 0x000000FF) << 8U
         | (uint(col.x) & 0x000000FF);
}

vec4 convRGBA8ToVec4(uint val) {
    uvec4 col = uvec4(
        (val & 0x000000FF),
        (val & 0x0000FF00) >> 8U,
        (val & 0x00FF0000) >> 16U,
        (val & 0xFF000000) >> 24U
    );
    return vec4(col) / 255.0;
}

