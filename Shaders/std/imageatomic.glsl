uint convVec4ToRGBA8(vec4 val) {
    val = clamp(val, 0.0, 1.0); // Ensure values are in the [0, 1] range
    return uint(val.w * 255.0) << 24U |
           uint(val.z * 255.0) << 16U |
           uint(val.y * 255.0) << 8U |
           uint(val.x * 255.0);
}

vec4 convRGBA8ToVec4(uint val) {
    return vec4(
        float((val & 0x000000FF) >> 0U) / 255.0,
        float((val & 0x0000FF00) >> 8U) / 255.0,
        float((val & 0x00FF0000) >> 16U) / 255.0,
        float((val & 0xFF000000) >> 24U) / 255.0
    );
}
