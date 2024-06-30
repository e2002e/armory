#version 450

#include "compiled.inc"
#include "std/gbuffer.glsl"

uniform sampler2D gbufferD;
uniform sampler2D gbuffer0;
uniform sampler2D tex;
uniform sampler2D gbufferS;

uniform vec2 dir;
uniform vec2 cameraProj;

in vec2 texCoord;
out vec4 fragColor;

const float SSSS_FOVY = 45.0;
const float SSSS_STRENGTH = 4.0; // Strength parameter to enhance SSS effect

float gaussian(float x, float sigma) {
    return exp(-0.5 * (x * x) / (sigma * sigma)) / (sigma * sqrt(2.0 * 3.141592653589793));
}

vec4 SSSSBlur(vec3 subsurface_color) {
    const int SSSS_N_SAMPLES = 31; // Increased sample count
    float kernel[SSSS_N_SAMPLES];
    float sigma = 0.2;
    float totalWeight = 0.0;

    // Generate Gaussian weights
    for (int i = 0; i < SSSS_N_SAMPLES; i++) {
        float x = (float(i) - float(SSSS_N_SAMPLES / 2)) / float(SSSS_N_SAMPLES / 2);
        kernel[i] = gaussian(x, sigma);
        totalWeight += kernel[i];
    }

    // Normalize weights and increase the influence
    for (int i = 0; i < SSSS_N_SAMPLES; i++) {
        kernel[i] /= totalWeight;
        kernel[i] *= SSSS_STRENGTH; // Increase kernel influence
    }

    vec4 colorM = textureLod(tex, texCoord, 0.0);

    // Fetch linear depth of current pixel
    float depth = textureLod(gbufferD, texCoord, 0.0).r * 2.0 - 1.0;
    float depthM = cameraProj.y / (depth - cameraProj.x);

    // Calculate the sssWidth scale (1.0 for a unit plane sitting on the projection window)
    float distanceToProjectionWindow = 1.0 / tan(0.5 * radians(SSSS_FOVY));
    float scale = distanceToProjectionWindow / depthM;

    // Calculate the adaptive step size
    float adaptiveStepSize = sssWidth * scale;
    adaptiveStepSize /= (SSSS_N_SAMPLES - 1); // Spread the kernel over the samples

    // Accumulate the center sample:
    vec4 colorBlurred = colorM * kernel[SSSS_N_SAMPLES / 2] * vec4(subsurface_color, 1.0);

    // Accumulate the other samples with some random jittering
    for (int i = 0; i < SSSS_N_SAMPLES; i++) {
        if (i != SSSS_N_SAMPLES / 2) {
            float offsetScale = (float(i) - float(SSSS_N_SAMPLES / 2)) / float(SSSS_N_SAMPLES / 2);
            vec2 offset = texCoord + offsetScale * adaptiveStepSize * dir;
            vec4 color = textureLod(tex, offset, 0.0);

			// Optional: Follow surface depth to avoid bleeding
            float sampleDepth = textureLod(gbufferD, offset, 0.0).r * 2.0 - 1.0;
            float sampleDepthM = cameraProj.y / (sampleDepth - cameraProj.x);
            float depthDifference = abs(depthM - sampleDepthM);
            float weight = exp(-depthDifference); // Adjust exponential factor for better surface following
            color.rgb = mix(color.rgb * subsurface_color, colorM.rgb, weight);

            // Accumulate
            colorBlurred.rgb += kernel[i] * color.rgb;
        }
    }

    return colorBlurred;
}

void main() {
    // SSS only for masked objects
    float metallic;
    uint matid;
    unpackFloatInt16(textureLod(gbuffer0, texCoord, 0.0).a, metallic, matid);

    vec4 finalColor = textureLod(tex, texCoord, 0.0); // Default to background texture
    if (matid == 2) {
        vec3 subsurface_color = textureLod(gbufferS, texCoord, 0.0).rgb;
        vec4 sssColor = clamp(SSSSBlur(subsurface_color), 0.0, 1.0);
        finalColor = mix(finalColor, sssColor, sssColor.a); // Blend SSS color with background
    }
    fragColor = finalColor;
}
