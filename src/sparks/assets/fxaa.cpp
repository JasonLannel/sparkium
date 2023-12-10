#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <math.h>

#include "sparks/assets/fxaa.h"

#define QUALITY(q) ((q) < 5 ? 1.0 : ((q) > 5 ? ((q) < 10 ? 2.0 : ((q) < 11 ? 4.0 : 8.0)) : 1.5))

namespace sparks {

FXAA::FXAA(uint32_t width, uint32_t height) {
    width_ = width;
    height_ = height;
    luminance_ = std::vector<float>(width*height, 0.0f);
}

uint32_t FXAA::idx(uint32_t i, uint32_t j) {
    return i + j* width_;
}

inline float saturate(float value) {
    if (value < 0.0f) {
        return 0.0f;
    } else if (value > 1.0f) {
        return 1.0f;
    } else {
        return value;
    }
}

inline float clamp(float x, float lowerlimit, float upperlimit) {
    if (x < lowerlimit) {
        return lowerlimit;
    } else if (x > upperlimit) {
        return upperlimit;
    } else {
        return x;
    }
}

inline float smoothstep(float edge0, float edge1, float x)
{
    // Scale, bias and saturate x to 0..1 range
    x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
    // Evaluate polynomial
    return x*x*(3 - 2 * x);
}

inline float rgb2luma(glm::vec4 rgb) {
    return 0.213f * rgb[0] + 0.715f * rgb[1] + 0.072f * rgb[2];
}

void FXAA::image_to_texture(std::vector<glm::vec4> & image) {
    texture_ = Texture(width_, height_, image.data(), SAMPLE_TYPE_LINEAR);
}

void FXAA::calcLuma(std::vector<glm::vec4> & image) {
    // calculate the luminance of each pixel, write it in luminance_
    //  L = 0.213 * R + 0.715 * G + 0.072 * B
    for (uint32_t i=0; i<width_; i++) {
        for (uint32_t j=0; j<height_; j++) {
            uint32_t ind = idx(i,j);
            luminance_[ind] = rgb2luma(image[ind]); 
        }
    }
}

glm::vec4 FXAA::ModifyPixel(uint32_t x, uint32_t y) {
    float x1 = (float(x)+ 0.5f) / float(width_);
    float y1 = (float(y)+ 0.5f) / float(height_);
    glm::vec2 texCoord = glm::vec2(x1, y1);
    uint32_t ind = idx(x,y);
    glm::vec4 color_M = texture_.Sample(texCoord);
    float M = luminance_[ind];
    float N = 0.0f;
    float E = 0.0f;
    float S = 0.0f;
    float W = 0.0f;
    float NE = 0.0f;
    float NW = 0.0f;
    float SE = 0.0f;
    float SW = 0.0f;
    if (x>0) {
        W = luminance_[idx(x-1,y)];
    }
    if (x<width_-1) {
        E = luminance_[idx(x+1,y)];
    }
    if (y>0) {
        N = luminance_[idx(x,y-1)];
    }
    if (y<height_-1) {
        S = luminance_[idx(x,y+1)];
    }
    if ((x>0) && (y>0)) {
        NW = luminance_[idx(x-1,y-1)];
    }
    if ((x>0) && (y<height_-1)) {
        SW = luminance_[idx(x-1,y+1)];
    }
    if ((x<width_-1) && (y>0)) {
        NE = luminance_[idx(x+1,y-1)];
    }
    if ((x<width_-1) && (y<height_-1)) {
        SE = luminance_[idx(x+1,y+1)];
    }

    // first calculate the contrast of the pixel with the N, E, S, W neighbors
    float maxLuminance = std::max(std::max(std::max(std::abs(N), std::abs(E)), std::abs(S)), std::abs(W));
    float minLuminance = std::min(std::min(std::min(std::abs(N), std::abs(E)), std::abs(S)), std::abs(W));
    float contrast = maxLuminance - minLuminance;
    if (contrast < std::max (_MinThreshold, _Threshold * maxLuminance)) {
        return color_M;
    }

    // decide the coefficient when we choose to mix the color of the pixel with its neighbors
    float filter = 0.0f;
    filter = 2.0f * (N + E + S + W) + NE + NW + SE + SW;
    filter /= 12.0f;
    filter = std::abs(filter - M);
    filter = saturate(filter / contrast);
    float PixelBlend = smoothstep(0, 1, filter);
    PixelBlend = PixelBlend * PixelBlend;

    // decide the edge direction, whether it is horizontal or vertical
    float Vertical = abs(N + S - 2 * M) * 2+ abs(NE + SE - 2 * E) + abs(NW + SW - 2 * W);
    float Horizontal = abs(E + W - 2 * M) * 2 + abs(NE + NW - 2 * N) + abs(SE + SW - 2 * S);
    bool isHorizontal = Vertical > Horizontal;
    glm::vec2 PixelStep = isHorizontal ? glm::vec2(0.0f, 1.0f / float(height_)) : glm::vec2(1.0f / float(width_), 0.0f);

    // calculate gradient in the edges
    float luma1 = isHorizontal ? S : E;
    float luma2 = isHorizontal ? N : W;
    float gradient1 = abs(luma1 - M);
    float gradient2 = abs(luma2 - M);
    bool isSteepest = gradient1 >= gradient2;
    float gradientScaled = 0.25f * std::max(gradient1, gradient2);

    float lumaLocalAverage = 0.0;
	if(isSteepest){
		// Switch the direction
		PixelStep = -PixelStep;
        lumaLocalAverage = 0.5 * (luma1 + M);
	} else {
		lumaLocalAverage = 0.5 * (luma2 + M);
	}

    glm::vec2 currentUv = texCoord;
    currentUv = currentUv + PixelStep * 0.5f;

    // compute offset and explore on each side of the edge.
    glm::vec2 offset = isHorizontal ? glm::vec2(1.0f / float(width_), 0.0f) : glm::vec2(0.0f, 1.0f / float(height_));
    glm::vec2 uv1 = currentUv - offset * glm::vec2(QUALITY(0));
    glm::vec2 uv2 = currentUv + offset * glm::vec2(QUALITY(0));

    float lumaEnd1;
	float lumaEnd2;
	bool reached1 = false;
	bool reached2 = false;
	bool reachedBoth = false;

    for(int i = 1; i < 12; i++){
		// If needed, read luma, compute delta. 
		if(!reached1){
			lumaEnd1 = rgb2luma(texture_.Sample(uv1));
			lumaEnd1 = lumaEnd1 - lumaLocalAverage;
			reached1 = abs(lumaEnd1) >= gradientScaled;
		}
		if(!reached2){
			lumaEnd2 = rgb2luma(texture_.Sample(uv2));
			lumaEnd2 = lumaEnd2 - lumaLocalAverage;
			reached2 = abs(lumaEnd2) >= gradientScaled;
		}
		reachedBoth = reached1 && reached2;

        // If not reached, continue to explore.
		if(!reached1){
			uv1 -= offset * glm::vec2(QUALITY(i));
		}
		if(!reached2){
			uv2 += offset * glm::vec2(QUALITY(i));
		}
		if(reachedBoth){ break;}
    }

	float distance1 = isHorizontal ? (texCoord[0] - uv1[0]) : (texCoord[1] - uv1[1]);
	float distance2 = isHorizontal ? (uv2[0] - texCoord[0]) : (uv2[1] - texCoord[1]);
	
	// In which direction is the side of the edge closer ?
	bool isDirection1 = distance1 < distance2;
	float distanceFinal = std::min(distance1, distance2);
	
	// Is the luma at center smaller than the local average ?
	bool isLumaCenterSmaller = M < lumaLocalAverage;
	bool correctVariation1 = (lumaEnd1 < 0.0) != isLumaCenterSmaller;
	bool correctVariation2 = (lumaEnd2 < 0.0) != isLumaCenterSmaller;
	
	bool correctVariation = isDirection1 ? correctVariation1 : correctVariation2;
	float edgeLength = (distance1 + distance2);
	float pixelOffset = - distanceFinal / edgeLength + 0.5;
	float finalOffset = correctVariation ? pixelOffset : 0.0;
	
	// Pick the biggest of the two offsets.
	finalOffset = std::max(finalOffset, PixelBlend);
	
	// Compute the final UV coordinates.
	glm::vec2 finalUv = texCoord;
	finalUv += finalOffset * PixelStep;
	
	glm::vec3 finalColor = texture_.Sample(finalUv);
	return glm::vec4(finalColor, 1.0f);
}

std::vector<glm::vec4> FXAA::Apply(std::vector<glm::vec4> & image) {
    calcLuma(image);
    image_to_texture(image);
    std::vector<glm::vec4> result(width_*height_, glm::vec4(0.0f));
    for (uint32_t i=0; i<width_; i++) {
        for (uint32_t j=0; j<height_; j++) {
            result[idx(i,j)] = ModifyPixel(i,j);
        }
    }
    return result;
}

} // namespace sparks
