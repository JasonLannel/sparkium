#pragma once
#include "grassland/util/util.h"

namespace sparks {
constexpr float PI = 3.14159265358979323f;
constexpr float INV_PI = 0.31830988618379067f;
constexpr float INFTY = std::numeric_limits<float>::infinity();
constexpr float EPSILON = std::numeric_limits<float>::epsilon() * 0.5;

std::string PathToFilename(const std::string &file_path);
}  // namespace sparks
