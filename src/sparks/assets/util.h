﻿#pragma once
#include "glm/glm.hpp"
#include "string"
#include "tinyxml2.h"
#include <sparks/util/util.h>

namespace sparks {
glm::vec3 DecomposeRotation(glm::mat3 R);

glm::mat4 ComposeRotation(glm::vec3 pitch_yaw_roll);

glm::vec2 StringToVec2(const std::string &s);

glm::vec3 StringToVec3(const std::string &s);

glm::vec4 StringToVec4(const std::string &s);

glm::mat4 XmlTransformMatrix(tinyxml2::XMLElement *transform_element);

glm::mat4 XmlComposeTransformMatrix(tinyxml2::XMLElement *object_element);

template <typename T>
inline T sigmoid(T x) {
  return 1.0 / (1 + exp(-(2 * x - 1))) * 2.164 - 0.582;
}

template <typename T>
inline T interpolate(const T &val1, const T &val2, const float &ratio) {
  // NOTICE: This interpolate function is different from the Lerp function in pbrt
  // Lerp is ratio, val1, val2! 
  return T((val2 - val1) * ratio + val1);
}

template <typename T>
inline T pow5(const T &val) {
  return val * val * val * val * val;
};

template <typename T>
inline T square(const T &val) {
  return val * val;
};

template <typename T, typename U, typename V>
inline T clamp(T val, U low, V high) {
  if (val < low)
    return low;
  else if (val > high)
    return high;
  else
    return val;
}

inline bool floatEq(const float &lhs, const float &rhs) {
  return std::abs(lhs - rhs) <= EPSILON * 1e2;
}

// NOTICE: below all assume vectors are in a y-up tangent space. 

inline float CosTheta(const glm::vec3 &w) {
  return w.y;
}

inline float calAbsCosTheta(const glm::vec3 &vec) {
  return std::abs(vec.y);
}

inline float Cos2Theta(const glm::vec3 &w) {
  return w.y * w.y;
}

inline float Sin2Theta(const glm::vec3 &w) {
  return std::max(0.0f, 1.0f - Cos2Theta(w));
}

inline float SinTheta(const glm::vec3 &w) {
  return sqrt(Sin2Theta(w));
}

inline float CosPhi(const glm::vec3 &w) {
  float sinTheta = SinTheta(w);
  return (sinTheta == 0) ? 1.0f : clamp(w.x / sinTheta, -1.0f, 1.0f);
}

inline float SinPhi(const glm::vec3 &w) {
  float sinTheta = SinTheta(w);
  return (sinTheta == 0) ? 1.0f : clamp(w.z / sinTheta, -1.0f, 1.0f);
}

inline float Cos2Phi(const glm::vec3 &w) {
  return CosPhi(w) * CosPhi(w);
}

inline float Sin2Phi(const glm::vec3 &w) {
  return SinPhi(w) * SinPhi(w);
}

inline float calAbsTanTheta(const glm::vec3 &vec) {
  if (vec.y == 0.0f)
    return INFTY;
  return std::abs(SinTheta(vec) / CosTheta(vec));
}

inline float calPhi(const glm::vec3 &vec) {
  float phi0_pi = std::acos(vec.x / std::sqrt(vec.x * vec.x + vec.y * vec.y));
  return vec.y >= 0 ? phi0_pi : PI * 2 - phi0_pi;
}

inline float calTheta(const glm::vec3 &vec) {
  return std::asin(vec.z / glm::length(vec));
}

inline glm::vec3 SphereToXYZ(const float &phi, const float &theta) {
  using namespace std;
  return glm::vec3(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
}

inline float FrDielectric(float cosThetaI, float etaI, float etaT) {
  cosThetaI = clamp(cosThetaI, -1, 1);
  // Potentially swap indices of refraction
  bool entering = cosThetaI > 0.f;
  if (!entering) {
    std::swap(etaI, etaT);
    cosThetaI = std::abs(cosThetaI);
  }

  // Compute _cosThetaT_ using Snell's law
  float sinThetaI = sqrt(std::max(0.0f, 1.0f - cosThetaI * cosThetaI));
  float sinThetaT = etaI / etaT * sinThetaI;

  // Handle total internal reflection
  if (sinThetaT >= 1)
    return 1;
  float cosThetaT = sqrt(std::max(0.0f, 1.0f - sinThetaT * sinThetaT));
  float Rparl = ((etaT * cosThetaI) - (etaI * cosThetaT)) /
                ((etaT * cosThetaI) + (etaI * cosThetaT));
  float Rperp = ((etaI * cosThetaI) - (etaT * cosThetaT)) /
                ((etaI * cosThetaI) + (etaT * cosThetaT));
  return (Rparl * Rparl + Rperp * Rperp) / 2;
}

inline bool Transmit(glm::vec3 wm, glm::vec3 wi, float n, glm::vec3 &wo) {
  float c = glm::dot(wi, wm);
  if (c < 0.0f) {
    c = -c;
    wm = -wm;
  }
  float root = 1.0f - n * n * (1.0f - c * c);
  if (root <= 0)
    return false;
  wo = (n * c - sqrt(root)) * wm - n * wi;
  return true;
}

}  // namespace sparks
