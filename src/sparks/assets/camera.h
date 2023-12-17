#pragma once
#include "glm/glm.hpp"
#include "sparks/assets/ray.h"

namespace sparks {
class Camera {
 public:
  Camera(float fov = 60.0f, float aperture = 0.0f, float focal_length = 3.0f, float time_0 = 0.0, float time_1 = 1.0, float focus_distance = 4.0f);
  [[nodiscard]] glm::mat4 GetProjectionMatrix(float aspect,
                                              float t_min,
                                              float t_max) const;
  void GenerateRay(float aspect,
                   glm::vec2 range_low,
                   glm::vec2 range_high,
                   Ray &ray, 
                   float rand_u = 0.0f,
                   float rand_v = 0.0f,
                   float rand_w = 0.0f,
                   float rand_r = 0.0f,
                   float rand_t = 0.0f) const;
  bool ImGuiItems();
  void UpdateFov(float delta);
  [[nodiscard]] float GetFov() const {
    return fov_;
  }
  [[nodiscard]] float GetAperture() const {
    return aperture_;
  }
  [[nodiscard]] float GetFocalLength() const {
    return focal_length_;
  }
  [[nodiscard]] float GetClamp() const {
    return clamp_;
  }
  [[nodiscard]] float GetGamma() const {
    return gamma_;
  }
  [[nodiscard]] float GetFrontDOF() const {
    return CoC_ * focus_distance_ * (focus_distance_-focal_length_) /
            (aperture_ * focal_length_ + CoC_ * (focus_distance_-focal_length_));
  }
  [[nodiscard]] float GetBackDOF() const {
    if (aperture_ * focal_length_ - CoC_ * (focus_distance_ - focal_length_) <
        1e-4)
      return CoC_ * focus_distance_ * (focus_distance_ - focal_length_) * 1e4;
    return CoC_ * focus_distance_ * (focus_distance_ - focal_length_) /
           (aperture_ * focal_length_ -
            CoC_ * (focus_distance_ - focal_length_));
  }

 private:
  float fov_{60.0f};
  float aperture_{0.0f};
  float focal_length_{3.0f};
  float clamp_{100.0f};
  float gamma_{2.2f};
  // motion blur parameters
  float time_0_{0.0};
  float time_1_{1.0};
  // field of view
  float focus_distance_{4.0f};
  float CoC_{0.2f};
};
}  // namespace sparks
