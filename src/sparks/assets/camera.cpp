#include "sparks/assets/camera.h"

#include "glm/gtc/matrix_transform.hpp"
#include "imgui.h"
#include "sparks/util/util.h"

namespace sparks {

glm::mat4 Camera::GetProjectionMatrix(float aspect,
                                      float t_min,
                                      float t_max) const {
  return glm::scale(glm::mat4{1.0f}, glm::vec3{1.0f, -1.0f, 1.0f}) *
         glm::perspectiveZO(glm::radians(fov_), aspect, t_min, t_max);
}

bool Camera::ImGuiItems() {
  bool value_changed = false;
  value_changed |= ImGui::SliderFloat("FOV", &fov_, 1.0f, 160.0f, "%.0f", 0);
  value_changed |=
      ImGui::SliderFloat("Aperture", &aperture_, 0.0f, 1.0f, "%.2f");
  value_changed |=
      ImGui::SliderFloat("Focal Length", &focal_length_, 0.1f, 10000.0f, "%.2f",
                         ImGuiSliderFlags_Logarithmic);
  value_changed |= ImGui::SliderFloat("Clamp", &clamp_, 1.0f, 1000000.0f,
                                      "%.2f", ImGuiSliderFlags_Logarithmic);
  ImGui::SliderFloat("Gamma", &gamma_, 0.1f, 10.0f);
  return value_changed;
}

void Camera::UpdateFov(float delta) {
  fov_ += delta;
  fov_ = glm::clamp(fov_, 1.0f, 160.0f);
}

void Camera::GenerateRay(float aspect,
                         glm::vec2 range_low,
                         glm::vec2 range_high,
                         Ray &ray, 
                         float rand_u,
                         float rand_v,
                         float rand_w,
                         float rand_r) const {
  auto pos = (range_high - range_low) * glm::vec2{rand_u, rand_v} + range_low;
  pos = pos * 2.0f - 1.0f;
  pos.y *= -1.0f;
  glm::vec3 origin = glm::vec3{0.0f};
  auto tan_fov = std::tan(glm::radians(fov_ * 0.5f));
  float theta = 2.0f * PI * rand_w;
  float sin_theta = std::sin(theta);
  float cos_theta = std::cos(theta);
  origin =
      glm::vec3{glm::vec2{sin_theta, cos_theta} * rand_r * aperture_, 0.0f};
  glm::vec3 direction = glm::normalize(
      glm::vec3{tan_fov * aspect * pos.x, tan_fov * pos.y, -1.0f} *
          focal_length_ -
      origin);
  double time = time_0_ + rand_u * (time_1_ - time_0_);
  ray = Ray{origin, direction, time};
}

Camera::Camera(float fov, float aperture, float focal_length, double time_0, double time_1)
    : fov_(fov), aperture_(aperture), focal_length_(focal_length), time_0_(time_0), time_1_(time_1) {
}
}  // namespace sparks
