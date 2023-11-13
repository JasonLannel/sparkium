#include "sparks/renderer/path_tracer.h"

#include "sparks/util/util.h"
#include <glm/ext/scalar_constants.hpp>

namespace sparks {
PathTracer::PathTracer(const RendererSettings *render_settings,
                       const Scene *scene) {
  render_settings_ = render_settings;
  scene_ = scene;
}

glm::vec3 PathTracer::RandomUnitVector(std::mt19937& rd) {
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  float theta = dist(rd) * 2.0f * glm::pi<float>();
  float phi = dist(rd) * glm::pi<float>();
  return glm::vec3{std::cos(theta) * std::sin(phi), std::sin(theta) * std::sin(phi), std::cos(phi)};
}
inline const glm::vec3 CosineWeightedSampleOnHemisphere(
    std::mt19937 &rd) noexcept {
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  float u1 = dist(rd);
  float u2 = dist(rd);
  const double cos_theta = sqrt(1.0 - u1);
  const double sin_theta = sqrt(u1);
  const double phi = 2.0 * PI * u2;
  return {std::cos(phi) * sin_theta, std::sin(phi) * sin_theta, cos_theta};
}
glm::vec3 PathTracer::SampleRay(glm::vec3 origin,
                                glm::vec3 direction,
                                int x,
                                int y,
                                int sample) const {
  glm::vec3 throughput{1.0f};
  glm::vec3 radiance{0.0f};
  HitRecord hit_record;
  const int max_bounce = render_settings_->num_bounces;
  std::mt19937 rd(sample ^ x ^ y);

  for (int bounce = 0; bounce < max_bounce; bounce++) {
    auto t = scene_->TraceRay(origin, direction, 1e-3f, 1e4f, &hit_record);
    if (t > 0.0f) {
      auto &material =
          scene_->GetEntity(hit_record.hit_entity_id).GetMaterial();
      if (material.material_type == MATERIAL_TYPE_EMISSION) {
        radiance += throughput * material.emission * material.emission_strength;
        break;
      } else {
        glm::vec3 albedo = material.albedo_color;
        if (material.albedo_texture_id >= 0) {
          albedo *= glm::vec3{
              scene_->GetTextures()[material.albedo_texture_id].Sample(
                  hit_record.tex_coord)};
        }
        throughput *= albedo;
        origin = hit_record.position;
        if (material.material_type == MATERIAL_TYPE_LAMBERTIAN){
          glm::vec3 n = hit_record.normal;
          glm::vec3 w = (glm::dot(n, direction) < 0) ? n : -n;
          glm::vec3 u =
              glm::normalize(glm::cross(w, 
                  (std::abs(w.x) > 0.1 ? glm::vec3(0.0, 1.0, 0.0): glm::vec3(1.0, 0.0, 0.0))));
          glm::vec3 v = glm::cross(u, w);
          glm::vec3 sample_d =
              CosineWeightedSampleOnHemisphere(rd);
          direction =
              glm::normalize(sample_d.x * u + sample_d.y * v + sample_d.z * w);

          float pdf = 0.5f / PI;
          float scatter_pdf = glm::dot(n, direction) / PI;
        } else if (material.material_type == MATERIAL_TYPE_SPECULAR) {
          direction = glm::reflect(direction, hit_record.normal);
        }
      }
    } else {
      radiance += throughput * glm::vec3{scene_->SampleEnvmap(direction)};
      break;
    }
  }
  return radiance;
}
}  // namespace sparks
