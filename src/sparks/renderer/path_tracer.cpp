#include "sparks/renderer/path_tracer.h"

#include "sparks/util/util.h"
#include "sparks/assets/pdf.h"
#include <glm/ext/scalar_constants.hpp>

#include <time.h>

namespace sparks {
PathTracer::PathTracer(const RendererSettings *render_settings,
                       const Scene *scene) {
  render_settings_ = render_settings;
  scene_ = scene;
}

glm::vec3 PathTracer::SampleRay(Ray ray,
                                int x,
                                int y,
                                int sample) const {
  glm::vec3 throughput{1.0f};
  glm::vec3 radiance{0.0f};
  glm::vec3 origin = ray.origin();
  glm::vec3 direction = ray.direction();
  double time = ray.time();
  HitRecord hit_record;
  const int max_bounce = render_settings_->num_bounces;
  const float RR = 0.9, INV_RR = 1.0 / RR;
  std::mt19937 rd(sample ^ x ^ y ^ std::time(0));
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  Pdf *emisgen = scene_->GetLightPdf();
  CosineHemispherePdf cosgen(direction);
  Pdf *gen;
  if (emisgen != nullptr) {
    gen = new MixturePdf(emisgen, &cosgen, 0.5f);
  } else {
    gen = &cosgen;
  }
  for (int bounce = 0; bounce < max_bounce; bounce++) {
    if (scene_->TraceRay(origin, direction, 1e-3f, 1e4f, &hit_record) > 0.0f) {
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
        origin = hit_record.position;
        glm::vec3 normal = hit_record.normal;
        if (dot(normal, direction) > 0)
          normal = -normal;
        if (material.material_type == MATERIAL_TYPE_LAMBERTIAN){
          cosgen = CosineHemispherePdf(normal);
          direction = gen->Generate(origin, rd);
          float pdf = gen->Value(origin, direction);
          float scatter = std::max(0.f, glm::dot(normal, direction) * INV_PI);
          throughput *= albedo * scatter / pdf;
        } else if (material.material_type == MATERIAL_TYPE_SPECULAR) {
          throughput *= albedo;
          direction = glm::reflect(direction, normal);
        } else if (material.material_type == MATERIAL_TYPE_PRINCIPLED) {
          cosgen = CosineHemispherePdf(normal);
          glm::vec3 dir_out = gen->Generate(origin, rd);
          throughput *= material.CookTorrance(normal, direction, dir_out);
          direction = dir_out;
        }
        if (dist(rd) > RR) {
          break;
        }
        albedo *= INV_RR;
      }
    } else {
      radiance += throughput * glm::vec3{scene_->SampleEnvmap(direction)};
      break;
    }
  }
  if (emisgen != nullptr)
    delete emisgen;
  return glm::min(radiance, glm::vec3(1));
}
}  // namespace sparks
