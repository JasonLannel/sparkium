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

glm::vec3 PathTracer::SampleRay(glm::vec3 origin,
                                glm::vec3 direction,
                                int x,
                                int y,
                                int sample) const {
  glm::vec3 throughput{1.0f};
  glm::vec3 radiance{0.0f};
  HitRecord hit_record;
  const int max_bounce = render_settings_->num_bounces;
  std::mt19937 rd(sample ^ x ^ y ^ std::time(0));
  Pdf *emisgen = scene_->GetLightPdf();
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
          Pdf *gen;
          Pdf *cosgen = new CosineHemispherePdf(normal);
          if (emisgen != nullptr){
            gen = new MixturePdf(emisgen, cosgen, 0.5f);
          } else {
            gen = cosgen;
          }
          direction = gen->Generate(origin, rd);
          float pdf = gen->Value(origin, direction);
          float scatter = std::max(0.f, glm::dot(normal, direction) * INV_PI);
          throughput *= albedo * scatter / pdf;
          delete gen;
          if (emisgen != nullptr) {
            delete cosgen;
          }
          /*
          if (emisgen != nullptr) {
            glm::vec3 dir = emisgen->Generate(origin, rd);
            HitRecord emisRec;
            if (scene_->TraceRay(origin, dir, 1e-2, 1e4f, &emisRec) > 0.0f) {
              float pdf = emisgen->Value(origin, dir);
              float scatter =
                  std::max(0.f, glm::dot(normal, dir) * INV_PI);
              auto &material =
                  scene_->GetEntity(emisRec.hit_entity_id).GetMaterial();
              if (material.material_type == MATERIAL_TYPE_EMISSION)
                  radiance += 0.5f * throughput * albedo * scatter * material.emission * material.emission_strength / pdf;
            }
          }
          Pdf *gen = new CosineHemispherePdf(normal);
          direction = gen->Generate(origin, rd);
          float pdf = gen->Value(origin, direction);
          float scatter = std::max(0.f, glm::dot(normal, direction) * INV_PI);
          throughput *= 0.5f * albedo * scatter / pdf;
          delete gen;
          */
        } else if (material.material_type == MATERIAL_TYPE_SPECULAR) {
          throughput *= albedo;
          direction = glm::reflect(direction, normal);
        }
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
