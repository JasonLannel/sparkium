#include "sparks/renderer/path_tracer.h"

#include "sparks/util/util.h"
#include "sparks/assets/pdf.h"
#include <glm/ext/scalar_constants.hpp>

#include <time.h>
#include <vector>

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
  struct Note {
    Ray ray;
    glm::vec3 throughput;
    int skip_id;
    int depth;
  };
  std::mt19937 rd(sample ^ x ^ y ^ std::time(0));
  std::uniform_real_distribution<float> RR_Validate(0.0f, 1.0f);
  Pdf *Light = scene_->GetLightPdf();
  std::vector<Note> q;
  q.push_back(Note{ray, glm::vec3(1), -1, render_settings_->num_bounces});
  glm::vec3 throughput;
  glm::vec3 radiance(0.0f);
  glm::vec3 origin, direction, normal, albedo;
  int skip_id;
  int depth;
  HitRecord hit_record;
  int idx = 0;
  while (idx < q.size()) {
    Note thead = q[idx];
    ++idx;
    ray = thead.ray;
    throughput = thead.throughput;
    skip_id = thead.skip_id;
    depth = thead.depth;
    if (depth < 0)
      continue;
    const float RR = 0.9, INV_RR = 1.0 / RR;
    if (scene_->TraceRay(ray.origin(), ray.direction(), 1e-3f, 1e4f,
                         &hit_record) > 0.0f) {
      auto &material =
          scene_->GetEntity(hit_record.hit_entity_id).GetMaterial();
      radiance += throughput * material.emission * material.emission_strength;
      if (material.material_type == MATERIAL_TYPE_EMISSION)
        continue;
      if (RR_Validate(rd) > RR) 
        continue;
      origin = hit_record.position;
      q.push_back(Note{Ray(origin, ray.direction(), ray.time()),
                  throughput * INV_RR * (1 - material.alpha),
                  hit_record.hit_entity_id, depth - 1});
      albedo = material.albedo_color;
      if (material.albedo_texture_id >= 0) {
        albedo *=
            glm::vec3{scene_->GetTextures()[material.albedo_texture_id].Sample(
                hit_record.tex_coord)};
      }
      normal = hit_record.normal;
      if (material.material_type == MATERIAL_TYPE_LAMBERTIAN) {
        if (glm::dot(normal, ray.direction()) > 0)
          normal = -normal;
          CosineHemispherePdf Dialectric(normal);
          Pdf *Gen;
          if (Light != nullptr) {
            Gen = new MixturePdf(&Dialectric, Light, 0.5f);
          } else {
            Gen = &Dialectric;
          }
          direction = Gen->Generate(origin, rd);
          float pdf = Gen->Value(origin, direction);
          float scatter = std::max(0.f, glm::dot(normal, direction) * INV_PI);
          throughput *= albedo * scatter / pdf;
          if (Light != nullptr)
            delete Gen;
        } else if (material.material_type == MATERIAL_TYPE_SPECULAR) {
          if (glm::dot(normal, ray.direction()) > 0)
            normal = -normal;
          throughput *= albedo;
          direction = glm::reflect(ray.direction(), normal);
        } else if (material.material_type == MATERIAL_TYPE_PRINCIPLED) {
            // 别忘记补上折射, 还没实现，这里是有问题的，只是用在 BRDF 时对
          if (glm::dot(normal, ray.direction()) > 0)
          normal = -normal;
          CosineHemispherePdf Dialectric(normal);
          Pdf *Gen;
          if (Light != nullptr) {
          Gen = new MixturePdf(&Dialectric, Light, 0.5f);
          } else {
            Gen = &Dialectric;
          }
          direction = Gen->Generate(origin, rd);
          Onb onb(normal);
          throughput *= material.DisneyPrincipled(
              normal, -ray.direction(), direction, onb.u(), onb.v(), albedo);
          if (Light != nullptr)
          delete Gen;
        }
        q.push_back(Note{Ray(origin, direction, ray.time()), material.alpha * INV_RR * throughput, -1, depth - 1});
    } else {
      radiance += throughput * glm::vec3{scene_->SampleEnvmap(ray.direction())};
    }
  }
  if (Light != nullptr)
    delete Light;
  return glm::min(radiance, glm::vec3(1));
}
}  // namespace sparks
