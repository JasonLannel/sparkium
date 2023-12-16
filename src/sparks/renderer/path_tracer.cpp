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
  struct Note {
    Ray ray;
    glm::vec3 throughput;
    int skip_id;
    int depth;
  };
  std::mt19937 rd(sample ^ x ^ y ^ std::time(0));
  std::uniform_real_distribution<float> RandomProb(0.0f, 1.0f);
  Pdf *Light = scene_->GetLightPdf();
  glm::vec3 throughput(1.0f);
  glm::vec3 radiance(0.0f);
  glm::vec3 origin, direction, normal, albedo;
  HitRecord hit_record;
  for (register int i = 0, max_depth = render_settings_->num_bounces; i < max_depth;
       ++i) {
    const float RR = 0.9, INV_RR = 1.0 / RR;
    if (scene_->TraceRay(ray.origin(), ray.direction(), 1e-3f, 1e4f,
                         &hit_record) > 0.0f) {
      auto &material =
          scene_->GetEntity(hit_record.hit_entity_id).GetMaterial();
      radiance += throughput * material.emission * material.emission_strength;
      if (material.material_type == MATERIAL_TYPE_EMISSION)
        break;
      if (RandomProb(rd) > RR)
        break;
      origin = hit_record.position;
      if (RandomProb(rd) > material.alpha) {
        //Alpha Shadow
        ray = Ray(hit_record.position, ray.direction(), ray.time());
        continue;
      }
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
      } else if (material.material_type == MATERIAL_TYPE_TRANSMISSIVE) {
          // 折射材料处理
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
      ray = Ray(origin, direction, ray.time());
      throughput *= INV_RR;
    } else {
      radiance += throughput * glm::vec3{scene_->SampleEnvmap(ray.direction())};
      break;
    }
  }
  if (Light != nullptr)
    delete Light;
  return glm::min(radiance, glm::vec3(1));
}
}  // namespace sparks
