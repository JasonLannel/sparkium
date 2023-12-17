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
  std::random_device rd_seed;
  std::mt19937 rd(x ^ y ^ sample ^ rd_seed());
  std::uniform_real_distribution<float> RandomProb(0.0f, 1.0f);
  Pdf *Light = scene_->GetLightPdf();
  glm::vec3 throughput(1.0f);
  glm::vec3 radiance(0.0f);
  glm::vec3 origin, direction, normal, tangent, albedo;
  HitRecord hit_record;
  for (register int i = 0, max_depth = render_settings_->num_bounces; i < max_depth;
       ++i) {
    const float RR = 0.9, INV_RR = 1.0 / RR;
    if (scene_->TraceRay(ray, 1e-3f, 1e4f,
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
      tangent = hit_record.tangent;
      if (material.use_normal_texture) {
        Onb onb(normal, tangent);
        glm::vec3 normalFromTex =
            glm::vec3{scene_->GetTextures()[material.normal_texture_id].Sample(
                hit_record.tex_coord)};
        normalFromTex = normalFromTex * 2.0f - glm::vec3(1.0f);
        normalFromTex[0] *= material.bumpScale;
        normalFromTex[1] *= material.bumpScale;
        normalFromTex[2] = 0;
        normalFromTex[2] = sqrt(1.f - glm::dot(normalFromTex, normalFromTex));
        normal = onb.local(normalFromTex);
      }
      if (material.material_type == MATERIAL_TYPE_LAMBERTIAN) {
        if (glm::dot(normal, ray.direction()) > 0)
          normal = -normal;
          CosineHemispherePdf Dialectric(normal, tangent);
          Pdf *Gen;
          if (Light != nullptr) {
            Gen = new MixturePdf(&Dialectric, Light, 0.5f);
          } else {
            Gen = &Dialectric;
          }
          direction = Gen->Generate(origin, ray.time(), rd);
          ray = Ray(origin, direction, ray.time());
          float pdf = Gen->Value(ray);
          float scatter = std::max(0.f, glm::dot(normal, direction) * INV_PI);
          float reflectance = material.reflectance;
          throughput *= albedo * reflectance * scatter / pdf;
          if (Light != nullptr)
            delete Gen;
      } else if (material.material_type == MATERIAL_TYPE_SPECULAR) {
          if (glm::dot(normal, ray.direction()) > 0)
            normal = -normal;
          direction = glm::reflect(ray.direction(), normal);
          // Fuzz
          UniformHemispherePdf Fuzz(direction);
          glm::vec3 fuzzDirection = Fuzz.Generate(origin, ray.time(), rd);
          direction = glm::normalize(direction + fuzzDirection * material.fuzz);
          if (glm::dot(normal, direction) < 0)
            break;  // Absorb Energy.
          ray = Ray(origin, direction, ray.time());
          throughput *= albedo;
      } else if (material.material_type == MATERIAL_TYPE_TRANSMISSIVE) {
          // 折射材料处理
      } else if (material.material_type == MATERIAL_TYPE_PRINCIPLED) {
          // 别忘记补上折射, 还没实现，这里是有问题的，只是用在 BRDF 时对
          if (glm::dot(normal, ray.direction()) > 0)
          normal = -normal;
          CosineHemispherePdf Dialectric(normal, tangent);
          Pdf *Gen;
          if (Light != nullptr) {
            Gen = new MixturePdf(&Dialectric, Light, 0.5f);
          } else {
            Gen = &Dialectric;
          }
          direction = Gen->Generate(origin, ray.time(), rd);
          Onb onb(normal);
          throughput *= material.DisneyPrincipled(
              normal, -ray.direction(), direction, onb.u(), onb.v(), albedo);
          ray = Ray(origin, direction, ray.time());
          float pdf = Gen->Value(ray);
          throughput /= pdf;
          if (Light != nullptr)
            delete Gen;
      }
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
