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
  const float RR = 0.9, INV_RR = 1.0 / RR;
  for (register int i = 0, max_depth = render_settings_->num_bounces; i < max_depth;
       ++i) {
    if (scene_->TraceRay(ray, 1e-3f, 1e4f,
                         &hit_record) > 0.0f) {
      auto &material =
          scene_->GetEntity(hit_record.hit_entity_id).GetMaterial(hit_record.material_id);
      radiance += throughput * material.emission * material.emission_strength;
      if (material.material_type == MATERIAL_TYPE_EMISSION)
        break;
      albedo = material.albedo_color;
      float alpha = material.alpha;
      if (material.albedo_texture_id >= 0) {
        auto tex_sample =
            scene_->GetTextures()[material.albedo_texture_id].Sample(
                hit_record.tex_coord);
        albedo *= glm::vec3(tex_sample);
        alpha *= tex_sample.w;
      }
      if (hit_record.material_id > 0 && RandomProb(rd) > alpha) {
        //Alpha Shadow
        ray = Ray(hit_record.position, ray.direction(), ray.time());
        --i;
        continue;
      }
      if (RandomProb(rd) > RR)
        break;
      origin = hit_record.position;
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
      glm::normalize(normal);
      glm::normalize(tangent);
      if (glm::dot(normal, ray.direction()) > 0)
        normal = -normal;
      if (material.material_type == MATERIAL_TYPE_LAMBERTIAN) {
          CosineHemispherePdf Dialectric(normal, tangent);
          Pdf *Gen = new MixturePdf(&Dialectric, Light, 1.f);
          direction = Gen->Generate(origin, ray.time(), rd);
          ray = Ray(origin, direction, ray.time());
          float pdf = Gen->Value(ray);
          float scatter = std::max(0.f, glm::dot(normal, direction) * INV_PI);
          float reflectance = material.reflectance;
          delete Gen;
          if (pdf < 1e-5)
            break;
          throughput *= albedo * reflectance * scatter / pdf;
          if (glm::dot(normal, direction) < 0)
            break;
      } else if (material.material_type == MATERIAL_TYPE_SPECULAR) {
          direction = glm::reflect(ray.direction(), normal);
          // Fuzz
          UniformHemispherePdf Fuzz(direction);
          glm::vec3 fuzzDirection = Fuzz.Generate(origin, ray.time(), rd);
          direction = glm::normalize(direction + fuzzDirection * material.fuzz);
          ray = Ray(origin, direction, ray.time());
          throughput *= albedo;
          if (glm::dot(normal, direction) < 0)
            break;
      } else if (material.material_type == MATERIAL_TYPE_TRANSMISSIVE) {
          // Assume all lights have same index of refraction
          float refract_ratio = hit_record.front_face ? (1.0 / material.IOR)
                                                      : material.IOR;
          if (material.thin) {
            float reflect_ratio = 1 - refract_ratio;
            reflect_ratio = (1 - reflect_ratio) * (1 - reflect_ratio) *
                            reflect_ratio / (1 - reflect_ratio * reflect_ratio);
            if (reflect_ratio > RandomProb(rd)) {
              direction = glm::reflect(ray.direction(), normal);
            } else {
              direction = ray.direction();
            }
          } else {
              float cos_theta = fmin(glm::dot(-ray.direction(), normal), 1.0f);
              direction = glm::refract(ray.direction(), normal, refract_ratio);
              if (direction.length() < 1e-4f ||
                  material.FresnelSchlick(refract_ratio, cos_theta) >
                      RandomProb(rd)) {
                direction = glm::reflect(ray.direction(), normal);
              }
          }
        glm::normalize(direction);
        ray = Ray(origin, direction, ray.time());
        throughput *= albedo;
      } //else if (material.material_type == MATERIAL_TYPE_PRINCIPLED) {
          // // 别忘记补上折射, 还没实现，这里是有问题的
          //CosineHemispherePdf Dialectric(normal, tangent);
          //Pdf *Gen = new MixturePdf(&Dialectric, Light, 0.5f);
          //direction = Gen->Generate(origin, ray.time(), rd);
          //Onb onb(normal);
          //throughput *= material.DisneyPrincipled(
              //normal, -ray.direction(), direction, onb.u(), onb.v(), albedo);
          //ray = Ray(origin, direction, ray.time());
          //float pdf = Gen->Value(ray);
          //throughput /= pdf;
          //delete Gen;
      //}
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
