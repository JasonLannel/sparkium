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
  LightPdf Light(scene_);
  glm::vec3 throughput(1.0f);
  glm::vec3 radiance(0.0f);
  glm::vec3 origin, direction, normal, tangent, albedo;
  HitRecord hit_record;
  Ray next_ray(ray), light_ray(ray);
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
      if (material.albedo_texture_id >= 0 &&
          material.material_type != MATERIAL_TYPE_MEDIUM) {
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
      if (material.use_normal_texture &&
          material.material_type != MATERIAL_TYPE_MEDIUM) {
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
        normal = glm::normalize(normal);
        tangent = glm::normalize(tangent - glm::dot(tangent, normal) * normal);
      }
      if (glm::dot(normal, ray.direction()) > 0)
        normal = -normal;
      if (material.material_type == MATERIAL_TYPE_LAMBERTIAN) {
          CosineHemispherePdf Lambert(normal, tangent);
          // Sample Light First.
          float Lamb_pdf, Light_pdf, w;
          direction = Light.Generate(origin, ray.time(), rd);
          light_ray = Ray(origin, direction, ray.time());
          Lamb_pdf = Lambert.Value(light_ray);
          Light_pdf = Light.Value(light_ray);
          w = square(Light_pdf) / (square(Lamb_pdf) + square(Light_pdf));
          HitRecord light_hit;
          if (scene_->TraceRay(light_ray, 1e-3f, 1e4f, &light_hit) > 0.0f) {
            auto &mat = scene_->GetEntity(light_hit.hit_entity_id)
                            .GetMaterial(light_hit.material_id);
            radiance += throughput * albedo *
                        fmax(0.f, glm::dot(normal, direction)) * w *
                        mat.emission * mat.emission_strength / Light_pdf;

          } else {
            radiance += throughput * albedo *
                        fmax(0.f, glm::dot(normal, direction)) *
                        w * glm::vec3{scene_->SampleEnvmap(ray.direction())} / Light_pdf;
          }
          direction = Lambert.Generate(origin, ray.time(), rd);
          next_ray = Ray(origin, direction, ray.time());
          Lamb_pdf = Lambert.Value(next_ray);
          Light_pdf = Light.Value(next_ray);
          w = square(Lamb_pdf) / (square(Lamb_pdf) + square(Light_pdf));
          throughput *= albedo * fmax(0.f, glm::dot(normal, direction)) * w / Lamb_pdf;
          /*
          MixturePdf Gen(&Lambert, &Light, 0.5f);
          direction = Gen.Generate(origin, ray.time(), rd);
          next_ray = Ray(origin, direction, ray.time());
          float pdf = Gen.Value(next_ray);
          float scatter = std::max(0.f, glm::dot(normal, direction) * INV_PI);
          if (pdf < 1e-9)
              break;
          throughput *= albedo * scatter / pdf;
          if (glm::dot(normal, direction) < 0)
            break;
          */
      } else if (material.material_type == MATERIAL_TYPE_MMD) {
          return glm::vec3(0.2f);
      } else if (material.material_type == MATERIAL_TYPE_SPECULAR) {
          // You should change this into Cook-Torrence.
          direction = glm::reflect(ray.direction(), normal);
          float cos_theta = fmin(glm::dot(-ray.direction(), normal), 1.0f);
          throughput *= material.FresnelSchlick(albedo, cos_theta);
          if (glm::dot(normal, direction) < 0)
            break;
          next_ray = Ray(origin, direction, ray.time());
      } else if (material.material_type == MATERIAL_TYPE_TRANSMISSIVE) {
          // Assume all lights have same index of refraction
          float refract_ratio = hit_record.front_face ? (1.0 / material.IOR)
                                                      : material.IOR;
          if (material.thin) {
            float reflect_ratio = 1 - refract_ratio;
            reflect_ratio = (1 - reflect_ratio) * (1 - reflect_ratio) *
                            reflect_ratio / (1 - reflect_ratio * reflect_ratio);
            float cos_theta = fmin(glm::dot(-ray.direction(), normal), 1.0f);
            if (reflect_ratio > RandomProb(rd)) {
              direction = glm::reflect(ray.direction(), normal);
              throughput *= material.FresnelSchlick(albedo, cos_theta);
            } else {
              direction = ray.direction();
              throughput *= albedo;
            }
          } else {
              float cos_theta = fmin(glm::dot(-ray.direction(), normal), 1.0f);
              float f0 = (1 - refract_ratio) / (1 + refract_ratio);
              f0 *= f0;
              direction = glm::refract(ray.direction(), normal, refract_ratio);
              if (glm::length(direction) == 0.f ||
                  material.FresnelSchlick(f0, cos_theta) >
                      RandomProb(rd)) {
                direction = glm::reflect(ray.direction(), normal);
                throughput *= material.FresnelSchlick(albedo, cos_theta);
              } else {
                throughput *= albedo;
              }
          }
        direction = glm::normalize(direction);
        next_ray = Ray(origin, direction, ray.time());
      } else if (material.material_type == MATERIAL_TYPE_MEDIUM) {
        UniformSpherePdf Medium(normal);
        MixturePdf Gen(&Medium, &Light, 0.1f);
        direction = glm::normalize(Gen.Generate(origin, ray.time(), rd));
        next_ray = Ray(origin, direction, ray.time());
        throughput *= albedo / Gen.Value(next_ray);
      } else if (material.material_type == MATERIAL_TYPE_PRINCIPLED) {
        const bool UseBSDFSampler = false;
        if (UseBSDFSampler) {
            // decide which lobe to sample
            float metallicBRDF = material.metallic;
            float specularBSDF = (1.0f - material.metallic) * material.specTrans;
            float dielectricBRDF =
                (1.0f - material.specTrans) * (1.0f - material.metallic);
            float specularWeight = metallicBRDF + dielectricBRDF;
            float transmissionWeight = specularBSDF;
            float diffuseWeight = dielectricBRDF;
            float clearcoatWeight = 1.0f * std::clamp(material.clearcoat, 0.0f, 1.0f);
            float norm = 1.0f / (specularWeight + transmissionWeight +
                                 diffuseWeight + clearcoatWeight);
            float pSpecular = specularWeight * norm;
            float pSpecTrans = transmissionWeight * norm;
            float pDiffuse = diffuseWeight * norm;
            float pClearcoat = clearcoatWeight * norm;
            SampleDisneyBRDFPdf sampleBRDF(normal, material, pSpecular);
            SampleDisneyClearCoatPdf sampleClearCoat(normal, material,
                                                    pClearcoat);
            SampleDisneyDiffusePdf sampleDiffuse(normal, material, pDiffuse);
            SampleDisneySpecTransPdf sampleSpecTrans(normal, material,
                                                    pSpecTrans);
            float pLobe = 0.0f;

                std::uniform_real_distribution<float> dist(0.0f, 1.0f);
                float p = dist(rd);
                if (p <= pSpecular) {
                direction = sampleBRDF.Generate(-ray.direction(), ray.time(), rd);
                  pLobe = pSpecular;
                } else if (p > pSpecular && p <= (pSpecular + pClearcoat)) {
                direction = sampleClearCoat.Generate(-ray.direction(), ray.time(), rd);
                  pLobe = pClearcoat;
                } else if (p > pSpecular + pClearcoat &&
                           p <= (pSpecular + pClearcoat + pDiffuse)) {
                direction = sampleDiffuse.Generate(-ray.direction(), ray.time(), rd);
                  pLobe = pDiffuse;
            } else if (pSpecTrans >= 0.0f) {
                direction  = sampleSpecTrans.Generate(-ray.direction(), ray.time(), rd);
                  pLobe = pSpecTrans;
            } else {
                // just break
                break;
            }
            
            if (glm::length(direction) < 1e-9)
                  break;
            float refract_ratio =
                hit_record.front_face ? (1.0 / material.IOR) : material.IOR;
            direction = glm::normalize(direction);
            next_ray = Ray(origin, direction, ray.time());
            throughput *= material.EvaluateDisney(-ray.direction(), next_ray.direction(),
                                        normal, albedo, refract_ratio) / pLobe;
        } else {
          UniformSpherePdf USP(normal, tangent);
          MixturePdf Gen(&USP, &Light, 0.5);
          direction = Gen.Generate(origin, ray.time(), rd);
          float refract_ratio =
              hit_record.front_face ? (1.0 / material.IOR) : material.IOR;
          direction = glm::normalize(direction);
          next_ray = Ray(origin, direction, ray.time());
          throughput *=
              material.EvaluateDisney(-ray.direction(), next_ray.direction(),
                                                normal, albedo, refract_ratio);
          throughput /= Gen.Value(next_ray);
        }
      }
      ray = next_ray;
      throughput *= INV_RR;
    } else {
      radiance += throughput * glm::vec3{scene_->SampleEnvmap(ray.direction())};
      break;
    }
  }
  return glm::min(radiance, glm::vec3(1));
}
}  // namespace sparks
