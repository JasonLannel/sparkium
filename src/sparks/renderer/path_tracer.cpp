#include "sparks/renderer/path_tracer.h"

#include "sparks/util/util.h"
#include "sparks/assets/bsdf.h"
#include "sparks/assets/pdf.h"
#include <glm/ext/scalar_constants.hpp>

#include <time.h>

constexpr auto USE_MIS = true;
constexpr auto USE_POWER_HEURISTIC = true;

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
  //Russian Roulette
  const float RR = 0.95, INV_RR = 1.0 / RR;
  //Color
  glm::vec3 throughput(1.0f);
  glm::vec3 radiance(0.0f);
  //Ray, Hit Record
  glm::vec3 origin, direction, normal, tangent, albedo;
  HitRecord hit_record;
  Ray next_ray(ray);
  //Current Medium
  Material medium_pre = Material();
  bsdf_handler bsdf_op;
  for (int bounce = 0, max_depth = render_settings_->num_bounces; bounce < max_depth;
       ++bounce) {
    bool medium_hit = false;
    float tmax = 1e5f;
    if (medium_pre.material_type == MATERIAL_TYPE_MEDIUM) {
      std::uniform_real_distribution<float> randomProb(0.0f, 1.0f);
      float hit_dis = -log(fmax(randomProb(rd), 1e-10)) / medium_pre.sigma;
      if (hit_dis < tmax) {
          medium_hit = true;
          hit_record.hit_entity_id = -1;
          hit_record.position = ray.at(hit_dis);
          hit_record.normal = -ray.direction();
          tmax = hit_dis;
      }
    }
    if (scene_->TraceRay(ray, 1e-3f, tmax,
                         &hit_record) > 0.0f || medium_hit) {
      Material mat;
      if (hit_record.hit_entity_id == -1) {
        mat = medium_pre;
      } else {
        mat = scene_->GetEntity(hit_record.hit_entity_id)
                  .GetMaterial(hit_record.material_id);
        scene_->LoadTextureForMaterial(mat, hit_record);
        medium_hit = false;
      }
      // Handle Medium Change
      if (!medium_hit && mat.material_type == MATERIAL_TYPE_MEDIUM) {
        if (hit_record.front_face) {
          medium_pre = mat;
        } else {
          medium_pre = Material();  //vaccum;
        }
        ray = Ray(hit_record.position, ray.direction(), ray.time());
        --bounce;
        continue;
      }
      radiance += throughput * mat.emission * mat.emission_strength;
      if (mat.material_type == MATERIAL_TYPE_EMISSION)
        break;
      if (hit_record.material_id > 0 && RandomProb(rd) > mat.alpha) {
        //Alpha Shadow
        ray = Ray(hit_record.position, ray.direction(), ray.time());
        --bounce;
        continue;
      }
      origin = hit_record.position;
      normal = hit_record.normal;
      tangent = hit_record.tangent;
      float light_pdf(0), bsdf_pdf(0), w(0);
      glm::vec3 reflectance(0); 
      glm::vec3 emission(0);
      // Light source sample
      if (mat.material_type != MATERIAL_TYPE_SPECULAR &&
          mat.material_type != MATERIAL_TYPE_TRANSMISSIVE) {
        direction =
            scene_->SampleLight(origin, ray.time(), rd, medium_pre, light_pdf, emission);
        if (light_pdf > 0.0f) {
          if constexpr (USE_MIS && USE_POWER_HEURISTIC) {
            reflectance = bsdf_op.evaluate(-ray.direction(), direction, origin,
                                         normal, tangent, mat, bsdf_pdf);
            w = square(light_pdf) / (square(light_pdf) + square(bsdf_pdf));
          } else if constexpr (USE_MIS) {
            reflectance = bsdf_op.evaluate(-ray.direction(), direction, origin,
                                         normal, tangent, mat, bsdf_pdf);
            w = square(light_pdf) / (square(light_pdf) + square(bsdf_pdf));
          } else {
            w = 0.5f;
          }
          radiance += throughput * reflectance * emission * w / light_pdf;
        }
      }
      // BSDF sample
      direction = bsdf_op.sample(-ray.direction(), origin, normal, tangent,
                                   mat, rd, bsdf_pdf, reflectance);
      if (bsdf_pdf == 0.0f)
        break;
      next_ray = Ray(origin, direction, ray.time());
      if (mat.material_type != MATERIAL_TYPE_SPECULAR &&
          mat.material_type != MATERIAL_TYPE_TRANSMISSIVE) {
        if constexpr (USE_MIS && USE_POWER_HEURISTIC) {
          light_pdf = scene_->LightValue(next_ray);
          w = square(bsdf_pdf) / (square(light_pdf) + square(bsdf_pdf));
        } else if constexpr (USE_MIS) {
          light_pdf = scene_->LightValue(next_ray);
          w = bsdf_pdf / (light_pdf + bsdf_pdf);
        } else {
          w = 0.5f;
        }
      } else {
        w = 1.0f;
      }
      throughput *= reflectance * w;
      ray = next_ray;
      if (bounce > 7) {
        if (RandomProb(rd) > RR) {
          break;
        }
        throughput *= INV_RR;
      }
    } else {
      radiance += throughput * glm::vec3{scene_->SampleEnvmap(ray.direction())};
      break;
    }
  }
  return glm::min(radiance, glm::vec3(1));
}
}  // namespace sparks
