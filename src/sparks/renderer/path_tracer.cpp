#include "sparks/renderer/path_tracer.h"

#include "sparks/util/util.h"
#include "sparks/assets/bsdf.h"
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
  glm::vec3 throughput(1.0f);
  glm::vec3 radiance(0.0f);
  glm::vec3 origin, direction, normal, tangent, albedo;
  HitRecord hit_record;
  Ray next_ray(ray);
  const float RR = 0.95, INV_RR = 1.0 / RR;
  Material medium_pre = Material();
  for (int bounce = 0, max_depth = render_settings_->num_bounces; bounce < max_depth;
       ++bounce) {
    bool medium_hit = false;
    float tmax = 1e10f;
    if (medium_pre.material_type == MATERIAL_TYPE_MEDIUM) {
      std::uniform_real_distribution<float> randomProb(0.0f, 1.0f);
      float hit_dis = -log(fmax(randomProb(rd), 1e-10)) / medium_pre.density;
      medium_hit = true;
      hit_record.hit_entity_id = -1;
      hit_record.position = ray.at(hit_dis);
      hit_record.normal = -ray.direction();
      tmax = fmin(tmax, hit_dis);
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
      if (mat.material_type == MATERIAL_TYPE_MMD)
          return glm::vec3(0.2f);
      if (hit_record.material_id > 0 && RandomProb(rd) > mat.alpha) {
        //Alpha Shadow
        ray = Ray(hit_record.position, ray.direction(), ray.time());
        --bounce;
        continue;
      }
      origin = hit_record.position;
      normal = hit_record.normal;
      tangent = hit_record.tangent;
      if (glm::dot(normal, ray.direction()) > 0)
        normal = -normal;
      std::unique_ptr<bsdf> mat_bsdf;
      switch (mat.material_type) { 
        case MATERIAL_TYPE_LAMBERTIAN:
          mat_bsdf = std::make_unique<Lambertian>();
          break;
        case MATERIAL_TYPE_SPECULAR:
          mat_bsdf = std::make_unique<Specular>();
          break;
        case MATERIAL_TYPE_TRANSMISSIVE:
          mat_bsdf = std::make_unique<Transmissive>();
          break;
        case MATERIAL_TYPE_MEDIUM:
          mat_bsdf = std::make_unique<Medium>();
          break;
        case MATERIAL_TYPE_PRINCIPLED:
        default:
          mat_bsdf = std::make_unique<Principled>();
          if (hit_record.front_face) {
            mat.IOR = 1.0f / mat.IOR;
          } else {
            normal = -normal;
          }
          break;
      }
      float light_pdf(0), bsdf_pdf(0), w(0);
      glm::vec3 reflectance(0); 
      glm::vec3 emission(0);
      if (mat.material_type != MATERIAL_TYPE_SPECULAR &&
          mat.material_type != MATERIAL_TYPE_TRANSMISSIVE) {
        direction =
            scene_->SampleLight(origin, ray.time(), rd, &light_pdf, &emission);
        if (light_pdf > 0.0f) {
          reflectance = mat_bsdf->evaluate(-ray.direction(), direction, origin,
                                            normal, tangent, mat, &bsdf_pdf);
          w = square(light_pdf) / (square(light_pdf) + square(bsdf_pdf));
          radiance += throughput * reflectance * emission * w / light_pdf;
        }
      }
      direction = mat_bsdf->sample(-ray.direction(), origin, normal, tangent,
                                   mat, rd, &bsdf_pdf);
      if (bsdf_pdf == 0.0f)
        break;
      next_ray = Ray(origin, direction, ray.time());
      if (mat.material_type != MATERIAL_TYPE_SPECULAR &&
          mat.material_type != MATERIAL_TYPE_TRANSMISSIVE) {
        light_pdf = scene_->LightValue(next_ray);
        w = square(bsdf_pdf) / (square(light_pdf) + square(bsdf_pdf));
        throughput *= mat_bsdf->evaluate(-ray.direction(), direction, origin,
                                         normal, tangent, mat, nullptr) *
                      w / bsdf_pdf;
      } else {
        throughput *= mat_bsdf->evaluate(-ray.direction(), direction, origin,
                                         normal, tangent, mat, nullptr);
      }
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
