#include "sparks/renderer/path_tracer.h"

#include "sparks/util/util.h"
#include <glm/ext/scalar_constants.hpp>

// TODO: implement a path tracing algorithm that could handle diffusive material and specular material correctly with a proper acceleration structure

namespace sparks {
PathTracer::PathTracer(const RendererSettings *render_settings,
                       const Scene *scene) {
  render_settings_ = render_settings;
  scene_ = scene;
  // initialize your acceleration structure here.
  // getMesh() function has not been implemented yet.
  // for (const auto& entity : scene_->GetEntities()) {
  //     if (entity.GetMesh().GetVertices().empty()) {
	//   continue;
	// }
	// accelerated_meshes_.emplace_back(entity.GetMesh());
	// accelerated_meshes_.back().BuildAccelerationStructure();
  // }
}

glm::vec3 PathTracer::RandomUnitVector(std::mt19937& rd) {
  std::uniform_real_distribution<float> dist(0.0f, 1.0f);
  float theta = dist(rd) * 2.0f * glm::pi<float>();
  float phi = dist(rd) * glm::pi<float>();
  return glm::vec3{std::cos(theta) * std::sin(phi), std::sin(theta) * std::sin(phi), std::cos(phi)};
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
      } else if (material.material_type == MATERIAL_TYPE_LAMBERTIAN) {
          // diffusive material
        glm::vec3 albedo = material.albedo_color;
        if (material.albedo_texture_id >= 0) {
          albedo *= glm::vec3{
              scene_->GetTextures()[material.albedo_texture_id].Sample(
                  hit_record.tex_coord)};
        }
        glm::vec3 new_direction = glm::normalize( hit_record.normal + glm::vec3{ RandomUnitVector(rd) });
        throughput *= albedo;
        origin = hit_record.position;
        direction = new_direction;
        radiance += throughput * scene_->GetEnvmapMinorColor();
        throughput *= std::abs(glm::dot(direction, hit_record.normal)); 
      } else if (material.material_type == MATERIAL_TYPE_SPECULAR) {
          // specular material
        glm::vec3 albedo = material.albedo_color;
        if (material.albedo_texture_id >= 0) {
          albedo *= glm::vec3{
              scene_->GetTextures()[material.albedo_texture_id].Sample(
                  hit_record.tex_coord)};
        }
        glm::vec3 new_direction = glm::reflect(direction, hit_record.normal);
        throughput *= albedo;
        origin = hit_record.position;
        direction = new_direction;
        radiance += throughput * scene_->GetEnvmapMinorColor();
        throughput *= std::abs(glm::dot(direction, hit_record.normal));
      } else {
        throughput *=
            material.albedo_color *
            glm::vec3{scene_->GetTextures()[material.albedo_texture_id].Sample(
                hit_record.tex_coord)};
        origin = hit_record.position;
        direction = scene_->GetEnvmapLightDirection();
        radiance += throughput * scene_->GetEnvmapMinorColor();
        throughput *=
            std::max(glm::dot(direction, hit_record.normal), 0.0f) * 2.0f;
        if (scene_->TraceRay(origin, direction, 1e-3f, 1e4f, nullptr) < 0.0f) {
          radiance += throughput * scene_->GetEnvmapMajorColor();
        }
      }

      // Russian roulette termination
      float rr_probability = 0.5f;
      assert (rr_probability >= 0.0f && rr_probability <= 1.0f);
      if (bounce > 5) {
        float rr_random = std::uniform_real_distribution<float>(0.0f, 1.0f)(rd);
        if (rr_random > rr_probability) {
          break;
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
