#pragma once
#include "random"
#include "sparks/assets/scene.h"
#include "sparks/renderer/renderer_settings.h"
#include "sparks/assets/accelerated_mesh.h"

// TODO: implement a path tracing algorithm that could handle diffusive material and specular material correctly with a proper acceleration structure

namespace sparks {
class PathTracer {
 public:
  PathTracer(const RendererSettings *render_settings, const Scene *scene);
  [[nodiscard]] glm::vec3 SampleRay(glm::vec3 origin,
                                    glm::vec3 direction,
                                    int x,
                                    int y,
                                    int sample) const;

 private:
  const RendererSettings *render_settings_{};
  const Scene *scene_{};
  // define a random unit vector generator here.
  static glm::vec3 RandomUnitVector(std::mt19937 &rd);

  // You can add your acceleration structure here.
  // std::vector<AcceleratedMesh> accelerated_meshes_{};
};
}  // namespace sparks
