#pragma once
#include "random"
#include "sparks/assets/scene.h"
#include "sparks/assets/ray.h"
#include "sparks/renderer/renderer_settings.h"
#include "sparks/assets/accelerated_mesh.h"

namespace sparks {
class PathTracer {
 public:
  PathTracer(const RendererSettings *render_settings, const Scene *scene);
  [[nodiscard]] glm::vec3 SampleRay(Ray ray,
                                    int x,
                                    int y,
                                    int sample) const;

 private:
  const RendererSettings *render_settings_{};
  const Scene *scene_{};
};
}  // namespace sparks
