#ifndef SPARKS_RAY_H
#define SPARKS_RAY_H

#include "glm/glm.hpp"

namespace sparks {

class Ray {
 public:

  Ray(const glm::vec3 &origin, const glm::vec3 &direction)
      : orig_(origin), dir_(direction), time_(0.0) {
  }

  Ray(const glm::vec3 &origin, const glm::vec3 &direction, double time)
      : orig_(origin), dir_(direction), time_(time) {}

  glm::vec3 origin() const {
    return orig_;
  }
  glm::vec3 direction() const {
    return dir_;
  }
  double time() const {
    return time_;
  }

  glm::vec3 at(double t) const {
    return orig_ + dir_ * glm::vec3(t);
  }

 private:
  glm::vec3 orig_{glm::vec3(0.0f)};
  glm::vec3 dir_{glm::vec3(0.0f)};
  double time_{0.0};
  // RayColor color_;
};

}  // namespace sparks

#endif  // SPARKS_RAY_H
