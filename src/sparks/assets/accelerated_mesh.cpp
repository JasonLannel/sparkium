#include "sparks/assets/accelerated_mesh.h"

#include "algorithm"

namespace sparks {
AcceleratedMesh::AcceleratedMesh(const Mesh &mesh) : Mesh(mesh) {
  BuildAccelerationStructure();
}

AcceleratedMesh::AcceleratedMesh(const std::vector<Vertex> &vertices,
                                 const std::vector<uint32_t> &indices)
    : Mesh(vertices, indices) {
  BuildAccelerationStructure();
}

float AcceleratedMesh::TraceRay(const glm::vec3 &origin,
                                const glm::vec3 &direction,
                                float t_min,
                                HitRecord *hit_record) const {
  return Mesh::TraceRay(origin, direction, t_min, hit_record);
}

void AcceleratedMesh::BuildAccelerationStructure() {
    // Step 1: Build the bounding volume hierarchy
    if (indices_.size() == 0 || indices_.size() == 1) {
		return;
	}
    if (bvh_nodes_.size() == 0) {
        bvh_nodes_.resize(indices_.size() * 2 - 1);
    }
    std::vector<int> entity_indices(indices_.size());
    for (int i = 0; i < entity_indices.size(); i++) {
	  entity_indices[i] = i;
	}
    BuildTree(entity_indices, 0, entity_indices.size() - 1, 0);
}

void AcceleratedMesh::BuildTree(std::vector<int>& entity_indices, int start_index, int end_index, int node_index) {
  TreeNode& node = bvh_nodes_[node_index];
  
  // Step 1: Set the range of entity indices for this node
  node.start_index = start_index;
  node.end_index = end_index;

  // Step 2: Check if this node is a leaf
  if (end_index - start_index <= 3) {
    node.is_leaf = true;
    return;
  }

  // Step 3: Calculate the bounding box of all entities in the range
  AxisAlignedBoundingBox aabb(
    std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(),
    std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(),
    std::numeric_limits<float>::max(), -std::numeric_limits<float>::max()
  );
  
  for (int i = start_index; i <= end_index; i++) {
    int entity_index = entity_indices[i];
    const Vertex& vertex = vertices_[indices_[entity_index]];
    const glm::vec3& position = vertex.position;

    // Expand the bounding box to include the vertex position
    aabb.x_low = std::min(aabb.x_low, position.x);
    aabb.x_high = std::max(aabb.x_high, position.x);
    aabb.y_low = std::min(aabb.y_low, position.y);
    aabb.y_high = std::max(aabb.y_high, position.y);
    aabb.z_low = std::min(aabb.z_low, position.z);
    aabb.z_high = std::max(aabb.z_high, position.z);
  }

  // Step 4: Split the range based on the largest axis of the bounding box
  glm::vec3 box_size = glm::vec3(
    aabb.x_high - aabb.x_low,
    aabb.y_high - aabb.y_low,
    aabb.z_high - aabb.z_low
  );

  int largest_axis = (box_size.x > box_size.y) ? ((box_size.x > box_size.z) ? 0 : 2) : ((box_size.y > box_size.z) ? 1 : 2);
  
  auto compare_fn = [largest_axis, &entity_indices, this](int a, int b) {
    const Vertex& vertex_a = vertices_[indices_[entity_indices[a]]];
    const Vertex& vertex_b = vertices_[indices_[entity_indices[b]]];
    return vertex_a.position[largest_axis] < vertex_b.position[largest_axis];
  };
  
  std::sort(entity_indices.begin() + start_index, entity_indices.begin() + end_index + 1, compare_fn);

  int mid_index = (start_index + end_index) / 2;
  node.left_child_index = node_index + 1;
  node.right_child_index = node_index + 2;

  // Step 5: Build left and right child nodes recursively
  BuildTree(entity_indices, start_index, mid_index, node.left_child_index);
  BuildTree(entity_indices, mid_index + 1, end_index, node.right_child_index);
}

}  // namespace sparks
