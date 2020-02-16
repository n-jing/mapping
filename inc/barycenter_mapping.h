#ifndef BARYCENTER_MAPPING_JJ_H
#define BARYCENTER_MAPPING_JJ_H

#include <Eigen/Core>
#include <unordered_set>

class HalfEdgeMesh;


Eigen::Vector2d get_boundary_parameter_coordinate(double t);

std::vector<Eigen::Vector2d> barycenter_mapping(const HalfEdgeMesh &he_mesh, const std::vector<size_t> &loop);

int write_parameter_domain(const HalfEdgeMesh &he_mesh, const std::vector<Eigen::Vector2d> &parameter_vert);


Eigen::VectorXd get_inner_vert_parameter_coordinate(const HalfEdgeMesh &he_mesh, const std::unordered_set<size_t> &boundary_vert, const std::unordered_map<size_t, size_t> &vert_to_inner_id, const std::vector<Eigen::Vector2d> &parameter_vert, int u);

#endif // BARYCENTER_MAPPING_JJ_H
