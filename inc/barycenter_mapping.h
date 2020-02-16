#ifndef BARYCENTER_MAPPING_JJ_H
#define BARYCENTER_MAPPING_JJ_H

#include <Eigen/Core>

class HalfEdgeMesh;


Eigen::Vector2d get_boundary_parameter_coordinate(double t);

std::vector<Eigen::Vector2d> barycenter_mapping(const HalfEdgeMesh &he_mesh, const std::vector<size_t> &loop);

int write_parameter_domain(const HalfEdgeMesh &he_mesh, const std::vector<Eigen::Vector2d> &parameter_vert);



#endif // BARYCENTER_MAPPING_JJ_H
