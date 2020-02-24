#ifndef LSCM_JJ_H
#define LSCM_JJ_H

#include <vector>
#include <tuple>
#include <Eigen/Core>


class HalfEdgeMesh;

std::tuple<std::vector<double>, std::vector<Eigen::MatrixXd>> get_triangle_area_and_transform(const HalfEdgeMesh &he_edge);




#endif // LSCM_JJ_H
