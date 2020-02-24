#ifndef LSCM_JJ_H
#define LSCM_JJ_H

#include <vector>
#include <tuple>
#include <Eigen/Core>
#include <array>
#include <Eigen/SparseCore>


class HalfEdgeMesh;

std::tuple<std::vector<double>, std::vector<Eigen::MatrixXd>> get_triangle_area_and_transform(const HalfEdgeMesh &he_mesh);


template<typename Derived>
std::tuple<double, Eigen::MatrixXd> get_triangle_transform(const Eigen::MatrixBase<Derived> &tri);

std::vector<Eigen::Vector2d> solver(int fix_v1, int fix_v2, const std::vector<double> &tri_area, const std::vector<Eigen::MatrixXd> &tri_transform, const HalfEdgeMesh &he_mesh);


template<typename Derived>
std::vector<Eigen::Vector2d> get_vert_parameter(int fix_v1, const std::array<double, 2> &uv_v1, int fix_v2, const std::array<double, 2> &uv_v2, const Eigen::MatrixBase<Derived> &v);


template<typename DerivedM, typename Derivedb>
void add_triangle_to_constrain(
  int fix_v1, int fix_v2, int &row, double area,
  const Eigen::MatrixBase<DerivedM> &trans_u,
  Eigen::MatrixBase<Derivedb> &b,
  const std::array<size_t, 3> &tri_v,
  const std::array<double, 2> &uv_v1,
  const std::array<double, 2> &uv_v2,
  std::vector<Eigen::Triplet<double>> &equation_coeff);

int convert_vert_id_to_var_id(int v, int fix_v1, int fix_v2);

void scale_parameter_to_unit(std::vector<Eigen::Vector2d> &para);

#endif // LSCM_JJ_H
