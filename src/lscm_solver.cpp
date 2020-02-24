#include "lscm.h"
#include "half_edge.h"


using namespace std;
using namespace Eigen;



std::tuple<std::vector<double>, std::vector<Eigen::MatrixXd>> get_triangle_area_and_transform(const HalfEdgeMesh &he_edge)
{
  const size_t tri_num = he_edge.get_face_num();
  vector<double> tri_area(tri_num);
  vector<MatrixXd> tri_transform(tri_num);
  for (size_t i = 0; i < tri_num; ++i)
  {
    Matrix3d tri = he_mesh.get_tri(i);
    tie{tri_area[i], tri_transform[i]} = get_triangle_transform(tri);
  }
  return {tri_area, tri_transform};
}

template<typename Derived>
std::tuple<double, Eigen::MatrixXd> get_triangle_transform(const Eigen::MatrixBase<Derived> &tri)
{
  MatrixXd basis_coord(2, 3);
  MatrixXd trans(2, 3);

  Vector3d X = (tri.col(1)-tri.col(0)).normalized();
  Vector3d norm = X.cross(tri.col(2)-tri.col(0)).normalized();
  Vector3d Y = norm.cross(X);
  basis_coord.col(0) = {0, 0};
  basis_coord.col(1) = {(tri.col(1)-tri.col(0)).norm(), 0};
  basis_coord(0, 2) = X.dot(tri.col(2) - tri.col(0));
  basis_coord(1, 2) = Y.dot(tri.col(2) - tri.col(0));
  double area = 0.5 * (tri.col(1)-tri.col(0)).cross(tri.col(2)-tri.col(0)).norm();
  trans(0, 0) = (basis_coord(1, 1) - basis_coord(1, 2)) / 2.0 / area;
  trans(0, 1) = (basis_coord(1, 2) - basis_coord(1, 0)) / 2.0 / area;
  trans(0, 2) = (basis_coord(1, 0) - basis_coord(1, 1)) / 2.0 / area;
  trans(1, 0) = (basis_coord(0, 2) - basis_coord(0, 1)) / 2.0 / area;
  trans(1, 1) = (basis_coord(0, 0) - basis_coord(0, 2)) / 2.0 / area;
  trans(1, 2) = (basis_coord(0, 1) - basis_coord(0, 0)) / 2.0 / area;

  return {area, trans};
}


Eigen::VectorXd solver(int fix_v1, int fix_v2, const std::vector<double> &tri_area, const std::vector<Eigen::MatrixXd> &tri_transform)
{
  array<double, 2> uv_v1 = {0.0, 0.0};
  array<double, 2> uv_v2 = {1.0, 1.0};
  vector<Triplet<double>> equation_coeff;

  const int var_num = 2 * (he_mesh.get_vert_num() - 2);
  VectorXd b = VectorXd::Zero(var_num);
  const size_t tri_num = he_edge.get_face_num();
  int row = 0;
  for (size_t i = 0; i < tri_num; ++i)
  {
    array<size_t, 3> tri_v = he_mesh.get_tri_vert_id(i);
    const MatrixXd &trans_u = tri_transform[i];
    add_triangle_to_constrain(fix_v1, fix_v2, tri_area[i], tri_transform[i], b, tri_v, uv_v1, uv_v2, equation_coeff);
  }  


  SparseMatrix<double> A;
  A.setFromTriplets(equation.begin(), equation.end());
  SparseQR<SparseMatrix<double>> qr_solver;
  solver.compute(A);
  if(solver.info()!=Success)
  {
    cerr << "decomposition failed" << endl;
    return VectorXd(0);
  }
  VectorXd x = solver.solve(b);
  if(solver.info()!=Success)
  {
    cerr << "solving failed" << endl;
    return VectorXd(0);
  }
  return x;
}

template<typename DerivedM, typename Derivedb>
void add_triangle_to_constrain(
  int fix_v1, int fix_v2, double area,
  const Eigen::MatrixBase<DerivedM> &trans_u,
  const Eigen::MatrixBase<Derivedb> &b,
  const std::array<size_t, 3> &tri_v,
  const std::array<double, 2> &uv_v1,
  const std::array<double, 2> &uv_v2,
  std::vector<Eigen::Triplet<double>> &equation_coeff)
{
  Matrix2d T;
  T << 0, -1, 1, 0;
  MatrixXd trans_v = T * trans_u;
  double c_area = sqrt(area);
  for (int r = 0; r < 2; ++r)
  {
    for (int v = 0; v < 3; ++v)
    {
      if (tri_v[v] == fix_v1)
        b[row] -= c_area * (trans_u(r, v) * uv_v1[0] - trans_v(r, v) * uv_v1[1]);
      else if (tri_v[v] == fix_v2)
        b[row] -= c_area * (trans_u(r, v) * uv_v2[0] - trans_v(r, v) * uv_v2[1]);
      else
      {
        int col = convert_vert_id_to_var_id(tri_v[v], fix_v1, fix_v2);
        equation_coeff.push_back(Triplet(row, 2 * col, c_area * trans_u(r, v)));
        equation_coeff.push_back(Triplet(row, 2*col+1, -c_area * trans_v(r, v)));
      }
    }
    ++row;
  }
}


inline int convert_vert_id_to_var_id(int v, int fix_v1, int fix_v2)
{
  if (v < min(fix_v1, fix_v2))
    return v;
  if (v < max(fix_v1, fix_v2))
    return v - 1;
  return v - 2;
}
