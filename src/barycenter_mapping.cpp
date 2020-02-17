#include "barycenter_mapping.h"
#include "half_edge.h"
#include "write_to_file.h"
#include "save_matrix_as_image.h"
#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>


#include <unsupported/Eigen/SparseExtra>


using namespace std;
using namespace Eigen;


Eigen::Vector2d get_boundary_parameter_coordinate(double t)
{
  if (t < 0.25)
    return {t / 0.25, 0};
  else if (t < 0.5)
    return {1, (t-0.25) / 0.25};
  else if (t < 0.75)
    return {(0.75 - t) / 0.25, 1};
  else
    return {0, (1 - t) / 0.25};
}

std::vector<Eigen::Vector2d> barycenter_mapping(const HalfEdgeMesh &he_mesh, const std::vector<size_t> &loop)
{
  const auto he_e = he_mesh.get_half_edge();
  double loop_length = accumulate(
    loop.begin(), loop.end(), 0.0, [&he_e, &he_mesh] (double b,size_t e)->double
    {
      size_t ne = he_e.at(he_e.at(e).next_e).next_e;
      array<size_t, 2> he_v = {he_e.at(e).he_v, he_e.at(ne).he_v};
      double le = (he_mesh.get_vert(he_v[0]) - he_mesh.get_vert(he_v[1])).norm();
      return b + le;
    });

  const size_t vert_num = he_mesh.get_vert_num();
  vector<Vector2d> parameter_vert(vert_num);
  parameter_vert[0] = {0, 0};

  unordered_set<size_t> boundary_vert;
  double length = 0;
  for (auto e : loop)
  {
    size_t ne = he_e.at(he_e.at(e).next_e).next_e;
    array<size_t, 2> he_v = {he_e.at(ne).he_v, he_e.at(e).he_v};
    double le = (he_mesh.get_vert(he_v[0]) - he_mesh.get_vert(he_v[1])).norm();
    length += le;
    double t = length / loop_length;
    parameter_vert[he_v[1]] = get_boundary_parameter_coordinate(t);
    boundary_vert.insert(he_v[0]);
  }

  unordered_map<size_t, size_t> vert_to_inner_id;
  unordered_map<size_t, size_t> inner_id_to_vert;
  for (size_t i = 0; i < vert_num; ++i)
  {
    if (boundary_vert.count(i))
      continue;
    vert_to_inner_id.emplace(i, vert_to_inner_id.size());
    inner_id_to_vert.emplace(vert_to_inner_id.size() - 1, i);
  }

  VectorXd u = get_inner_vert_parameter_coordinate(he_mesh, boundary_vert, vert_to_inner_id, parameter_vert, 0);
  VectorXd v = get_inner_vert_parameter_coordinate(he_mesh, boundary_vert, vert_to_inner_id, parameter_vert, 1);

  const size_t inner_vert_num = vert_to_inner_id.size();
  for (size_t i = 0; i < inner_vert_num; ++i)
  {
    size_t vert_id = inner_id_to_vert.at(i);
    parameter_vert[vert_id][0] = u[i];
    parameter_vert[vert_id][1] = v[i];
  }

  return parameter_vert;
}

Eigen::VectorXd get_inner_vert_parameter_coordinate(const HalfEdgeMesh &he_mesh, const std::unordered_set<size_t> &boundary_vert, const std::unordered_map<size_t, size_t> &vert_to_inner_id, const std::vector<Eigen::Vector2d> &parameter_vert, int u)
{
  const size_t vert_num = he_mesh.get_vert_num();
  const size_t inner_vert_num = vert_to_inner_id.size();

  VectorXd b = VectorXd::Zero(inner_vert_num, 1);
  vector<Triplet<double>> triplet_list;
  for (size_t i = 0; i < vert_num; ++i)
  {
    if (boundary_vert.count(i))
      continue;
    const size_t inner_id = vert_to_inner_id.at(i);
    triplet_list.emplace_back(inner_id, inner_id, 1);
    vector<size_t> neigh_vert = he_mesh.get_vert_one_ring(i);
    const size_t neigh_vert_num = neigh_vert.size();
    for (auto v : neigh_vert)
    {
      if (boundary_vert.count(v))
      {
        b[inner_id] += parameter_vert[v][u] / neigh_vert_num;
        continue;
      }
      const size_t neigb_inner_id = vert_to_inner_id.at(v);
      triplet_list.emplace_back(inner_id, neigb_inner_id, -1.0 / neigh_vert_num);
    }
  }

  SparseMatrix<double> A(inner_vert_num, inner_vert_num);
  A.setFromTriplets(triplet_list.begin(), triplet_list.end());
  A.makeCompressed();

  SparseLU<SparseMatrix<double>> solver;
  solver.compute(A);
  if (solver.info() != Success)
  {
    cerr << "decomposition failed" << endl;
  }
  VectorXd x = solver.solve(b);
  if (solver.info() != Success)
  {
    cerr << "solving failed" << endl;
  }
  
  return x;
}

int write_parameter_domain(const HalfEdgeMesh &he_mesh, const std::vector<Eigen::Vector2d> &parameter_vert)
{
  const size_t tri_num = he_mesh.get_face_num();
  vector<MatrixXd> parameter_tri;
  for (size_t i = 0; i < tri_num; ++i)
  {
    auto tri = he_mesh.get_tri_vert_id(i);
    MatrixXd triangle(2, 3);
    for (int v = 0; v < 3; ++v)
      triangle.col(v) = parameter_vert.at(tri[v]);
    parameter_tri.push_back(triangle);
  }

  write_to_obj<2>(parameter_tri, "parameter.obj");
  return 0;
}

int write_model_with_texture(const HalfEdgeMesh &he_mesh, const std::vector<Eigen::Vector2d> &parameter_vert, const char* const path)
{
  ofstream f_out(path);
  if (!f_out)
  {
    cerr << "error in file open" << endl;
    return 1;
  }

  const size_t vert_num = he_mesh.get_vert_num();
  for (size_t i = 0; i < vert_num; ++i)
  {
    Vector3d v = he_mesh.get_vert(i);
    f_out << "v " << v[0] << " " << v[1] << " " << v[2] << endl;
  }

  const size_t pv_num = parameter_vert.size();
  for (size_t i = 0; i < pv_num; ++i)
    f_out << "vt " << parameter_vert[i][0] << " " << parameter_vert[i][1] <<endl;

  const size_t tri_num = he_mesh.get_face_num();
  for (size_t i = 0; i < tri_num; ++i)
  {
    auto tri = he_mesh.get_tri_vert_id(i);
    f_out << "f " << tri[0]+1 << "/" << tri[0]+1 << " " << tri[1]+1 << "/" << tri[1]+1 << " " << tri[2]+1 << "/" << tri[2]+1 << endl;
  }

  f_out.close();
  return 0;
}
