#include "barycenter_mapping.h"
#include "half_edge.h"
#include "write_to_file.h"
#include <iostream>
#include <vector>




using namespace std;
using namespace Eigen;


Eigen::Vector2d get_boundary_parameter_coordinate(double t)
{
  if (t < 0.25)
    return {0, t / 0.25};
  else if (t < 0.5)
    return {1, (t-0.25) / 0.25};
  else if (t < 0.75)
    return {(0.75 - t) / 0.25, 1};
  else
    return {0, (1 - t) / 0.25};
}

std::vector<Eigen::Vector2d> barycenter_mapping(const HalfEdgeMesh &he_mesh, const std::vector<size_t> &loop)
{
  const auto he_e = he_mesh.he_e_;
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
  fill(parameter_vert.begin(), parameter_vert.end(), Vector2d(0, 0));
  // parameter_vert[0] = {0, 0};
  double length = 0;
  for (auto e : loop)
  {
    size_t ne = he_e.at(he_e.at(e).next_e).next_e;
    array<size_t, 2> he_v = {he_e.at(ne).he_v, he_e.at(e).he_v};
    double le = (he_mesh.get_vert(he_v[0]) - he_mesh.get_vert(he_v[1])).norm();
    length += le;
    double t = length / loop_length;
    parameter_vert[he_v[1]] = get_boundary_parameter_coordinate(t);
  }

  
  return parameter_vert;
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

  write_to_vtk<2>(parameter_tri, "parameter.vtk");
  return 0;
}
