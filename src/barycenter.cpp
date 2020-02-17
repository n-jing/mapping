#include "is_model_patch.h"
#include "half_edge.h"
#include "write_to_file.h"
#include "barycenter_mapping.h"
#include <iostream>
#include <Eigen/Core>


using namespace std;
using namespace Eigen;


int main (int argc, char *argv[])
{
  bool is_patch = is_model_patch(argv[1]);
  if (is_patch == false)
  {
    cerr << "input model is not a patch" << endl;
    return 1;
  }

  HalfEdgeMesh he_mesh(argv[1]);
  const size_t edge_num = he_mesh.get_edge_num();
  size_t boundary_he = numeric_limits<size_t>::max();
  for (size_t i = 0; i < edge_num; ++i)
    if (he_mesh.is_boundary_half_edge(i))
    {
      boundary_he = i;
      break;
    }

  const auto he_e = he_mesh.get_half_edge();
  vector<size_t> loop;
  loop.push_back(boundary_he);
  size_t next_he = he_e.at(boundary_he).next_e;
  while (he_e.at(next_he).pair_e != numeric_limits<size_t>::max())
    next_he = he_e.at(he_e.at(next_he).pair_e).next_e;
  while (next_he != boundary_he)
  {
    loop.push_back(next_he);
    next_he = he_e.at(next_he).next_e;
    while (he_e.at(next_he).pair_e != numeric_limits<size_t>::max())
      next_he = he_e.at(he_e.at(next_he).pair_e).next_e;
  }

  auto parameter_vert = barycenter_mapping(he_mesh, loop);
  write_model_with_texture(he_mesh, parameter_vert, "patch.obj");

  return 0;
}
