#include "lscm.h"
#include "half_edge.h"
#include "is_model_patch.h"
#include "barycenter_mapping.h"
#include <iostream>
#include <random>

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
  const int vert_num = he_mesh.get_vert_num();
  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<int> vert_u(0, vert_num-1);
  array<int, 2> fix_v = {vert_u(gen), vert_u(gen)};
  while (fix_v[1] == fix_v[0])
    fix_v[1] = vert_u(gen);
  
  cerr << he_mesh.get_vert_num() <<" " << fix_v[0] << " " << fix_v[1] << endl;
  vector<double> tri_area;
  vector<MatrixXd> tri_transform;
  tie(tri_area, tri_transform) = get_triangle_area_and_transform(he_mesh);
  cerr << tri_area.size() << " " << tri_transform.size() << endl;
  vector<Vector2d> para = solver(fix_v[0], fix_v[1], tri_area, tri_transform, he_mesh);

  scale_parameter_to_unit(para);
  write_parameter_domain(he_mesh, para);
  write_model_with_texture(he_mesh, para, "patch.obj");
  cerr << "Hello World!" << endl;
  return 0;
}

