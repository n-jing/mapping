#ifndef FWRITE_TO_FILE_JJ_H
#define FWRITE_TO_FILE_JJ_H


#include <vector>
#include <Eigen/Core>
#include <map>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>

// #define REMOVE_DUPLICATE_VERT_IN_VTK
// #define TWO_COLOR

//*****************************************************************************
//! \brief vert comparision struct used for map and set
//*****************************************************************************
template<typename Derived>
struct EigenVertComp
{
  bool operator()(const Eigen::MatrixBase<Derived> &v1,
                  const Eigen::MatrixBase<Derived> &v2) const
    {
      for (int i = 0; i < v1.cols(); ++i)
      {
        if (v1(i) < v2(i)) return true;
        if (v1(i) > v2(i)) return false;
      }
      return false;
    }
};

//*****************************************************************************
//
//! \brief write geometry to vtk, with template parameter dimension, element and container
//!
//! \param D dimension
//! \param Derived Eigen element tyep
//! \param C container
//! \note stl containers have two template parameters, type and allocate
//
//*****************************************************************************

template<size_t D = 3,
         typename Derived = Eigen::MatrixXd,
         template<typename...> class C = std::vector>
int write_to_vtk(
  const C<Derived> &table_element,
  const char *path = "demo.vtk", const double *ptr_table_color = nullptr)
{
  using T = typename Derived::Scalar;
  using VectorXst = Eigen::Matrix<size_t, Eigen::Dynamic, 1>;
  std::ofstream f_out(path);
  if (!f_out)
  {
    std::cerr << "error:: file open" << std::endl;
    return -1;
  }
  
  std::map<Eigen::Matrix<T, D, 1>, size_t,
           EigenVertComp<Eigen::Matrix<T, D, 1>>> map_vert;
  std::vector<Eigen::Matrix<T, D, 1>> table_vert;
  std::vector<VectorXst> table_face_id;
  size_t cell_sum = 0;
  for (const auto &itr : table_element)
  {
    if (itr.rows() != 3)
      std::cerr << "[ warning ]: element dimension is " <<itr.rows() << std::endl;

    size_t num = itr.cols();
    VectorXst face(num);
    for (size_t id_v = 0; id_v < num; ++id_v)
    {
      const auto &v = itr.col(id_v);
      auto is_insert = map_vert.emplace(itr.col(id_v), map_vert.size());
      bool is_s = is_insert.second;
#ifndef REMOVE_DUPLICATE_VERT_IN_VTK
      is_s = true;
#endif
      if (is_s == true)
      {
        table_vert.push_back(itr.col(id_v));
        face(id_v) = table_vert.size() - 1;
      } 
      else
        face(id_v) = map_vert.at(v);
    }
    table_face_id.push_back(face);
    cell_sum += (num + 1);
  }
  
  f_out << "# vtk DataFile Version 2.8\nUnstructured Grid Example\nASCII\nDATASET UNSTRUCTURED_GRID\n";
  f_out << "POINTS " << table_vert.size() << " double\n";
  for (const auto &v : table_vert)
  {
    for (size_t i = 0; i < D; ++i)
      f_out << v(i) << " ";
    for (size_t i = D; i < 3; ++i)
      f_out << " 0";

    f_out << "\n";
  } 


  f_out << "CELLS " << table_face_id.size() << " " << cell_sum << "\n";
  for (const auto &itr : table_face_id)
  {
    f_out << itr.size();
    size_t num = itr.size();
    for (size_t id_v = 0; id_v < num; ++id_v)
      f_out << " " << itr(id_v);
    f_out << "\n";
  }

  f_out << "CELL_TYPES " << table_face_id.size() << "\n";
  for (const auto &itr : table_face_id)
  {
    size_t num = itr.size();
    switch (num)
    {
    case 2:
      f_out << "3\n";
      break;
    case 1:
      f_out << "1\n";
      break;
    case 3:
      f_out << "5\n";
      break;
    default:
      f_out << "7\n";
    }
  }

  f_out << "CELL_DATA " << table_face_id.size() << "\n";
  f_out << "SCALARS edge double 1\n";
  f_out << "LOOKUP_TABLE my_table\n";

  size_t num_face = table_face_id.size();
  for (size_t itr = 0; itr < num_face; ++itr)
  {
    if (ptr_table_color != nullptr)
      f_out << ptr_table_color[itr] << "\n";
    else
#ifndef TWO_COLOR
      f_out << fmod(itr * 0.1, 1.0999999999999999) << "\n";
#else
    f_out << itr % 2 << "\n";
#endif
  }

#ifndef TWO_COLOR
  f_out << "LOOKUP_TABLE my_table 11" << "\n";

  f_out << "0.890     0.090    0.051    1" << "\n"; // red
  f_out << "0.4       1        0.3      1" << "\n"; // blue
  f_out << "0.26666   0.8      0        1" << "\n"; // blue
  f_out << "0         0.302    0.902    1" << "\n"; // blue
  f_out << "0.902     0.75     0        1" << "\n"; // blue
  f_out << "0.98      1        0.94     1" << "\n"; // white
  f_out << "226       158      186      1" << "\n";
  f_out << "141       230      104      1" << "\n";
  f_out << "85        173      168      1" << "\n";
  f_out << "44        230      185      1" << "\n";
  f_out << "165       85       173      1" << "\n";
#else
  f_out << "LOOKUP_TABLE my_table 2\n";

  f_out << "0.890     0.090    0.051   1" << "\n"; // red
  f_out << "0.4       1        0.3      1" << "\n"; // blue
#endif

  f_out.close();
  return 0;
}

template<size_t D = 3,
         typename Derived = Eigen::MatrixXd,
         template<typename...> class container = std::vector>
int write_to_obj(
  const container<Derived> &table_element,
  const char *const path = "demo.obj")
{
  std::ofstream f_out(path);
  if (!f_out)
  {
    std::cerr << "error:: file open" << std::endl;
    return -1;
  }
  
  std::stringstream s_str_tri;

  size_t v_count = 1;
  for (auto &p : table_element)
  {
    switch(p.cols())
    {
    case 3:
      s_str_tri << "f";
      break;
    case 2:
      s_str_tri << "l";
      break;
    }
    for (int i = 0; i < p.cols(); ++i)
    {
      f_out << "v";
      for (size_t r = 0; r < D; ++r)
        f_out << " " << p(r, i);
      for (size_t r = D; r < 3; ++r)
        f_out << " 0";
      f_out << "\n";
      s_str_tri << " " << v_count;
      ++v_count;
    }
    s_str_tri << "\n";
  }

  f_out << s_str_tri.str();
  f_out.close();
  return 0;
}

#endif // FWRITE_TO_FILE_JJ_H
