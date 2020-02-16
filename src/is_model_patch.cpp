#include "is_model_patch.h"
#include "hash_key.h"
#include <iostream>
#include <vector>
#include <string>
#include <list>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <Eigen/Core>


using namespace std;
using namespace Eigen;



bool is_model_patch(const char* const path)
{
  ifstream f_in(path);
  if (!f_in)
  {
    cout << "error in file open!" << endl;
    return false;
  }

  unordered_map<Vector3d, int, HashFunc<Vector3d, 3>, EqualKey<Vector3d, 3>>
    vert_map;
  vector<Vector3d> mesh_vert;
  vector<Vector3i> mesh_tri;
  unordered_map<int, int> old_id_to_new;
  string str_line;
  while (getline(f_in, str_line))
  {
    if (str_line.size() < 2)
      continue;
    if (str_line.substr(0, 2) == "v ")
    {
      Vector3d vert;
      sscanf(str_line.c_str(), "%*s%lf%lf%lf", &vert[0], &vert[1], &vert[2]);
      auto ret = vert_map.emplace(vert, mesh_vert.size());
      old_id_to_new.emplace(old_id_to_new.size()+1, ret.first->second);
      if (ret.second == true)
        mesh_vert.push_back(vert);
    }
    else if (str_line.substr(0, 2) == "f ")
    {
      Vector3i tri;
      if (str_line.find("/") == string::npos)
        sscanf(str_line.c_str(), "%*s%d%d%d", &tri[0], &tri[1], &tri[2]);
      else
        sscanf(str_line.c_str(), "%*s%d%*s%d%*s%d", &tri[0], &tri[1], &tri[2]);
      tri = {old_id_to_new.at(tri[0]),
             old_id_to_new.at(tri[1]),
             old_id_to_new.at(tri[2])};
      mesh_tri.push_back(tri);
    }
  }

  // exclude three half edge on the same edge
  unordered_set<Vector2i, HashFunc<Vector2i, 2>, EqualKey<Vector2i, 2>> edge_set;
  for (const auto &tri : mesh_tri)
  {
    for (int i = 0; i < 3; ++i)
    {
      Vector2i e = {tri[i], tri[(i+1) % 3]};
      auto ret = edge_set.emplace(e);
      if (ret.second == false)
        return false;
    }
  }

  list<Vector2i> boundary_edge;
  for (const auto &e : edge_set)
  {
    Vector2i ce = {e[1], e[0]};
    if (edge_set.count(ce))
      continue;
    boundary_edge.push_back(e);
  }
  if (boundary_edge.empty())
    return false;

  list<Vector2i> loop;
  loop.push_back(boundary_edge.front());
  boundary_edge.pop_front();
  while (loop.front()[0] != loop.back()[1])
  {
    bool is_edge_add = false;
    for (auto itr = boundary_edge.begin(); itr != boundary_edge.end();)
    {
      if ((*itr)[0] == loop.back()[1])
      {
        loop.push_back(*itr);
        itr = boundary_edge.erase(itr);
        is_edge_add = true;
        continue;
      } 
      if ((*itr)[1] == loop.front()[0])
      {
        loop.push_front(*itr);
        itr = boundary_edge.erase(itr);
        is_edge_add = true;
        continue;
      }
      itr++;
    }
    if (is_edge_add == false)
      return false;
  }

  return true;
}

