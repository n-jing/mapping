#ifndef SAVE_MATRIX_AS_IMAGE_JJ_H
#define SAVE_MATRIX_AS_IMAGE_JJ_H


#include <string>
#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>
#include <Eigen/SparseCore>
#include <unsupported/Eigen/SparseExtra>



template<class Derived>
void save_matrix_as_image(const Derived &mat, const char *const file = "image.png")
{
  Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic> M = mat.cwiseAbs();
  typename Derived::Scalar max_coeff = M.maxCoeff();
  M = (255.0 / max_coeff) * M;
  cv::Mat image;
  cv::eigen2cv(M, image);
  cv::imwrite(file, image);
}

template<typename DerivedA, typename DerivedB>
void save_matrix_as_mtx(const DerivedA &A, const DerivedB &b, const char *const path = "solver_test")
{
  bool ret = std::is_same<Eigen::SparseMatrix<typename DerivedA::Scalar>, DerivedA>::value;
  if (ret)
  {
    Eigen::Matrix<typename DerivedA::Scalar, Eigen::Dynamic, Eigen::Dynamic> a(A);
    Eigen::saveMarket(a, std::string(std::string(path) + ".mtx").c_str());
  }
  else
    Eigen::saveMarket(A, std::string(std::string(path) + ".mtx").c_str());

  Eigen::saveMarketVector(b, std::string(std::string(path) + "_b.mtx").c_str());
}

template<typename DerivedA>
void save_matrix_as_mtx(const DerivedA &A, const char *const path = "solver_test")
{
  bool ret = std::is_same<Eigen::SparseMatrix<typename DerivedA::Scalar>, DerivedA>::value;
  if (ret)
  {
    Eigen::Matrix<typename DerivedA::Scalar, Eigen::Dynamic, Eigen::Dynamic> a(A);
    Eigen::saveMarket(a, std::string(std::string(path) + ".mtx").c_str());
  }
  else
    Eigen::saveMarket(A, std::string(std::string(path) + ".mtx").c_str());
}

#endif // SAVE_MATRIX_AS_IMAGE_JJ_H
