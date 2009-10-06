#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>
#include <glog/logging.h>
#include "utils.h"

namespace graph {
  typedef struct Edge {
    int node;
    double weight;
    Edge (int node, double weight)
    : node(node), weight(weight) {}
  } Edge;
  typedef std::vector<Edge> Array;
  typedef std::vector<Array> Matrix;
  typedef std::vector<int> Labels;
  typedef std::vector<std::vector<double> > LabelMatrix;

  void normalize (Matrix& trans_mat, Matrix& norm_trans_mat);
  void row_normalize (Matrix& mat);
  void load_mat(Matrix& trans_mat, Matrix& norm_trans_mat, const std::string& input);
  int  load_lab(Labels& lab, const std::string& intpu);
  void load_submatrix(const Matrix& mat, Matrix& mat_uu, Matrix& mat_ul,
                      const int U, const int L);
  void show_normalized_trans(const Matrix& norm);
  void show_normalized_trans_u(const Matrix& norm, const int L);
}

#endif
