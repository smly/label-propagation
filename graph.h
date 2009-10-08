#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <fstream>
#include <boost/algorithm/string.hpp>
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
  void load_mat(Matrix& trans_mat, Matrix& norm_trans_mat, const std::string& input);
  int  load_lab(Labels& lab, const std::string& intpu);
  void load_submatrix(const Matrix& mat, Matrix& mat_uu, Matrix& mat_ul,
                      const int U, const int L,
                      const std::vector<int> unlabeled_nodes,
                      const std::vector<int> labeled_nodes,
                      const Labels& lab);
}

#endif
