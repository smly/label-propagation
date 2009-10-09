#include "graph.h"

namespace ssl_lprop {

void normalize(Matrix& trans_mat,
               Matrix& norm_trans_mat) {
  const int V = trans_mat.size();
  std::vector<double> in_weight(V, 0.0);
  for (int i = 0; i < V; i++) {
    const int edges_sz = trans_mat[i].size();
    double edge_weight_sum = 0.0;
    for (int j = 0; j < edges_sz; j++) {
      edge_weight_sum += trans_mat[i][j].weight;
    }
    for (int j = 0; j < edges_sz; j++) {
      trans_mat[i][j].weight /= edge_weight_sum;
      // calc in_weight for row-normalize
      const int in_node = trans_mat[i][j].node;
      assert(in_node > 0);
      in_weight[in_node - 1] += trans_mat[i][j].weight;
    }
  }
  norm_trans_mat = Matrix(V);
  // calc row-normalized transition matrix
  for (int i = 0; i < V; i++) {
    const int edges_sz = trans_mat[i].size();
    for (int j=0; j<edges_sz; j++) {
      const int dst_index = trans_mat[i][j].node - 1;
      assert(dst_index >= 0);
      const double w = trans_mat[i][j].weight / in_weight[dst_index];
      norm_trans_mat[dst_index].push_back(Edge(i+1, w));
    }
  }
}
int load_lab(Labels& lab,
             const std::string& input) {
  //    if (! fs::exists(input)) {
  //        throw std::runtime_error("input file not found.");
  //    }
  std::ifstream ifs(input.c_str(), std::ios::in);
  lab = Labels();
  int max_label = 0;
  while (! ifs.eof()) {
    std::string line;
    std::getline(ifs, line);
    if (line == "?") {
      lab.push_back(-1);
    } else {
      const char* s = line.c_str();
      const int id = atoi(s);
      max_label = std::max(id, max_label);
      lab.push_back(id);
    }
    ifs.peek();
  }

  return max_label;
}
void load_mat (Matrix& trans_mat,
               Matrix& norm_trans_mat,
               const std::string& input)
{
  //    if (! fs::exists(input)) {
  //        throw std::runtime_error("input file not found.");
  //    }
  try {
    std::ifstream ifs(input.c_str(), std::ios::in);
    trans_mat = Matrix();
    while (! ifs.eof()) {
      std::string line;
      std::getline(ifs, line);
      const char *s = line.c_str();
      unsigned int elemnum = 0;
      unsigned int len = strlen(s);
      unsigned int pos = 0;
      for (unsigned int i = 0; i < len; i++) {
        if (s[i] == ':') elemnum++;
      }
      Array arry(elemnum, Edge());
      for (int i = 0; i < elemnum; i++) {
        int dst;
        double w;
        while (pos < len && isspace(s[pos])) pos++; // skip
        dst = atoi(s + pos);
        while (pos + 1 < len && s[pos] != ':') pos++;
        w = atof(s + pos + 1);
        while (pos < len && !isspace(s[pos])) pos++;
        arry[i] = Edge(dst, w);
      }
      trans_mat.push_back(arry);
      ifs.peek();
    }
  } catch (const std::exception& e) {
    std::stringstream ss;
    ss << "load_mat() says: " << e.what();
    throw std::runtime_error (ss.str());
  }
  // normalize weights
  normalize(trans_mat, norm_trans_mat);
}
void load_submatrix(const Matrix& mat,
                    Matrix& mat_uu,
                    Matrix& mat_ul,
                    const int U, const int L,
                    const std::vector<int> unlabeled_nodes,
                    const std::vector<int> labeled_nodes,
                    const Labels& lab)
{
  for (unsigned int i=0; i<unlabeled_nodes.size(); i++) {
    const int node_index = unlabeled_nodes[i];
    const int edges_sz = mat[node_index].size();
    for (int j=0; j<edges_sz; j++) {
      const int src_index = mat[node_index][j].node - 1;
      if (lab[src_index] >= 0) {
        mat_ul[i].push_back( mat[node_index][j] );
      } else {
        mat_uu[i].push_back( mat[node_index][j] );
      }
    }
  }
}

}// end of namespace ssl_lprop
