#include <iostream>
#include <sstream>
#include <ctime>
#include <boost/random.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/algorithm/string.hpp>
#include <glog/logging.h>
#include <google/gflags.h>
#include "utils.h"
#include "graph.h"
DEFINE_string(inputMatrix, "", "input graph file");
DEFINE_string(inputLabels, "", "input label file");
using namespace std;
int rnd = 0;
bool prob(double event)
{
    using boost::mt19937;
    using boost::uniform_real;
    using boost::variate_generator;

    mt19937 gen(static_cast<unsigned long>(time(0)+rnd++));
    uniform_real<> dst(0, 1);
    typedef variate_generator<mt19937&, uniform_real<> > generator;
    generator rand( gen, dst );
    if (event > rand()) return true;
    else return false;
}
int dwalk(int node, graph::Matrix& mat)
{
    // sink
    return 1;
}
void debugPring (const graph::Matrix& mat)
{
    const int N = mat.size();
    for(int i=0; i<N; i++){
        const int E = mat[i].size();
        cout << i << ":" << endl;
        for(int j=0; j<E; j++){
            cout << mat[i][j].node << ":" << mat[i][j].weight << endl;
        }
    }
    LOG(INFO) << "mat size: " << mat.size() << endl;
}
int main (int argc, char** argv)
{
    google::ParseCommandLineFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);

    graph::Matrix trans_mat;
    graph::Matrix norm_trans_mat;
    graph::Labels labels;

    try {
        graph::load_mat(trans_mat, norm_trans_mat, FLAGS_inputMatrix);
        graph::load_lab(labels, FLAGS_inputLabels);
    } catch (const exception &e) {
        cout << "ERROR: " << e.what() << endl;
    }
    assert(trans_mat.size() == labels.size());
    const int N = trans_mat.size(); // 1 .. N
    int labeled_cnt = 0;
    vector<int> unlabeled;
    for (int i=0; i<N; i++) {
        if (labels[i] < 0) {
            unlabeled.pb(i+1);
        } else {
            labeled_cnt++;
        }
    }
    const int L = labeled_cnt;
    const int U = N-L;
    assert(U == static_cast<int>(unlabeled.size()));
    assert(L > 0);
    assert(U > 0);
    LOG(INFO) << "# of nodes: " << N << endl;
    LOG(INFO) << "# of labeled nodes: " << L << endl;
    LOG(INFO) << "# of unlabeled nodes: " << U << endl;

    graph::Matrix norm_trans_uu(N);
    graph::Matrix norm_trans_ul(N);
    load_submatrix(norm_trans_mat, norm_trans_uu, norm_trans_ul, U, L);
    show_normalized_trans(norm_trans_uu);
    show_normalized_trans(norm_trans_ul);
    const int C = 2;
    vector<vector<int> > y(N, vector<int>(C, 0));
    for (int i=0; i<N; i++) {
        assert(labels[i] != 0);
        if (labels[i] > 0) {
            // labeled
            const int label_index = labels[i] - 1;
            y[i][label_index] = 1;
        }
    }
    // N x C matrix は adjcency matrix
}

