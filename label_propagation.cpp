#include <iostream>
#include <iomanip>
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
DEFINE_int32(iteration, 1000, "iteration maximum");
DEFINE_double(eps, 1.0e-9, "error precision");
DEFINE_bool(showWeight, false, "show label weight matrix");
DEFINE_int32(precision, 10, "precision");
using namespace std;

void showLabels(const graph::LabelMatrix& y_l, const graph::LabelMatrix& y_u,
                int L, int U, int C,
                const vector<int>& labeled_nodes,
                const vector<int>& unlabeled_nodes,
                const graph::Labels& lab)
{
    const int N = L + U;
    cout.setf(std::ios::fixed, std::ios::floatfield);
    cout.precision(FLAGS_precision);
    if (FLAGS_showWeight) {
        for (int i=0; i<N; i++) {
            if (lab[i] >= 0) {
                // labeled
                cout << "L: ";
                for (int c=0; c<C; c++) {
                    if (c!=0) cout << " ";
                    cout << y_l[c][i];
                }
                cout << endl;
            } else {
                // unlabeled
                cout << "U: ";
                for (int c=0; c<C; c++) {
                    if (c!=0) cout << " ";
                    cout << y_u[c][i];
                }
                cout << endl;
            }
        }
    } else {
        for (int i=0; i<N; i++) {
            int argmax = 0;
            double argmax_val = -1.0;

            for (int c=0; c<C; c++) {
                // argmax
                if (lab[i] >= 0) {
                    for (int c=0; c<C; c++) {
                        if (argmax_val < y_l[c][i]) {
                            argmax_val = y_l[c][i];
                            argmax = c;
                        }
                    }
                } else {
                    for (int c=0; c<C; c++) {
                        if (argmax_val < y_u[c][i]) {
                            argmax_val = y_u[c][i];
                            argmax = c;
                        }
                    }
                }
            }

            cout << argmax+1 << endl;
        }
    }
}
int main (int argc, char** argv)
{
    google::ParseCommandLineFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);

    graph::Matrix trans_mat;
    graph::Matrix norm_trans_mat;
    graph::Labels labels;

    int max_lab;
    try {
        graph::load_mat(trans_mat, norm_trans_mat, FLAGS_inputMatrix);
        max_lab =  graph::load_lab(labels, FLAGS_inputLabels);
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
    vector<int> labeled_nodes, unlabeled_nodes; // collection of node index
    for (int i=0; i<N; i++) {
        if (labels[i] < 0) {
            unlabeled_nodes.pb(i);
        } else {
            labeled_nodes.pb(i);
        }
    }
    assert(static_cast<int>(unlabeled_nodes.size()) == U);
    assert(static_cast<int>(labeled_nodes.size()) == L);
    assert(U == static_cast<int>(unlabeled.size()));
    assert(L > 0);
    assert(U > 0);
    LOG(INFO) << "# of nodes: " << N << endl;
    LOG(INFO) << "# of labeled nodes: " << L << endl;
    LOG(INFO) << "# of unlabeled nodes: " << U << endl;

    graph::Matrix norm_trans_uu(N, graph::Array());
    graph::Matrix norm_trans_ul(N, graph::Array()); // start from L
    load_submatrix(norm_trans_mat, norm_trans_uu, norm_trans_ul,
                   U, L, unlabeled_nodes, labeled_nodes,
                   labels);
/*
    LOG(INFO) << "- show_normalized_trans: UU";
    show_normalized_trans_u(norm_trans_uu, 0);
    LOG(INFO) << "- show_normalized_trans: UL";
    show_normalized_trans_u(norm_trans_ul, 0);
    LOG(INFO) << "- show_normalized_trans: MAT";
    show_normalized_trans(norm_trans_mat);
*/
    const int C = max_lab;
    graph::LabelMatrix y_u(C, vector<double>(N, 0.0));
    graph::LabelMatrix y_l(C, vector<double>(N, 0.0));
    foreach (const int i, labeled_nodes) {
        assert(labels[i] != 0);
        if (labels[i] > 0) {
            // labeled
            const int label_index = labels[i] - 1;
            y_l[label_index][i] = 1.0;
        }
    }
//    showLabels(y_l, y_u, L, U, C, labeled_nodes, unlabeled_nodes, labels);
    const int max_iteration = FLAGS_iteration;
    int iteration;
    for (iteration = 0; iteration<max_iteration; iteration++) {
        vector<vector<double> > y_ret(C, vector<double>(N, 0.0)); // y_u'
        // calc y_u' = T_ul y_l
        // calc y_u' = T_uu y_u
        double err = 0.0;
        foreach (const int i, unlabeled_nodes) {
            const int edge_sz = norm_trans_ul[i].size();
            for (int c=0; c<C; c++) {
                for (int j=0; j<edge_sz; j++) {
                    const double w = norm_trans_ul[i][j].weight;
                    const int src_index = norm_trans_ul[i][j].node - 1;
                    const double y_u_w = y_l[c][src_index];
                    if (y_u_w < 1.0e-200 || w < 1.0e-200) continue;
                    y_ret[c][i] += w * y_l[c][src_index];
                }
            }
        }

        foreach (const int i, unlabeled_nodes) {
            const int edge_sz = norm_trans_uu[i].size();
            for (int c=0; c<C; c++) {
                for (int j=0; j<edge_sz; j++) {
                    const double w = norm_trans_uu[i][j].weight;
                    const int src_index = norm_trans_uu[i][j].node - 1;
                    const double y_u_w = y_u[c][src_index];
                    if (y_u_w < 1.0e-200 || w < 1.0e-200) continue;
                    y_ret[c][i] += w * y_u[c][src_index];
                }
                err += std::abs(y_ret[c][i] - y_u[c][i]);
            }
        }
        y_u = y_ret;
//        LOG(INFO) << "iteration: " << iteration << " done. error: " << err << endl;
        if (FLAGS_eps > err) break;
    }
    LOG(INFO) << "iteration: " << iteration << " done." << endl;
    showLabels(y_l, y_u, L, U, C, labeled_nodes, unlabeled_nodes, labels);

}
