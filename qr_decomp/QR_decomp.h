#ifndef QR_DECOMP_H
#define QR_DECOMP_H

#include <vector>
using namespace std;

pair<vector<vector<double>>, vector<vector<double>>> gs_QR(vector<vector<double>>& in_matrix, const int n);

pair<vector<vector<double>>, vector<vector<double>>> householder_QR(vector<vector<double>>& in_matrix, int n);

pair<vector<vector<double>>, vector<vector<double>>> givens_QR(vector<vector<double>>& in_matrix, int n);

#endif
