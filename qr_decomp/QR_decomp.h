#ifndef QR_DECOMP_H
#define QR_DECOMP_H

#include <vector>

using namespace std;

pair<vector<vector<float>>, vector<vector<float>>> gs_QR(vector<vector<float>>& in_matrix, const int n);

pair<vector<vector<float>>, vector<vector<float>>> householder_QR(vector<vector<float>>& in_matrix, int n);

#endif