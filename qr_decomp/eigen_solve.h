#ifndef EIGEN_SOLVE_H
#define EIGEN_SOLVE_H

#include <vector>
using namespace std;

pair<vector<vector<double>>, vector<double>> qr_algorithm(const vector<vector<double>>& mat);

#endif
