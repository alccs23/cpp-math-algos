#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "../qr_decomp/QR_decomp.h"
#include "../utils/matrix_utils.h" 

pair<vector<vector<double>>, vector<double>> qr_algorithm(const vector<vector<double>>& mat) {

    size_t n = mat.size();

    vector<double> result(n);

    vector<double> residual(n, 1.0);
    vector<vector<double>> a_k = mat;
    vector<vector<double>> Q_last = identity_matrix(n);

    while (get_vect_max(residual) > 1e-8) {
        vector<double> cur_diag = get_diag(a_k);
        auto [Q, R] = givens_QR(a_k, n);
        Q_last = matrix_mult(Q_last, Q);
        a_k = matrix_mult(R, Q);
        residual = subtract_vectors(cur_diag, get_diag(a_k));
    }

    return {Q_last, get_diag(a_k)};
}

