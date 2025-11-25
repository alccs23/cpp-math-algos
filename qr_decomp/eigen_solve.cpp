#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "QR_decomp.h"
#include "../utils/matrix_utils.h" 

vector<float> qr_algorithm(const vector<vector<float>>& mat){

    size_t n = mat.size();

    //list of eigenvalues
    vector<float> result(n);

    vector<float> residual(n, 1.0f);
    vector<vector<float>> a_k = mat;

    //1e-5 precision is good enough, dont have small eigenvalues please!!!!
    while (get_vect_max(residual) > 1e-5){
        vector<float> cur_diag = get_diag(a_k);
        auto [Q, R] = givens_QR(a_k, n);
        a_k = matrix_mult(R, Q);
        residual = subtract_vectors(cur_diag, get_diag(a_k));
    }

    return get_diag(a_k);
}