#include <iostream>
#include <vector>
#include <cmath>
#include "../utils/matrix_utils.h"

using namespace std;

// Since I only care about real life stuff, we assume all square matrix systems


int sgn(float val) {
    return (0.0f < val) - (val < 0.0f);
}

pair<vector<vector<float>>, vector<vector<float>>>
 gs_QR(vector<vector<float>>& in_matrix, const int n) {
    transpose_matrix(in_matrix);

    // Gram-Schmidt Process for QR Decomposition
    vector<vector<float>> q_matrix(n, vector<float>(n));
    vector<vector<float>> r_matrix(n, vector<float>(n, 0.0f));
    vector<float> e_vect(n);
    vector<float> u_vect(n);
    vector<float> proj_vect(n);
    float proj_int;

    for (int i = 0; i < n; i++){
        u_vect = in_matrix[i];
        e_vect = u_vect;
        for (int j = 0; j < i; j++){
            proj_int = dot_product(q_matrix[j],u_vect)/dot_product(q_matrix[j],q_matrix[j]);
            proj_vect = scalar_multiply(q_matrix[j], proj_int);
            e_vect = subtract_vectors(e_vect, proj_vect);
        }
        proj_int = 1/sqrt(dot_product(e_vect, e_vect));
        q_matrix[i] = scalar_multiply(e_vect, proj_int);

        for (int j = i; j < n; j++){
            r_matrix[j][i] = dot_product(q_matrix[i], in_matrix[j]);
        }
    }

    transpose_matrix(q_matrix);
    transpose_matrix(r_matrix);
    return {q_matrix, r_matrix};
}

pair<vector<vector<float>>, vector<vector<float>>>
 householder_QR(vector<vector<float>>& in_matrix, int n) {

    vector<float> y_k(n);
    float norm_y;
    float norm_w;
    vector<vector<float>> identity = identity_matrix(n);
    vector<vector<float>> h_k(n, vector<float>(n));
    vector<vector<float>> q_k = identity;


    // Householder QR Decomposition math
    for (int k = 0; k < n; k++) {

        int v_k_len = n - k;
        int zero_pad = k;

        y_k = get_column(in_matrix, k);

        vector<float> w_k(v_k_len);
        vector<float> v_k(n, 0.0f);

        for (int i = 0; i < v_k_len; i++) {
            w_k[i] = y_k[k + i];
        }

        norm_y = sqrt(dot_product(w_k, w_k));
        w_k[0] = w_k[0] + (sgn(w_k[0]) *norm_y);
        norm_w = sqrt(dot_product(w_k, w_k));

        w_k = scalar_multiply(w_k, 1 / norm_w);
        
        for (int i = zero_pad; i < n; i++) {
            v_k[i] = w_k[i - zero_pad];
        }

        h_k = subtract_matrix(identity,
              scalar_multiply_matrix(outer_product(v_k, v_k), 2));

        q_k = matrix_mult(q_k, h_k);

        in_matrix = matrix_mult(h_k, in_matrix);
    }

    return {q_k, in_matrix};
}
