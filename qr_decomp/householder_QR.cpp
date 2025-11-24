#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "../utils/matrix_utils.h"
using namespace std;


// Since I only care about real life stuff, we assume all square matrix systems

int sgn(float val) {
    return (0.0f < val) - (val < 0.0f);
}


void householder_QR(vector<vector<float>>& in_matrix, int n) {

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

    cout << "Q Matrix:" << endl;
    print_matrix(q_k);
    cout << " " << endl;
    cout << "R Matrix:" << endl;
    print_matrix(in_matrix);
    cout << " " << endl;
    cout << "QR Matrix Validation:" << endl;
    print_matrix(matrix_mult(q_k,in_matrix));
}