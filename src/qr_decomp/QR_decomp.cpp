#include <iostream>
#include <vector>
#include <cmath>
#include "../utils/matrix_utils.h"

using namespace std;

// Since I only care about real life stuff, we assume all square matrix systems

int sgn(double val) {
    return (0.0 < val) - (val < 0.0);
}

pair<vector<vector<double>>, vector<vector<double>>>
gs_QR(vector<vector<double>>& in_matrix, const int n) {
    transpose_matrix(in_matrix);

    vector<vector<double>> q_matrix(n, vector<double>(n));
    vector<vector<double>> r_matrix(n, vector<double>(n, 0.0));
    vector<double> e_vect(n);
    vector<double> u_vect(n);
    vector<double> proj_vect(n);
    double proj_int;

    for (int i = 0; i < n; i++){
        u_vect = in_matrix[i];
        e_vect = u_vect;
        for (int j = 0; j < i; j++){
            proj_int = dot_product(q_matrix[j], u_vect)/dot_product(q_matrix[j], q_matrix[j]);
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

pair<vector<vector<double>>, vector<vector<double>>>
householder_QR(vector<vector<double>>& in_matrix, int n) {

    vector<double> y_k(n);
    double norm_y;
    double norm_w;
    vector<vector<double>> identity = identity_matrix(n);
    vector<vector<double>> h_k(n, vector<double>(n));
    vector<vector<double>> q_k = identity;

    for (int k = 0; k < n; k++) {

        int v_k_len = n - k;
        int zero_pad = k;

        y_k = get_column(in_matrix, k);

        vector<double> w_k(v_k_len);
        vector<double> v_k(n, 0.0);

        for (int i = 0; i < v_k_len; i++) {
            w_k[i] = y_k[k + i];
        }

        norm_y = sqrt(dot_product(w_k, w_k));
        w_k[0] = w_k[0] + (sgn(w_k[0]) * norm_y);
        norm_w = sqrt(dot_product(w_k, w_k));

        w_k = scalar_multiply(w_k, 1 / norm_w);
        
        for (int i = zero_pad; i < n; i++) {
            v_k[i] = w_k[i - zero_pad];
        }

        h_k = subtract_matrix(identity, scalar_multiply_matrix(outer_product(v_k, v_k), 2));

        q_k = matrix_mult(q_k, h_k);

        in_matrix = matrix_mult(h_k, in_matrix);
    }

    return {q_k, in_matrix};
}

pair<vector<vector<double>>, vector<vector<double>>>
givens_QR(vector<vector<double>>& in_matrix, int n) {

    vector<vector<double>> identity = identity_matrix(n);
    vector<vector<double>> q_k = identity;
    vector<vector<double>> g_k = identity;

    for (int i = 0; i < n; i++){
        for (int j = i + 1; j < n; j++){
            if (in_matrix[j][i] != 0){
                double r = sqrt(pow(in_matrix[i][i], 2.0) + pow(in_matrix[j][i], 2.0));
                double c = in_matrix[i][i]/r;
                double s = -in_matrix[j][i]/r;
                g_k[i][i] = c;
                g_k[j][i] = s;
                g_k[i][j] = -s;
                g_k[j][j] = c;
                q_k = matrix_mult(g_k, q_k);
                in_matrix = matrix_mult(g_k, in_matrix);
                g_k = identity;
            }
        }
    }
    transpose_matrix(q_k);

    return {q_k, in_matrix};
}
