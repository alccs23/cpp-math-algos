#include <iostream>
#include <vector>
#include <cmath>
#include "../utils/matrix_utils.h"

using namespace std;

// Since I only care about real life stuff, we assume all square matrix systems


void gs_QR(vector<vector<float>>& in_matrix, const int n) {
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
    cout << "" << endl;
    cout << "R Matrix is: " << endl;
    print_matrix_transpose(r_matrix);
    cout << "" << endl;
    cout << "Q Matrix is: " << endl;
    print_matrix_transpose(q_matrix);
    cout << "" << endl;
    cout << "QR Matrix validation: " << endl;
    // It's reversed cuz they are actually transposed in this code (WE HATE COLUMN VECTS IN CPP!!!)
    print_matrix_transpose(matrix_mult(r_matrix, q_matrix));
}
