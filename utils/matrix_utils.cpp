#include "matrix_utils.h"
#include <cmath>
#include <iostream>

using namespace std;

vector<vector<float>> scalar_multiply_matrix(const vector<vector<float>>& A, float scalar) {
    int rows = A.size();
    int cols = A[0].size();
    vector<vector<float>> result(rows, vector<float>(cols));
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result[i][j] = A[i][j] * scalar;
    return result;
}

vector<vector<float>> subtract_matrix(const vector<vector<float>>& A, const vector<vector<float>>& B) {
    int rows = A.size();
    int cols = A[0].size();
    vector<vector<float>> result(rows, vector<float>(cols));
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result[i][j] = A[i][j] - B[i][j];
    return result;
}

vector<float> scalar_multiply(const vector<float>& vec, float scalar) {
    vector<float> result(vec.size());
    for (int i = 0; i < vec.size(); i++)
        result[i] = vec[i] * scalar;
    return result;
}

float dot_product(const vector<float>& vec1, const vector<float>& vec2) {
    float dot = 0.0f;
    for (int i = 0; i < vec1.size(); i++)
        dot += vec1[i] * vec2[i];
    return dot;
}

vector<float> add_vectors(const vector<float>& vec1, const vector<float>& vec2) {
    vector<float> result(vec1.size());
    for (int i = 0; i < vec1.size(); i++)
        result[i] = vec1[i] + vec2[i];
    return result;
}

vector<float> subtract_vectors(const vector<float>& vec1, const vector<float>& vec2) {
    size_t size = vec1.size();
    vector<float> result(size);
    for (size_t i = 0; i < size; i++) {
        result[i] = vec1[i] - vec2[i];
    }
    return result;
}

vector<float> get_column(const vector<vector<float>>& mat, int col) {
    vector<float> result(mat.size());
    for (int i = 0; i < mat.size(); i++)
        result[i] = mat[i][col];
    return result;
}

void set_column(vector<vector<float>>& mat, const vector<float>& vect, int col) {
    for (int i = 0; i < mat.size(); i++)
        mat[i][col] = vect[i];
}

vector<vector<float>> outer_product(const vector<float>& v1, const vector<float>& v2) {
    int n = v1.size();
    vector<vector<float>> result(n, vector<float>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            result[i][j] = v1[i] * v2[j];
    return result;
}

vector<vector<float>> identity_matrix(int n) {
    vector<vector<float>> I(n, vector<float>(n, 0));
    for (int i = 0; i < n; i++)
        I[i][i] = 1;
    return I;
}

vector<vector<float>> matrix_mult(const vector<vector<float>>& A, const vector<vector<float>>& B) {
    int n = A.size();
    vector<vector<float>> result(n, vector<float>(n, 0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                result[i][j] += A[i][k] * B[k][j];
    return result;
}

void print_vector(const vector<float>& vec) {
    for (float v : vec)
        cout << v << " ";
    cout << endl;
}

void print_matrix(const vector<vector<float>>& mat) {
    const float eps = 1e-6f;
    for (const auto& row : mat) {
        for (float val : row) {
            if (abs(val) < eps) val = 0.0f;
            cout << val << " ";
        }
        cout << endl;
    }
}

void print_matrix_transpose(const vector<vector<float>>& mat) {
    if (mat.empty()) return;

    size_t n_rows = mat.size();
    size_t n_cols = mat[0].size();

    for (size_t i = 0; i < n_cols; i++) {
        for (size_t j = 0; j < n_rows; j++) {
            cout << mat[j][i] << " ";
        }
        cout << endl;
    }
}

void transpose_matrix(vector<vector<float>>& mat) {
    size_t n = mat.size();
    if (n == 0) return;

    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            swap(mat[i][j], mat[j][i]);
        }
    }
}


