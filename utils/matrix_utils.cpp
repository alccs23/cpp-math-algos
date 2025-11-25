#include "matrix_utils.h"
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

vector<vector<double>> scalar_multiply_matrix(const vector<vector<double>>& A, double scalar) {
    int rows = A.size();
    int cols = A[0].size();
    vector<vector<double>> result(rows, vector<double>(cols));
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result[i][j] = A[i][j] * scalar;
    return result;
}

vector<vector<double>> subtract_matrix(const vector<vector<double>>& A, const vector<vector<double>>& B) {
    int rows = A.size();
    int cols = A[0].size();
    vector<vector<double>> result(rows, vector<double>(cols));
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result[i][j] = A[i][j] - B[i][j];
    return result;
}

vector<double> scalar_multiply(const vector<double>& vec, double scalar) {
    vector<double> result(vec.size());
    for (int i = 0; i < vec.size(); i++)
        result[i] = vec[i] * scalar;
    return result;
}

double dot_product(const vector<double>& vec1, const vector<double>& vec2) {
    double dot = 0.0;
    for (int i = 0; i < vec1.size(); i++)
        dot += vec1[i] * vec2[i];
    return dot;
}

vector<double> add_vectors(const vector<double>& vec1, const vector<double>& vec2) {
    vector<double> result(vec1.size());
    for (int i = 0; i < vec1.size(); i++)
        result[i] = vec1[i] + vec2[i];
    return result;
}

vector<double> subtract_vectors(const vector<double>& vec1, const vector<double>& vec2) {
    size_t size = vec1.size();
    vector<double> result(size);
    for (size_t i = 0; i < size; i++) {
        result[i] = vec1[i] - vec2[i];
    }
    return result;
}

vector<double> get_column(const vector<vector<double>>& mat, int col) {
    vector<double> result(mat.size());
    for (int i = 0; i < mat.size(); i++)
        result[i] = mat[i][col];
    return result;
}

void set_column(vector<vector<double>>& mat, const vector<double>& vect, int col) {
    for (int i = 0; i < mat.size(); i++)
        mat[i][col] = vect[i];
}

vector<vector<double>> outer_product(const vector<double>& v1, const vector<double>& v2) {
    int n = v1.size();
    vector<vector<double>> result(n, vector<double>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            result[i][j] = v1[i] * v2[j];
    return result;
}

vector<vector<double>> identity_matrix(int n) {
    vector<vector<double>> I(n, vector<double>(n, 0));
    for (int i = 0; i < n; i++)
        I[i][i] = 1.0;
    return I;
}

vector<vector<double>> matrix_mult(const vector<vector<double>>& A, const vector<vector<double>>& B) {
    int n = A.size();
    vector<vector<double>> result(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                result[i][j] += A[i][k] * B[k][j];
    return result;
}

vector<double> matrix_vector_mult(const vector<vector<double>>& A, const vector<double>& v) {
    int n = A.size();
    vector<double> result(n, 0.0);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            result[i] += A[i][j] * v[j];

    return result;
}

void print_vector(const vector<double>& vec) {
    for (double v : vec)
        cout << v << " ";
    cout << endl;
}

void print_matrix(const vector<vector<double>>& mat) {
    const double eps = 1e-12;

    const int width = 12;
    const int precision = 12;

    for (const auto& row : mat) {
        for (double val : row) {
            if (fabs(val) < eps) val = 0.0;

            cout << setw(width) << fixed << setprecision(precision) << val;
        }
        cout << "\n";
    }
}

void print_matrix_transpose(const vector<vector<double>>& mat) {
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

void transpose_matrix(vector<vector<double>>& mat) {
    size_t n = mat.size();
    if (n == 0) return;

    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            swap(mat[i][j], mat[j][i]);
        }
    }
}

vector<double> get_diag(const vector<vector<double>>& mat){
    size_t n = mat.size();
    vector<double> result(n);
    for (int i = 0; i < n; i++){
        result[i] = mat[i][i];
    }
    return result;
}

double get_vect_max(const vector<double>& vec) {
    if (vec.empty()) {
        throw invalid_argument("Vector is empty.");
    }

    double max_val = vec[0];
    for (double v : vec) {
        if (v > max_val) max_val = v;
    }
    return max_val;
}
