#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

#include <vector>
using namespace std;

int sgn(float val);

vector<vector<float>> scalar_multiply_matrix(const vector<vector<float>>& A, float scalar);
vector<vector<float>> subtract_matrix(const vector<vector<float>>& A, const vector<vector<float>>& B);
vector<float> scalar_multiply(const vector<float>& vec, float scalar);
float dot_product(const vector<float>& vec1, const vector<float>& vec2);
vector<float> add_vectors(const vector<float>& vec1, const vector<float>& vec2);
vector<float> subtract_vectors(const vector<float>& vec1, const vector<float>& vec2);

vector<float> get_column(const vector<vector<float>>& mat, int col);
void set_column(vector<vector<float>>& mat, const vector<float>& vect, int col);

vector<vector<float>> outer_product(const vector<float>& vect1, const vector<float>& vect2);
vector<vector<float>> identity_matrix(int n);
vector<vector<float>> matrix_mult(const vector<vector<float>>& mat1, const vector<vector<float>>& mat2);
void transpose_matrix(vector<vector<float>>& mat);

void print_vector(const vector<float>& vec);
void print_matrix(const vector<vector<float>>& mat);
void print_matrix_transpose(const vector<vector<float>>& mat);

#endif
