#ifndef MATRIX_UTILS_H
#define MATRIX_UTILS_H

#include <vector>
using namespace std;

int sgn(double val);

vector<vector<double>> scalar_multiply_matrix(const vector<vector<double>>& A, double scalar);
vector<vector<double>> subtract_matrix(const vector<vector<double>>& A, const vector<vector<double>>& B);
vector<double> scalar_multiply(const vector<double>& vec, double scalar);
double dot_product(const vector<double>& vec1, const vector<double>& vec2);
vector<double> add_vectors(const vector<double>& vec1, const vector<double>& vec2);
vector<double> subtract_vectors(const vector<double>& vec1, const vector<double>& vec2);
double get_vect_max(const vector<double>& vec);

vector<double> get_column(const vector<vector<double>>& mat, int col);
void set_column(vector<vector<double>>& mat, const vector<double>& vect, int col);

vector<vector<double>> outer_product(const vector<double>& vect1, const vector<double>& vect2);
vector<vector<double>> identity_matrix(int n);
vector<vector<double>> matrix_mult(const vector<vector<double>>& mat1, const vector<vector<double>>& mat2);
vector<double> matrix_vector_mult(const vector<vector<double>>& A, const vector<double>& v);
void transpose_matrix(vector<vector<double>>& mat);
vector<double> get_diag(const vector<vector<double>>& mat);

void print_vector(const vector<double>& vec);
void print_matrix(const vector<vector<double>>& mat);
void print_matrix_transpose(const vector<vector<double>>& mat);

#endif
