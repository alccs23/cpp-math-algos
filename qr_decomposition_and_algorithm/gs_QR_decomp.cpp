#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// Since I only care about real life stuff, we assume all square matrix systems

float dot_product(const vector<float>& vec1, const vector<float>& vec2)
{
    float dot_product = 0.0;
    int size = vec1.size();

    for (int i = 0; i < size; i++){
        dot_product += vec1[i] * vec2[i];
    }

    return dot_product;
}

vector<float> subtract_vectors(const vector<float>& vec1, const vector<float>& vec2) {
    size_t size = vec1.size();
    vector<float> result(size);
    for (size_t i = 0; i < size; i++) {
        result[i] = vec1[i] - vec2[i];
    }
    return result;
}

vector<float> scalar_multiply(const vector<float>& vec, float scalar) {
    size_t size = vec.size();
    vector<float> result(size);
    for (size_t i = 0; i < size; i++) {
        result[i] = vec[i] * scalar;
    }
    return result;
}

void print_vector(const vector<float>& vec) {
    for (float val : vec) {
        cout << val << " ";
    }
    cout << endl;
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

void print_matrix(const vector<vector<float>>& mat) {
    for (const auto& row : mat) {
        for (float val : row) {
            cout << val << " ";
        }
        cout << endl;
    }
}

int main() {
    int n;

    cout << "Enter the size of the square matrix: ";
    cin >> n;

    vector<vector<float>> in_matrix(n, vector<float>(n));

    cout << "Enter the matrix values row by row:\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> in_matrix[j][i]; 
        }
    }

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
    return 0;
}
