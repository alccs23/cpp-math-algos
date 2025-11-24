#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

using namespace std;

// Since I only care about real life stuff, we assume all square matrix systems

int sgn(float val) {
    return (0.0f < val) - (val < 0.0f);
}

vector<vector<float>> scalar_multiply_matrix(const vector<vector<float>>& A, float scalar) {
    int rows = A.size();
    int cols = A[0].size();

    vector<vector<float>> result(rows, vector<float>(cols));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[i][j] = A[i][j] * scalar;
        }
    }

    return result;
}


vector<vector<float>> subtract_matrix(const vector<vector<float>>& A, const vector<vector<float>>& B) {
    int rows = A.size();
    int cols = A[0].size();

    vector<vector<float>> result(rows, vector<float>(cols));

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[i][j] = A[i][j] - B[i][j];
        }
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

float dot_product(const vector<float>& vec1, const vector<float>& vec2)
{
    float dot_product = 0.0;
    int size = vec1.size();

    for (int i = 0; i < size; i++){
        dot_product += vec1[i] * vec2[i];
    }

    return dot_product;
}

vector<float> add_vectors(const vector<float>& vec1, const vector<float>& vec2) {
    size_t size = vec1.size();
    vector<float> result(size);
    for (size_t i = 0; i < size; i++) {
        result[i] = vec1[i] + vec2[i];
    }
    return result;
}

vector<float> get_column(const vector<vector<float>>& mat, int col){
    vector<float> result(mat.size());
    for (int i = 0; i < mat.size(); i++){
        result[i] = mat[i][col];
    }
    
    return result;
}

void set_column(vector<vector<float>>& mat, const vector<float>& vect, int col){
    vector<float> result(mat.size());
    for (int i = 0; i < mat.size(); i++){
        mat[i][col] = vect[i];
    }
}

vector<vector<float>> outer_product(const vector<float>& vect1, const vector<float>& vect2){
    int n = vect1.size();

    vector<vector<float>> result(n, vector<float>(n));

    for (int i = 0; i < n; i ++){
        for (int j = 0; j < n; j++){
            result[i][j] = vect1[i]*vect2[j];
        }
    }
    return result;
}

vector<vector<float>> identity_matrix(int n) {
    vector<vector<float>> I(n, vector<float>(n, 0));

    for (int i = 0; i < n; i++) {
        I[i][i] = 1;
    }

    return I;
}

vector<vector<float>> matrix_mult(const vector<vector<float>>& mat1, const vector<vector<float>>& mat2) {
    // we assume that the matrices are made up of column vectors
    int n = mat1.size();

    vector<vector<float>> result(n, vector<float>(n, 0.0f));
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            float sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += mat1[i][k] * mat2[k][j];
            }
            result[i][j] = sum;
        }
    }
    return result;
}

void print_vector(const vector<float>& vec) {
    for (float val : vec) {
        cout << val << " ";
    }
    cout << endl;
}

void print_matrix(const vector<vector<float>>& mat) {
    const float eps = 1e-6f;
    for (const auto& row : mat) {
        for (float val : row) {
            if (abs(val) < eps) val = 0.0f;  // remove tiny artifacts
            cout << val << " ";
        }
        cout << endl;
    }
}

int main() {
    string filename = "matrix.txt";
    ifstream infile(filename);

    if (!infile.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return 1;
    }

    int n;
    infile >> n;

    vector<vector<float>> in_matrix(n, vector<float>(n));

    // read matrix values
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            infile >> in_matrix[i][j];
        }
    }

    infile.close();

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
    
    return 0;
}