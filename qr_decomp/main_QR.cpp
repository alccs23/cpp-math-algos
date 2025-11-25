#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "QR_decomp.h"
#include "eigen_solve.h"
#include "../utils/matrix_utils.h" 

using namespace std;

int main() {

    string filename = "matrix.txt";
    ifstream infile(filename);

    if (!infile.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return 1;
    }

    int n;
    infile >> n;

    vector<vector<double>> in_matrix(n, vector<double>(n));

    // read matrix values
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            infile >> in_matrix[i][j];
        }
    }

    infile.close();

    auto [eigenvects, eigenvals] = qr_algorithm(in_matrix);
    cout << "Av: " << endl;
    print_vector(matrix_vector_mult(in_matrix, get_column(eigenvects, 0)));
    cout << endl;
    print_vector(scalar_multiply(get_column(eigenvects, 0), eigenvals[0]));

}
