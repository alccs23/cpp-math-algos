#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "QR_decomp.h"
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

    vector<vector<float>> in_matrix(n, vector<float>(n));

    // read matrix values
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            infile >> in_matrix[i][j];
        }
    }

    infile.close();

    auto [Q, R] = gs_QR(in_matrix, n);
    print_matrix(Q);
    print_matrix(R);
    print_matrix(matrix_mult(Q, R));
    return 0;
}
