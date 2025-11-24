#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "gs_QR.h"
#include "householder_QR.h"
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

    // Make a copy for the second run
    vector<vector<float>> matrix_copy = in_matrix;

    cout << "Running Gram-Schmidt QR..." << endl;
    gs_QR(in_matrix, n);
    cout << " " << endl;

    cout << "Running Householder QR..." << endl;
    householder_QR(matrix_copy, n);

    return 0;
}
