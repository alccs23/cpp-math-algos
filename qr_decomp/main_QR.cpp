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

    vector<vector<float>> in_matrix(n, vector<float>(n));

    // read matrix values
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            infile >> in_matrix[i][j];
        }
    }

    infile.close();

    print_vector(qr_algorithm(in_matrix));
    
}
